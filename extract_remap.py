#!/usr/bin/env python3
"""Excise a reference to the regions covered by one or more samples, then remap each sample.

Input: a reference FASTA and a sample manifest (TSV, header required):

    name    r1              r2              bam
    s1      s1_R1.fq.gz     s1_R2.fq.gz     s1.bam
    s2      s2.fq.gz                        -

  * name: letters, digits, '.', '_' or '-'; unique.
  * r1: required. r2: optional (paired if set, single-end if empty or '-').
  * bam: optional (empty or '-' = none). Coordinate-sorted, indexed, aligned to --ref.

Pipeline:
  1. Coverage BAM per sample (given, or first-pass mapping to the full reference).
  2. Originally mapped reads per sample (primary, mapped records).
  3. Kept contigs = contigs with a primary mapped record in any sample.
  4. Covered positions per sample (depth >= --min-depth, kept contigs only), unioned,
     padded by --flank, clamped and merged -> subset reference (one for all samples).
  5. Each sample's originally mapped reads are remapped to the subset reference.

Layers (kept separate so the mapper can be swapped later):
  * Utilities      - FASTA/FASTQ/BAM helpers, samtools calls, interval logic. No mapper code.
  * Mapper         - abstract interface (index / align). Concrete: BwaMemMapper.
  * Pipeline       - main(): wires utilities and a mapper together.
"""
import argparse
import csv
import gzip
import os
import re
import subprocess
import sys
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Dict, Iterable, Iterator, List, Optional, Set, Tuple

PRIMARY_MAPPED_EXCLUDE = "0x904"  # unmapped | secondary | supplementary
NAME_RE = re.compile(r"[A-Za-z0-9._-]+")


# ---------------------------------------------------------------- utilities

def log(msg: str) -> None:
    print(f"[extract_remap] {msg}", file=sys.stderr)


def run(cmd: List[str]) -> None:
    log("run: " + " ".join(cmd))
    subprocess.run(cmd, check=True)


def open_text(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


def read_fasta(path: str) -> Dict[str, str]:
    seqs: Dict[str, str] = {}
    name: Optional[str] = None
    chunks: List[str] = []
    with open_text(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if name is not None:
        seqs[name] = "".join(chunks)
    if not seqs:
        raise ValueError(f"no sequences in {path}")
    return seqs


def write_fasta(path: str, records: Iterable[Tuple[str, str]], width: int = 60) -> None:
    with open(path, "w") as fh:
        for name, seq in records:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i:i + width] + "\n")


def fasta_lengths(path: str) -> Dict[str, int]:
    """Contig lengths from one streaming pass; sequence is not kept in memory."""
    lengths: Dict[str, int] = {}
    name: Optional[str] = None
    n = 0
    with open_text(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    lengths[name] = n
                name = line[1:].split()[0]
                n = 0
            else:
                n += len(line)
    if name is not None:
        lengths[name] = n
    if not lengths:
        raise ValueError(f"no sequences in {path}")
    return lengths


def _wrap_write(fh, seq: str, col: int, width: int) -> int:
    """Write seq wrapped at `width`, continuing from column `col`. Returns new column."""
    i = 0
    while i < len(seq):
        piece = seq[i:i + (width - col)]
        fh.write(piece)
        i += len(piece)
        col += len(piece)
        if col == width:
            fh.write("\n")
            col = 0
    return col


def group_by_contig(intervals: List[Tuple[str, int, int]]) -> List[Tuple[str, List[Tuple[int, int]]]]:
    """[(contig, s, e)] in reference order -> [(contig, [(s, e), ...])], one entry per contig."""
    groups: List[Tuple[str, List[Tuple[int, int]]]] = []
    for c, s, e in intervals:
        if groups and groups[-1][0] == c:
            groups[-1][1].append((s, e))
        else:
            groups.append((c, [(s, e)]))
    return groups


def spacer_sequence(spacer: int, first: bool) -> str:
    """N run placed between two joined intervals; empty for the first interval or spacer 0."""
    return "" if first or spacer <= 0 else "N" * spacer


def extract_intervals(ref_path: str, intervals: List[Tuple[str, int, int]], out_fa: str,
                      width: int = 60, spacer: int = 0) -> None:
    """Stream the reference once and write one record per contig.

    Each contig's covered intervals (1-based, inclusive) are concatenated, in order,
    into a single record named after the original contig. With spacer > 0, a run of
    `spacer` N bases is inserted between consecutive intervals of the same contig.
    `intervals` must be in reference order, which pad_and_merge guarantees.
    """
    idx = 0
    n = len(intervals)
    cur: Optional[str] = None
    pos = 0          # 1-based position of the last base consumed on `cur`
    col = 0          # output column for wrapping
    active = False   # current contig has pending intervals
    begun = -1       # index of the interval whose start (and spacer) was already written
    emitted = False  # current contig already has at least one interval written
    with open_text(ref_path) as fin, open(out_fa, "w") as fout:
        for raw in fin:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if active and col:
                    fout.write("\n")
                cur = line[1:].split()[0]
                pos = 0
                col = 0
                emitted = False
                active = idx < n and intervals[idx][0] == cur
                if active:
                    fout.write(f">{cur}\n")
                continue
            if not active:
                continue
            L = len(line)
            while idx < n and intervals[idx][0] == cur:
                _, s, e = intervals[idx]
                if s > pos + L:
                    break  # starts on a later line
                if begun != idx:
                    gap = spacer_sequence(spacer, not emitted)
                    if gap:
                        col = _wrap_write(fout, gap, col, width)
                    begun = idx
                    emitted = True
                lo = max(s, pos + 1)
                hi = min(e, pos + L)
                col = _wrap_write(fout, line[lo - pos - 1:hi - pos], col, width)
                if e <= pos + L:  # interval ends on this line
                    idx += 1
                else:
                    break  # continues on the next line
            pos += L
            if idx < n and intervals[idx][0] != cur:  # this contig is finished
                if col:
                    fout.write("\n")
                col = 0
                active = False
        if active and col:
            fout.write("\n")


def read_name(header: str) -> str:
    """Normalise a FASTQ header or BAM QNAME to a shared read name."""
    name = header.lstrip("@>").split()[0]
    if name.endswith("/1") or name.endswith("/2"):
        name = name[:-2]
    return name


def filter_fastq(src: str, dst: str, keep: Set[str]) -> int:
    """Write records whose name is in `keep`. Pairs share names, so both mates are kept."""
    kept = 0
    with open_text(src) as fin, gzip.open(dst, "wt") as fout:
        while True:
            rec = [fin.readline() for _ in range(4)]
            if not rec[0]:
                break
            if any(not line for line in rec):
                raise ValueError(f"truncated FASTQ record in {src}")
            if read_name(rec[0]) in keep:
                fout.writelines(rec)
                kept += 1
    return kept


def samtools_faidx(fasta: str) -> None:
    run(["samtools", "faidx", fasta])


def samtools_index(bam: str) -> None:
    run(["samtools", "index", bam])


def bam_primary_mapped(bam: str) -> Tuple[Set[str], Set[str]]:
    """Return (read names, contigs) over primary mapped records of a BAM."""
    proc = subprocess.Popen(
        ["samtools", "view", "-F", PRIMARY_MAPPED_EXCLUDE, bam],
        stdout=subprocess.PIPE, text=True,
    )
    names: Set[str] = set()
    contigs: Set[str] = set()
    assert proc.stdout is not None
    for line in proc.stdout:
        fields = line.split("\t", 3)
        names.add(read_name(fields[0]))
        contigs.add(fields[2])
    if proc.wait() != 0:
        raise RuntimeError("samtools view failed")
    return names, contigs


def bam_depth_rows(bam: str, keep: Optional[Set[str]] = None) -> Iterator[Tuple[str, int, int]]:
    """Yield (contig, 1-based pos, depth) for positions with depth > 0, streamed.

    If `keep` is given, only contigs in that set are yielded.
    """
    proc = subprocess.Popen(["samtools", "depth", bam], stdout=subprocess.PIPE, text=True)
    assert proc.stdout is not None
    for line in proc.stdout:
        contig, pos, depth = line.rstrip("\n").split("\t")
        if keep is not None and contig not in keep:
            continue
        yield contig, int(pos), int(depth)
    if proc.wait() != 0:
        raise RuntimeError("samtools depth failed")


def covered_runs(rows: Iterable[Tuple[str, int, int]], min_depth: int) -> List[Tuple[str, int, int]]:
    """Collapse per-base depth into contiguous (contig, start, end) runs, 1-based inclusive."""
    runs: List[List] = []
    cur: Optional[List] = None
    for contig, pos, depth in rows:
        if depth < min_depth:
            continue
        if cur is not None and cur[0] == contig and pos == cur[2] + 1:
            cur[2] = pos
        else:
            if cur is not None:
                runs.append(cur)
            cur = [contig, pos, pos]
    if cur is not None:
        runs.append(cur)
    return [(c, s, e) for c, s, e in runs]


def pad_and_merge(runs: List[Tuple[str, int, int]], lengths: Dict[str, int],
                  flank: int) -> List[Tuple[str, int, int]]:
    """Pad runs by `flank`, clamp to contig bounds, merge overlapping or adjacent intervals.

    Runs from several samples can be passed together; the result is their union.
    Only contigs present in `lengths` are emitted, in that order.
    """
    by_contig: Dict[str, List[Tuple[int, int]]] = {}
    for contig, s, e in runs:
        by_contig.setdefault(contig, []).append((max(1, s - flank), min(lengths[contig], e + flank)))
    out: List[Tuple[str, int, int]] = []
    for contig in lengths:  # keep reference order
        merged: List[List[int]] = []
        for s, e in sorted(by_contig.get(contig, [])):
            if merged and s <= merged[-1][1] + 1:
                merged[-1][1] = max(merged[-1][1], e)
            else:
                merged.append([s, e])
        out.extend((contig, s, e) for s, e in merged)
    return out


def subset_sequence(seq: str, ivs: List[Tuple[int, int]], spacer: int) -> str:
    """Concatenate 1-based inclusive intervals of `seq`, with N spacers between them."""
    parts: List[str] = []
    for i, (s, e) in enumerate(ivs):
        parts.append(spacer_sequence(spacer, i == 0))
        parts.append(seq[s - 1:e])
    return "".join(parts)


def write_bed(path: str, intervals: List[Tuple[str, int, int]], spacer: int = 0) -> None:
    """Coordinate map from subset to original reference, one line per interval.

    Columns (0-based, half-open): subset contig, subset start, subset end,
    original contig, original start, original end. Subset contigs are the original
    contig names; each contig's intervals are concatenated in the order given here,
    with `spacer` N bases between consecutive intervals (not counted in any row).
    """
    with open(path, "w") as fh:
        fh.write("#subset_contig\tsubset_start\tsubset_end\torig_contig\torig_start\torig_end\n")
        for contig, ivs in group_by_contig(intervals):
            offset = 0
            for i, (s, e) in enumerate(ivs):
                if i > 0:
                    offset += max(spacer, 0)
                length = e - s + 1
                fh.write(f"{contig}\t{offset}\t{offset + length}\t{contig}\t{s - 1}\t{e}\n")
                offset += length


def write_empty_bam(out_bam: str, fai: str, read_group: str) -> None:
    """Write a header-only, coordinate-sorted BAM with @SQ lines from a .fai file."""
    lines = ["@HD\tVN:1.6\tSO:coordinate"]
    with open(fai) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            lines.append(f"@SQ\tSN:{fields[0]}\tLN:{fields[1]}")
    lines.append(read_group.replace("\\t", "\t"))
    log(f"writing empty BAM: {out_bam}")
    subprocess.run(["samtools", "view", "-b", "-o", out_bam, "-"],
                   input="\n".join(lines) + "\n", text=True, check=True)


# --------------------------------------------------------------- manifest

@dataclass
class Sample:
    name: str
    r1: str
    r2: Optional[str]
    bam: Optional[str]

    @property
    def reads(self) -> List[str]:
        return [self.r1, self.r2] if self.r2 else [self.r1]

    @property
    def paired(self) -> bool:
        return self.r2 is not None


def _cell(value: Optional[str]) -> Optional[str]:
    value = (value or "").strip()
    return None if value in ("", "-") else value


def load_manifest(path: str) -> List[Sample]:
    samples: List[Sample] = []
    seen: Set[str] = set()
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames or "name" not in reader.fieldnames or "r1" not in reader.fieldnames:
            raise ValueError("manifest needs a tab-separated header with at least: name, r1 (optional: r2, bam)")
        for lineno, row in enumerate(reader, start=2):
            name = (row.get("name") or "").strip()
            if not NAME_RE.fullmatch(name):
                raise ValueError(f"manifest line {lineno}: invalid sample name {name!r}")
            if name in seen:
                raise ValueError(f"manifest line {lineno}: duplicate sample name {name!r}")
            seen.add(name)
            r1 = _cell(row.get("r1"))
            r2 = _cell(row.get("r2"))
            bam = _cell(row.get("bam"))
            if r1 is None:
                raise ValueError(f"manifest line {lineno} ({name}): r1 is required")
            for label, p in (("r1", r1), ("r2", r2), ("bam", bam)):
                if p is not None and not os.path.isfile(p):
                    raise ValueError(f"manifest line {lineno} ({name}): {label} not found: {p}")
            samples.append(Sample(name, r1, r2, bam))
    if not samples:
        raise ValueError("manifest has no samples")
    return samples


def read_group_arg(sample: str) -> str:
    """bwa -R argument (bwa expands the literal \\t escapes)."""
    return f"@RG\\tID:{sample}\\tSM:{sample}\\tPL:ILLUMINA"


# ------------------------------------------------------------ mapper layer

class Mapper(ABC):
    """Interface any aligner must implement. Pipeline only talks to this."""

    name: str = "abstract"

    @abstractmethod
    def index(self, ref_fasta: str) -> None:
        """Build the aligner index for ref_fasta."""

    @abstractmethod
    def align(self, ref_fasta: str, reads: List[str], out_bam: str,
              read_group: str, threads: int) -> None:
        """Align reads (1 file = single-end, 2 files = paired) and write a sorted BAM."""


class BwaMemMapper(Mapper):
    name = "bwa-mem"

    def index(self, ref_fasta: str) -> None:
        run(["bwa", "index", ref_fasta])

    def align(self, ref_fasta: str, reads: List[str], out_bam: str,
              read_group: str, threads: int) -> None:
        mem = subprocess.Popen(
            ["bwa", "mem", "-t", str(threads), "-R", read_group, ref_fasta, *reads],
            stdout=subprocess.PIPE,
        )
        log("run: bwa mem ... " + " ".join(reads))
        sort = subprocess.Popen(
            ["samtools", "sort", "-@", str(threads), "-o", out_bam, "-"],
            stdin=mem.stdout,
        )
        assert mem.stdout is not None
        mem.stdout.close()
        sort_rc = sort.wait()
        mem_rc = mem.wait()
        if mem_rc != 0 or sort_rc != 0:
            raise RuntimeError(f"alignment failed (bwa rc={mem_rc}, sort rc={sort_rc})")


MAPPERS: Dict[str, type] = {BwaMemMapper.name: BwaMemMapper}
DEFAULT_MAPPER = BwaMemMapper.name


# ---------------------------------------------------------------- pipeline

def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--ref", required=True, help="full reference FASTA")
    p.add_argument("--manifest", required=True, help="sample manifest TSV (see above)")
    p.add_argument("--outdir", required=True, help="output directory")
    p.add_argument("--min-depth", type=int, default=1,
                   help="minimum depth per sample to count as covered (default 1)")
    p.add_argument("--flank", type=int, default=20, help="bp padding around covered regions (default 20)")
    p.add_argument("--spacer", type=int, default=0,
                   help="N bases inserted between excised intervals of the same contig "
                        "(default 0 = joined directly)")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--stream-threshold", type=float, default=500,
                   help="reference size in MB above which the FASTA is streamed instead of "
                        "loaded into memory (default 500; 0 = always stream)")
    args = p.parse_args(argv)

    if args.min_depth < 1:
        p.error("--min-depth must be >= 1")
    if args.flank < 0:
        p.error("--flank must be >= 0")
    if args.spacer < 0:
        p.error("--spacer must be >= 0")
    return args


def main(argv: Optional[List[str]] = None) -> int:
    args = parse_args(argv)
    samples = load_manifest(args.manifest)
    os.makedirs(args.outdir, exist_ok=True)
    mapper: Mapper = MAPPERS[DEFAULT_MAPPER]()
    log(f"{len(samples)} sample(s) in manifest")

    ref_mb = os.path.getsize(args.ref) / (1024 * 1024)
    streaming = ref_mb > args.stream_threshold
    if streaming:
        lengths = fasta_lengths(args.ref)
        ref = None
        log(f"reference {ref_mb:.1f} MB > threshold: streaming mode")
    else:
        ref = read_fasta(args.ref)
        lengths = {name: len(seq) for name, seq in ref.items()}
        log(f"reference {ref_mb:.1f} MB <= threshold: in-memory mode")

    # 1. coverage BAM per sample: given, or first-pass mapping to the full reference
    if any(s.bam is None for s in samples):
        log("samples without bam: mapping them to the full reference first")
        mapper.index(args.ref)
    sample_bams: Dict[str, str] = {}
    for s in samples:
        sdir = os.path.join(args.outdir, s.name)
        os.makedirs(sdir, exist_ok=True)
        if s.bam is not None:
            sample_bams[s.name] = s.bam
            log(f"{s.name}: coverage source {s.bam}")
        else:
            bam = os.path.join(sdir, "first_pass.bam")
            mapper.align(args.ref, s.reads, bam, read_group_arg(s.name), args.threads)
            samtools_index(bam)
            sample_bams[s.name] = bam
            log(f"{s.name}: first-pass mapping -> {bam}")

    # 2. originally mapped reads and contigs per sample
    mapped: Dict[str, Set[str]] = {}
    kept: Set[str] = set()
    for s in samples:
        names, contigs = bam_primary_mapped(sample_bams[s.name])
        mapped[s.name] = names
        kept |= contigs
        log(f"{s.name}: {len(names)} originally mapped reads on {len(contigs)} contigs")
    if not kept:
        raise SystemExit("no sample has primary mapped reads; nothing to do")
    unknown = sorted(c for c in kept if c not in lengths)
    if unknown:
        raise SystemExit(f"BAM contigs not found in reference: {', '.join(unknown[:5])}")
    kept_order = [c for c in lengths if c in kept]
    contigs_file = os.path.join(args.outdir, "contigs_kept.txt")
    with open(contigs_file, "w") as fh:
        fh.write("".join(f"{c}\n" for c in kept_order))
    log(f"contigs kept (mapped by at least one sample): {len(kept_order)} -> {contigs_file}")

    # 3. covered positions per sample (min-depth per sample), unioned on kept contigs
    runs: List[Tuple[str, int, int]] = []
    for s in samples:
        sample_runs = covered_runs(bam_depth_rows(sample_bams[s.name], kept), args.min_depth)
        log(f"{s.name}: {len(sample_runs)} covered runs on kept contigs")
        runs.extend(sample_runs)
    intervals = pad_and_merge(runs, {c: lengths[c] for c in kept_order}, args.flank)
    if not intervals:
        raise SystemExit("no covered positions; nothing to extract")
    total_bp = sum(e - s + 1 for _, s, e in intervals)
    log(f"union intervals (padded, merged): {len(intervals)}, {total_bp} bp")

    # 4. subset reference, one for all samples
    bed = os.path.join(args.outdir, "coverage_intervals.bed")
    write_bed(bed, intervals, args.spacer)
    subset_fa = os.path.join(args.outdir, "subset_ref.fa")
    if streaming:
        extract_intervals(args.ref, intervals, subset_fa, spacer=args.spacer)
    else:
        write_fasta(subset_fa, ((c, subset_sequence(ref[c], ivs, args.spacer))
                                for c, ivs in group_by_contig(intervals)))
    samtools_faidx(subset_fa)
    mapper.index(subset_fa)

    # 5. remap each sample's originally mapped reads to the subset reference
    for s in samples:
        sdir = os.path.join(args.outdir, s.name)
        if s.paired:
            out_reads = [os.path.join(sdir, f"reads_subset_R{i}.fastq.gz") for i in (1, 2)]
        else:
            out_reads = [os.path.join(sdir, "reads_subset.fastq.gz")]
        n_kept = 0
        for src, dst in zip(s.reads, out_reads):
            n = filter_fastq(src, dst, mapped[s.name])
            log(f"{s.name}: filtered {src} -> {dst}: {n} records")
            n_kept = n
        out_bam = os.path.join(sdir, "remapped.sorted.bam")
        if n_kept == 0:
            log(f"{s.name}: no reads to remap; writing empty BAM")
            write_empty_bam(out_bam, subset_fa + ".fai", read_group_arg(s.name))
        else:
            mapper.align(subset_fa, out_reads, out_bam, read_group_arg(s.name), args.threads)
        samtools_index(out_bam)

    log(f"done: {args.outdir} ({len(samples)} sample(s), subset {subset_fa})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
