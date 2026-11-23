"""Tests for extract_remap.

Run from the repository root:

    python -m unittest discover -s tests -v

End-to-end tests need `samtools` and `bwa` on PATH and are skipped otherwise.
"""
import csv
import gzip
import os
import random
import shutil
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

import extract_remap as xr  # noqa: E402

HAVE_TOOLS = shutil.which("samtools") is not None and shutil.which("bwa") is not None


def read_fasta_dict(path):
    out, name = {}, None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                name = line[1:].split()[0]
                out[name] = []
            else:
                out[name].append(line)
    return {k: "".join(v) for k, v in out.items()}


def read_bed(path):
    """Parse coverage_intervals.bed into rows of (subset_contig, subset_start, subset_end,
    orig_contig, orig_start, orig_end); numeric columns are ints."""
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            sc, ss, se, oc, os_, oe = line.split()
            rows.append((sc, int(ss), int(se), oc, int(os_), int(oe)))
    return rows


def rand_seq(rng, n):
    return "".join(rng.choice("ACGT") for _ in range(n))


class TempDirTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="extract_remap_test_")

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def path(self, name):
        return os.path.join(self.tmp, name)


class TestCoverageLogic(unittest.TestCase):
    def test_covered_runs_merges_consecutive_and_respects_min_depth(self):
        rows = [("c1", 1, 3), ("c1", 2, 3), ("c1", 3, 1), ("c1", 4, 2), ("c1", 6, 5), ("c2", 1, 9)]
        self.assertEqual(xr.covered_runs(rows, 2), [("c1", 1, 2), ("c1", 4, 4), ("c1", 6, 6), ("c2", 1, 1)])
        self.assertEqual(xr.covered_runs(rows, 1)[0], ("c1", 1, 4))

    def test_covered_runs_does_not_join_across_contigs(self):
        rows = [("c1", 5, 1), ("c2", 6, 1)]
        self.assertEqual(xr.covered_runs(rows, 1), [("c1", 5, 5), ("c2", 6, 6)])

    def test_pad_and_merge_pads_clamps_and_merges(self):
        lengths = {"c1": 100, "c2": 50}
        runs = [("c1", 5, 10), ("c1", 30, 40), ("c2", 45, 50)]
        out = xr.pad_and_merge(runs, lengths, flank=10)
        # c1: [1,20], [20,50] -> merged [1,50]; c2: [35,50] clamped
        self.assertEqual(out, [("c1", 1, 50), ("c2", 35, 50)])

    def test_pad_and_merge_unions_samples_and_keeps_reference_order(self):
        lengths = {"c1": 1000, "c2": 1000}
        s1 = [("c2", 100, 200)]
        s2 = [("c1", 500, 600), ("c2", 150, 300)]
        out = xr.pad_and_merge(s1 + s2, lengths, flank=0)
        self.assertEqual(out, [("c1", 500, 600), ("c2", 100, 300)])

    def test_pad_and_merge_adjacent_intervals_join(self):
        out = xr.pad_and_merge([("c1", 1, 10), ("c1", 11, 20)], {"c1": 50}, flank=0)
        self.assertEqual(out, [("c1", 1, 20)])

    def test_group_by_contig(self):
        ivs = [("c1", 1, 5), ("c1", 9, 12), ("c2", 3, 4)]
        self.assertEqual(xr.group_by_contig(ivs),
                         [("c1", [(1, 5), (9, 12)]), ("c2", [(3, 4)])])


class TestSubsetSequence(TempDirTest):
    def test_spacer_zero_joins_directly(self):
        self.assertEqual(xr.subset_sequence("ACGTACGT", [(1, 2), (5, 6)], 0), "ACAC")

    def test_spacer_inserts_n_between_not_at_ends(self):
        self.assertEqual(xr.subset_sequence("ACGTACGT", [(1, 2), (5, 6), (7, 8)], 3),
                         "AC" + "NNN" + "AC" + "NNN" + "GT")

    def test_extract_intervals_matches_in_memory_for_random_layouts(self):
        rng = random.Random(7)
        seqs = {"a": rand_seq(rng, 523) + "acgtN", "b": rand_seq(rng, 40), "c": "G", "d": rand_seq(rng, 301)}
        ref = self.path("messy.fa")
        with open(ref, "w") as fh:
            for k, (name, s) in enumerate(seqs.items()):
                width = [7, 60, 1000][k % 3]
                fh.write(f">{name}\n")
                for i in range(0, len(s), width):
                    fh.write(s[i:i + width] + "\n")
            # no trailing newline on the last record
        with open(ref) as fh:
            data = fh.read().rstrip("\n")
        with open(ref, "w") as fh:
            fh.write(data)

        order = list(seqs)
        for trial in range(200):
            ivs = []
            for c in order:
                L = len(seqs[c])
                for _ in range(rng.randint(0, 4)):
                    s = rng.randint(1, L)
                    e = rng.randint(s, L)
                    ivs.append((c, s, e))
            ivs.sort(key=lambda t: (order.index(t[0]), t[1]))
            # pad_and_merge produces merged, sorted intervals; mirror that
            merged = []
            for c, s, e in ivs:
                if merged and merged[-1][0] == c and s <= merged[-1][2] + 1:
                    merged[-1] = (c, merged[-1][1], max(merged[-1][2], e))
                else:
                    merged.append((c, s, e))
            if not merged:
                continue
            spacer = [0, 1, 5][trial % 3]
            out = self.path("out.fa")
            xr.extract_intervals(ref, merged, out, width=60, spacer=spacer)
            expected = self.path("exp.fa")
            with open(expected, "w") as fh:
                for c, ivs_c in xr.group_by_contig(merged):
                    seq = xr.subset_sequence(seqs[c], ivs_c, spacer)
                    fh.write(f">{c}\n")
                    for i in range(0, len(seq), 60):
                        fh.write(seq[i:i + 60] + "\n")
            with open(out) as a, open(expected) as b:
                self.assertEqual(a.read(), b.read(), f"trial {trial} spacer {spacer} {merged}")

    def test_extract_intervals_spacer_positions_on_wrapped_lines(self):
        # spacer that straddles the 60-column wrap must still be correct
        seq = "A" * 58 + "C" * 10
        ref = self.path("ref.fa")
        with open(ref, "w") as fh:
            fh.write(">c1\n" + seq[:60] + "\n" + seq[60:] + "\n")
        out = self.path("out.fa")
        xr.extract_intervals(ref, [("c1", 1, 3), ("c1", 60, 62)], out, width=60, spacer=4)
        got = read_fasta_dict(out)["c1"]
        self.assertEqual(got, "AAA" + "NNNN" + seq[59:62])


class TestBed(TempDirTest):
    def test_bed_offsets_without_spacer(self):
        bed = self.path("b.bed")
        xr.write_bed(bed, [("c1", 10, 20), ("c1", 30, 40), ("c2", 1, 5)], spacer=0)
        rows = [list(map(str, r)) for r in read_bed(bed)]
        self.assertEqual(rows, [
            ["c1", "0", "11", "c1", "9", "20"],
            ["c1", "11", "22", "c1", "29", "40"],
            ["c2", "0", "5", "c2", "0", "5"],
        ])

    def test_bed_offsets_with_spacer_skip_gap(self):
        bed = self.path("b.bed")
        xr.write_bed(bed, [("c1", 10, 20), ("c1", 30, 40)], spacer=5)
        rows = read_bed(bed)
        self.assertEqual(list(rows[1][1:3]), [16, 27])

    def test_bed_maps_subset_back_to_reference(self):
        rng = random.Random(1)
        ref_seq = rand_seq(rng, 500)
        ivs = [("c1", 20, 60), ("c1", 200, 230), ("c1", 400, 401)]
        for spacer in (0, 9):
            subset = xr.subset_sequence(ref_seq, [(s, e) for _, s, e in ivs], spacer)
            bed = self.path(f"b{spacer}.bed")
            xr.write_bed(bed, ivs, spacer)
            for sc, ss, se, oc, os_, oe in read_bed(bed):
                self.assertEqual(subset[ss:se], ref_seq[os_:oe])


class TestManifest(TempDirTest):
    def write_manifest(self, rows, header=("name", "r1", "r2", "bam")):
        path = self.path("m.tsv")
        with open(path, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(header)
            w.writerows(rows)
        return path

    def touch(self, name):
        p = self.path(name)
        open(p, "w").close()
        return p

    def test_valid_manifest_with_dash_and_empty_cells(self):
        r1 = self.touch("s1_R1.fq")
        r2 = self.touch("s1_R2.fq")
        r3 = self.touch("s2.fq")
        m = self.write_manifest([["s1", r1, r2, "-"], ["s2", r3, "", ""]])
        samples = xr.load_manifest(m)
        self.assertEqual([s.name for s in samples], ["s1", "s2"])
        self.assertTrue(samples[0].paired)
        self.assertIsNone(samples[0].bam)
        self.assertFalse(samples[1].paired)
        self.assertEqual(samples[1].reads, [r3])

    def test_duplicate_name_rejected(self):
        r1 = self.touch("a.fq")
        m = self.write_manifest([["s1", r1, "", ""], ["s1", r1, "", ""]])
        with self.assertRaisesRegex(ValueError, "duplicate"):
            xr.load_manifest(m)

    def test_invalid_name_rejected(self):
        r1 = self.touch("a.fq")
        m = self.write_manifest([["bad name/x", r1, "", ""]])
        with self.assertRaisesRegex(ValueError, "invalid sample name"):
            xr.load_manifest(m)

    def test_missing_file_rejected(self):
        m = self.write_manifest([["s1", self.path("nope.fq"), "", ""]])
        with self.assertRaisesRegex(ValueError, "not found"):
            xr.load_manifest(m)

    def test_missing_r1_rejected(self):
        m = self.write_manifest([["s1", "", "", ""]])
        with self.assertRaisesRegex(ValueError, "r1 is required"):
            xr.load_manifest(m)

    def test_missing_header_rejected(self):
        path = self.path("bad.tsv")
        with open(path, "w") as fh:
            fh.write("sample\tfile\n")
        with self.assertRaisesRegex(ValueError, "header"):
            xr.load_manifest(path)


class TestArgs(unittest.TestCase):
    def base(self, *extra):
        return ["--ref", "r.fa", "--manifest", "m.tsv", "--outdir", "o", *extra]

    def test_spacer_default_zero(self):
        self.assertEqual(xr.parse_args(self.base()).spacer, 0)

    def test_spacer_accepted(self):
        self.assertEqual(xr.parse_args(self.base("--spacer", "100")).spacer, 100)

    def test_negative_spacer_rejected(self):
        with self.assertRaises(SystemExit):
            xr.parse_args(self.base("--spacer", "-1"))


class TestFastaIO(TempDirTest):
    def test_fasta_lengths_matches_read_fasta(self):
        path = self.path("r.fa")
        with open(path, "w") as fh:
            fh.write(">a desc\nACGT\nAC\n>b\nGGGGG")
        self.assertEqual(xr.fasta_lengths(path), {"a": 6, "b": 5})
        self.assertEqual(read_fasta_dict(path), {"a": "ACGTAC", "b": "GGGGG"})

    def test_read_name_normalises_suffixes(self):
        self.assertEqual(xr.read_name("@r1/1 extra"), "r1")
        self.assertEqual(xr.read_name("r1/2"), "r1")
        self.assertEqual(xr.read_name("@r1"), "r1")


@unittest.skipUnless(HAVE_TOOLS, "samtools and bwa required")
class TestEndToEnd(TempDirTest):
    """Full pipeline on synthetic data: streaming vs in-memory, with and without spacer."""

    @classmethod
    def setUpClass(cls):
        cls._tmp_root = tempfile.mkdtemp(prefix="extract_remap_e2e_")
        d = cls._tmp_root
        rng = random.Random(42)
        c1 = rand_seq(rng, 6000)
        c2 = rand_seq(rng, 3000)
        c3 = rand_seq(rng, 2000)  # receives no reads: must be dropped
        ref = os.path.join(d, "ref.fa")
        with open(ref, "w") as fh:
            for name, s in (("c1", c1), ("c2", c2), ("c3", c3)):
                fh.write(f">{name}\n")
                for i in range(0, len(s), 70):
                    fh.write(s[i:i + 70] + "\n")

        def rc(s):
            return s[::-1].translate(str.maketrans("ACGT", "TGCA"))

        def pairs(src, spans, tag, n_per):
            r1, r2, i = [], [], 0
            for lo, hi in spans:
                for _ in range(n_per):
                    p = rng.randint(lo, hi - 300)
                    frag = src[p:p + 300]
                    r1.append(f"@{tag}{i}\n{frag[:100]}\n+\n{'I' * 100}\n")
                    r2.append(f"@{tag}{i}\n{rc(frag)[:100]}\n+\n{'I' * 100}\n")
                    i += 1
            return "".join(r1), "".join(r2)

        a1, a2 = pairs(c1, [(500, 1500), (4000, 5000)], "a", 120)
        b1, b2 = pairs(c2, [(1000, 1800)], "b", 120)
        for name, data in (("s1_R1.fq.gz", a1), ("s1_R2.fq.gz", a2),
                           ("s2_R1.fq.gz", b1), ("s2_R2.fq.gz", b2)):
            with gzip.open(os.path.join(d, name), "wt") as fh:
                fh.write(data)

        subprocess_run = subprocess.run
        subprocess_run(["bwa", "index", ref], check=True, capture_output=True)
        s1_bam = os.path.join(d, "s1.bam")
        mem = subprocess.Popen(
            ["bwa", "mem", ref, os.path.join(d, "s1_R1.fq.gz"), os.path.join(d, "s1_R2.fq.gz")],
            stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        subprocess_run(["samtools", "sort", "-o", s1_bam, "-"], stdin=mem.stdout, check=True)
        mem.stdout.close()
        mem.wait()
        subprocess_run(["samtools", "index", s1_bam], check=True)

        manifest = os.path.join(d, "m.tsv")
        with open(manifest, "w") as fh:
            fh.write("name\tr1\tr2\tbam\n")
            fh.write(f"s1\t{os.path.join(d, 's1_R1.fq.gz')}\t{os.path.join(d, 's1_R2.fq.gz')}\t{s1_bam}\n")
            fh.write(f"s2\t{os.path.join(d, 's2_R1.fq.gz')}\t{os.path.join(d, 's2_R2.fq.gz')}\t-\n")
        cls.ref = ref
        cls.manifest = manifest
        cls.seqs = {"c1": c1, "c2": c2, "c3": c3}

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls._tmp_root, ignore_errors=True)

    def run_pipeline(self, outname, extra):
        out = self.path(outname)
        rc = xr.main(["--ref", self.ref, "--manifest", self.manifest, "--outdir", out,
                      "--threads", "2", *extra])
        self.assertEqual(rc, 0)
        return out

    def check_run(self, out, spacer):
        with open(os.path.join(out, "contigs_kept.txt")) as fh:
            kept = [l.strip() for l in fh if l.strip()]
        self.assertEqual(kept, ["c1", "c2"])  # c3 has no mapped reads

        subset = read_fasta_dict(os.path.join(out, "subset_ref.fa"))
        self.assertEqual(sorted(subset), ["c1", "c2"])

        # every BED row maps subset bases back to the original reference
        rows = read_bed(os.path.join(out, "coverage_intervals.bed"))
        for sc, ss, se, oc, os_, oe in rows:
            self.assertEqual(subset[sc][ss:se], self.seqs[oc][os_:oe])

        # spacer: N runs appear exactly between intervals of the same contig
        for contig in ("c1", "c2"):
            crows = [r for r in rows if r[0] == contig]
            if spacer == 0:
                self.assertNotIn("N", subset[contig])  # random ACGT reference, no N
                gaps = []
            else:
                gaps = [int(crows[i + 1][1]) - int(crows[i][2]) for i in range(len(crows) - 1)]
                self.assertTrue(all(g == spacer for g in gaps), gaps)
                for i in range(len(crows) - 1):
                    gap_start = int(crows[i][2])
                    self.assertEqual(subset[contig][gap_start:gap_start + spacer], "N" * spacer)
            self.assertEqual(len(subset[contig]),
                             sum(int(r[2]) - int(r[1]) for r in crows) + spacer * len(gaps))

        for s in ("s1", "s2"):
            bam = os.path.join(out, s, "remapped.sorted.bam")
            self.assertTrue(os.path.exists(bam + ".bai"))
            subprocess.run(["samtools", "quickcheck", bam], check=True)
            mapped = int(subprocess_output(["samtools", "view", "-c", "-F", "4", bam]))
            self.assertGreater(mapped, 0, s)
            header_sq = subprocess_output(["samtools", "view", "-H", bam])
            self.assertIn("SN:c1", header_sq)
            self.assertNotIn("SN:c3", header_sq)

    def test_spacer_zero_streaming_and_in_memory_identical(self):
        a = self.run_pipeline("s0_stream", ["--stream-threshold", "0"])
        b = self.run_pipeline("s0_mem", ["--stream-threshold", "1e9"])
        with open(os.path.join(a, "subset_ref.fa"), "rb") as fa, \
                open(os.path.join(b, "subset_ref.fa"), "rb") as fb:
            self.assertEqual(fa.read(), fb.read())
        self.check_run(a, 0)

    def test_spacer_nonzero_streaming_and_in_memory_identical(self):
        a = self.run_pipeline("s9_stream", ["--stream-threshold", "0", "--spacer", "9"])
        b = self.run_pipeline("s9_mem", ["--stream-threshold", "1e9", "--spacer", "9"])
        with open(os.path.join(a, "subset_ref.fa"), "rb") as fa, \
                open(os.path.join(b, "subset_ref.fa"), "rb") as fb:
            self.assertEqual(fa.read(), fb.read())
        self.check_run(a, 9)


def subprocess_output(cmd):
    return subprocess.run(cmd, check=True, capture_output=True, text=True).stdout.strip()


if __name__ == "__main__":
    unittest.main()
