# extract_remap

Excise a reference genome down to the regions covered by one or more samples,
then remap each sample's reads onto that smaller reference.

Only reference contigs that at least one sample mapped to are kept. Within those
contigs, the covered regions of all samples are unioned, padded, and cut out. The
excised pieces of each contig are pasted back together, in order, into one record
with the original contig name.

Given:
- a reference FASTA,
- a sample manifest (TSV): one row per sample, with paired or single-end reads and
  optionally a BAM aligned to the reference,

it writes:
- `subset_ref.fa` (+ `.fai`, bwa index): the excised reference, shared by all samples,
  one record per kept contig,
- `coverage_intervals.bed`: coordinate map from subset positions to original positions,
- `contigs_kept.txt`: reference contigs that at least one sample mapped to,
- `<sample>/remapped.sorted.bam` (+ `.bai`): the sample's originally mapped reads aligned to the subset reference,
- `<sample>/reads_subset_*.fastq.gz`: the reads used for that sample's remap,
- `<sample>/first_pass.bam` (+ `.bai`): only for samples without a BAM.

## Requirements

- Python 3 (standard library only)
- `samtools`
- `bwa` (the `bwa mem` mapper)

## Usage

```
python extract_remap.py --ref REF.fa --manifest samples.tsv --outdir OUT
```

### Manifest

A tab-separated file with a header line. Columns:

| Column | Required | Meaning |
|---|---|---|
| `name` | yes | Sample name. Letters, digits, `.`, `_`, `-`. Unique. Used for the output subdirectory and read group. |
| `r1` | yes | Reads R1 (or single-end reads). Plain or gzip. |
| `r2` | no | Reads R2. If set, the sample is paired. Empty or `-` = single-end. |
| `bam` | no | BAM aligned to `--ref`. Coordinate-sorted and indexed. Empty or `-` = no BAM. |

Example:

```
name    r1              r2              bam
s1      s1_R1.fq.gz     s1_R2.fq.gz     s1.bam
s2      s2.fq.gz        -               -
s3      s3_R1.fq.gz     s3_R2.fq.gz     -
```

All paths are checked when the manifest is loaded. Relative paths are resolved from
the directory you run the command in.

### Options

| Option | Default | Meaning |
|---|---|---|
| `--ref` | required | Full reference FASTA (plain or `.gz`). |
| `--manifest` | required | Sample manifest TSV (see above). |
| `--outdir` | required | Output directory (created if missing). |
| `--min-depth` | 1 | Minimum depth, per sample, for a position to count as covered. |
| `--flank` | 20 | Bases added on each side of a covered region (clamped to contig ends). |
| `--spacer` | 0 | Number of `N` bases inserted between excised intervals of the same contig. `0` joins them directly. |
| `--threads` | 4 | Threads for `bwa mem` and `samtools sort`. |
| `--stream-threshold` | 500 | Reference size in MB above which the FASTA is streamed instead of loaded into memory. `0` always streams. |

## How it works

1. **Coverage BAM per sample.** A sample with a `bam` uses it. A sample without one
   is mapped to the full reference first, and the result is written to
   `OUT/<sample>/first_pass.bam`. The full-reference index is built only if needed.
2. **Originally mapped reads and contigs.** For each sample, primary mapped records
   are read (unmapped, secondary and supplementary records are ignored). This gives
   the sample's read names and the contigs it mapped to.
3. **Kept contigs.** A reference contig is kept if at least one sample has a primary
   mapped record on it. Other contigs are dropped. The list is written to
   `contigs_kept.txt`.
4. **Covered intervals.** For each sample, `samtools depth` is streamed on the kept
   contigs. Positions with depth >= `--min-depth` form runs. The runs of all samples
   are combined, padded by `--flank`, clamped to contig bounds, and merged when they
   overlap or touch. The result is one set of intervals for all samples.
5. **Subset reference.** For each kept contig, its intervals are concatenated in
   coordinate order into one record with the original contig name. With `--spacer N`,
   a run of N `N` bases is inserted between consecutive intervals (none before the
   first or after the last). With the default `--spacer 0` they are joined directly.
   The mapper builds one index on this subset reference.
6. **Remap per sample.** Each sample's reads that were primary-mapped in step 2 are
   written to `reads_subset_*.fastq.gz` (pairs are kept if either mate qualifies, and
   both mates are written). They are aligned to the subset reference, and the result
   is written as `<sample>/remapped.sorted.bam` (+ `.bai`).

### Edge cases

- **No reads for a sample after filtering:** the sample gets a header-only BAM
  (with `@SQ` lines for the subset contigs). It is not skipped.
- **No sample mapped anywhere:** the run stops with an error.
- **A BAM contig not in the reference:** the run stops with an error.

## Output coordinates

`coverage_intervals.bed` has one line per excised interval, with a `#` header line:

```
#subset_contig  subset_start  subset_end  orig_contig  orig_start  orig_end
c1              0             1032        c1           1480        2512
c1              1032          2067        c1           8982        10017
```

Columns: subset contig, subset start, subset end, original contig, original start,
original end. Coordinates are 0-based, half-open (BED convention). Each interval
occupies `[subset_start, subset_end)` on the subset contig, and the same-length
block `[orig_start, orig_end)` on the original contig. Spacer `N` bases are not in
any row: with `--spacer S`, the next interval on the same contig starts `S` bases
after the end of the previous one.

To map a subset position `p` (0-based) on contig `C` back to the original reference,
find the row where `subset_contig == C` and `subset_start <= p < subset_end`, then
use `orig_start + (p - subset_start)`.

With `--spacer 0` (default), intervals are joined without spacers, so a read that
spans a junction between two intervals can align across unrelated sequence. Use
`--spacer N` (for example 100 or more) so that such junction-spanning alignments
have to pass through an `N` run instead. Reads that map entirely within one interval
are not affected.

To inspect junctions in a browser, load the output in the separate mapping_tracker
app (Browser tab). It marks each junction and the reads that cross one, in subset
coordinates. The app's Run tab can also run `map_references.py` and then this
script for a chosen reference, writing the output to `<output-dir>/extract/<pathogen>/`.

## Large references

The reference is loaded fully into memory when it is at or below `--stream-threshold`.
Above that, the tool reads it in two streaming passes: one for contig lengths, and
one that writes only the intervals. Both modes produce byte-identical `subset_ref.fa`.

Streaming only reduces memory for building the subset FASTA. `bwa index` still loads
the whole subset reference, and a run with any sample lacking a BAM also runs
`bwa index` on the full reference. Those are the real memory limits for large genomes.

## Design: swappable mapper

The code has three layers:

- **Utilities**: FASTA/FASTQ I/O, BAM/samtools calls, interval logic, manifest loading.
  No aligner code.
- **Mapper**: the abstract `Mapper` class (`index`, `align`). The only implementation
  is `BwaMemMapper`, registered in `MAPPERS`.
- **Pipeline**: `main()`, which connects the utilities and the mapper.

To add another aligner, subclass `Mapper`, implement `index()` and `align()` (a sorted
BAM output), register it in `MAPPERS`, and set `DEFAULT_MAPPER`. The rest of the
pipeline does not change.

## Limitations

- Input BAMs must be coordinate-sorted, indexed, and aligned to the same `--ref`
  (contig names must match).
- Read names are matched after stripping a leading `@`/`>` and any `/1` or `/2`
  suffix. Paired mates must share a name.
- Read groups are per sample: `ID` and `SM` are the sample name, `PL:ILLUMINA`.
- The `--min-depth` and `--flank` settings apply to all samples in the run.
- No command-line option selects the mapper yet. Only `bwa mem` is available.

## Tests

The test suite uses the standard library `unittest` and lives in `tests/`:

```
python -m unittest discover -s tests -v
```

Coverage:

- Coverage logic: run merging and `--min-depth`, contig boundaries, padding, clamping,
  merging, and the union of several samples.
- Subset sequence: spacer insertion, including a spacer across a 60-column wrap.
- Streaming extraction against slicing on random interval sets, with mixed line
  widths, lowercase, N, and no trailing newline, for spacer 0, 1 and 5.
- BED: offsets with and without spacer, and a check that each row maps back to the
  original reference.
- Manifest validation: valid input, `-`/empty cells, duplicate or invalid names,
  missing files, missing `r1`, missing header.
- Argument validation for `--spacer`.
- End-to-end (needs `samtools` and `bwa`, otherwise skipped): two samples, three
  contigs with one unmapped. Checks kept contigs, BED mapping, spacer placement,
  quickcheck and `@SQ` lines, and that streaming and in-memory subset FASTA files are
  byte-identical, with spacer 0 and 9.

Earlier manual checks on synthetic data (not shipped):

- Three samples on a three-contig reference: one with a BAM, one without, and one
  with only unmappable reads. The unmapped contig was dropped, the union gave two
  intervals, the two mapped samples remapped all their reads, and the unmappable
  sample got a header-only BAM. All BAMs passed `samtools quickcheck`.
- Two intervals on one contig: the subset contig was rebuilt from the two intervals
  and matched the reconstruction from the BED and the original reference. Streaming
  and in-memory modes produced identical subset FASTA files.
- Streamed extraction compared with plain slicing on random interval sets over a
  reference with mixed line widths, lowercase, N, and no trailing newline.

## License

GPL-3.0. See [LICENSE](LICENSE).
