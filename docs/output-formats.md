# Output formats

This page describes every file that `cansnper2 type` writes, with the columns of each table. It is a
reference for users who read the results with scripts or in a pipeline. How to read a call is
described in [interpreting-results.md](interpreting-results.md), and the files that the `db`,
`scheme` and `ref` verbs write are described in [schemes-and-databases.md](schemes-and-databases.md).

The excerpts below come from typing three *Francisella tularensis* assemblies and one simulated
paired-end read set against the `francisella_tularensis` scheme with the default aligner,
minimap2; paths are shortened.

## Directory layout

```text
results/                              --outdir (default results)
  summary.html                        overview of all samples (every run with a sample)
  snps_summary.txt                    table of all samples (--summary)
  summary_tree.html                   tree with the final SNPs of all samples (--summary)
  samples/
    <sample>_snps.txt                 derived canSNPs and the call
    <sample>_not_called.txt           all other canSNPs
    <sample>_site_counts.tsv          base counts per canSNP (reads only)
    <sample>_tree.html                tree of one sample (--save-tree)
logs/                                 --log-dir (default logs)
  cansnper2-type-<Mon-DD-YYYY>.log    log of the run
  <reference>_<sample>.<aligner>.log  output of each alignment
```

`--outdir`, `--log-dir` and `--tmp-dir` are resolved against `--workdir` (default the current
folder). The `samples` folder is created when the first per-sample file is written.
`--no-snp-files` leaves out `<sample>_snps.txt`, `<sample>_not_called.txt` and
`<sample>_site_counts.tsv`; the call is still printed and written to `summary.html`.

`summary.html` links to the files in `samples/` by relative paths; keep the folder together
when it is moved.

## Terminal output

For each sample, `cansnper2 type` prints one line with the sample name and the final SNP, or
`NA` when the sample was not called:

```text
GCA_000018925.1: M.1
GCA_000746455.1: A.I.14
GCA_000807905.1: B.136
```

Warnings and errors are printed on standard error and written to the log. Warnings about canSNP
positions without alignment data are written only to the log, and one line at the end gives
their number and the log file (see [troubleshooting.md](troubleshooting.md)).

## Per-sample canSNP files

`<sample>_snps.txt` and `<sample>_not_called.txt` are tab-separated tables with the same
columns. Each row is the typed canSNP of one node (the canSNP that CanSNPer2 calls); the
alternative canSNPs of a node are not read. Together the two files list every node once:

- `<sample>_snps.txt` holds the canSNPs at which the sample has the derived base. This includes
  derived canSNPs outside the SNP path; the path itself is given in the `# SNP path:` line.
- `<sample>_not_called.txt` holds all other canSNPs: ancestral, other base, gaps and `N`, and,
  at the end, the canSNPs without any alignment data.

| Column | Content |
|---|---|
| `Name` | Name of the node whose typed canSNP the row gives. |
| `Reference` | Reference genome the canSNP is defined on. |
| `Sequence` | Sequence of that reference genome (first word of its FASTA header). |
| `Pos` | 1-based position of the canSNP on that sequence. |
| `QuerySeq` | Assemblies: the query sequence (contig) the base was read from. Reads: `reads`. `N/A` when there was no alignment data. |
| `QueryPosition` | Assemblies: position of the base in that query sequence, `N/A` for a gap. Reads: the number of reads counted at the site (the depth). `N/A` when there was no alignment data. |
| `Ancestral base` | Ancestral base recorded by the scheme. |
| `Derived base` | Derived base recorded by the scheme. |
| `Target base` | Base of the sample, on the forward strand of the reference. `-` is a gap, `N` an ambiguous base or, for reads, a site whose counts did not decide a base. `N/A` when there was no alignment data. |
| `Edge flag` | `edge` when the canSNP lies within 25 bases of either end of the alignment block (assemblies only); otherwise empty. |

After the table follow lines that start with `# `, so that a reader that skips comment lines
(for example `pandas.read_csv(path, sep="\t", comment="#")`) reads the table alone.
`<sample>_snps.txt` ends with these lines, in this order:

| Line | Content |
|---|---|
| `# SNP path:` | The SNP path, `;`-separated, from the root to the final SNP. For a sample that was not called, all derived canSNPs, or `NA`. |
| `# Final SNP:` | The final SNP and `found/depth:`, the number of derived canSNPs on the path and the depth of the final node in the tree (the number of nodes from the top of the tree to it, both included, label nodes counted). `# Final SNP: NA` for a sample that was not called. See the note below the table. |
| `# Reads:` | Read samples only; see [Reads line](#reads-line). |
| `# Query file:` | Absolute path of the query file (both files of a read pair, space-separated) and `sha256:` of the query. For a read pair this is the sha256 of the two file sha256 values joined by a space. |
| `# Typing settings:` | `database sha256:` of the database file and the settings that change a call; see [Typing settings line](#typing-settings-line). |

`<sample>_not_called.txt` ends with the `# Query file:` and `# Typing settings:` lines.

The support of a call, `found/expected` (see
[interpreting-results.md](interpreting-results.md#support)), counts the canSNPs on the path, and
label nodes have none (see [method.md](method.md#label-nodes)). In a tree without label nodes,
which includes every verified scheme of CanSNPer2-data, the depth equals `expected`, so the
`# Final SNP:` line gives the support. In a tree with label nodes the depth is larger; the
support is then the `Support` column of `summary.html` and the `support` column of
`snps_summary.txt`.

Excerpt of `samples/GCA_000807905.1_snps.txt`:

```tsv
Name	Reference	Sequence	Pos	QuerySeq	QueryPosition	Ancestral base	Derived base	Target base	Edge flag
T/N.1	SCHUS4.2	AJ749949.2	83976	JUJV01000058.1	9850	G	A	A	
T.1	SCHUS4.2	AJ749949.2	1165688	JUJV01000046.1	5353	A	G	G	
B.1	OSU18	CP000437.1	1710718	JUJV01000058.1	9619	C	T	T	
...
B.136	FSC200	CP003862.1	350028	JUJV01000122.1	8867	C	T	T	
# SNP path: T/N.1;T.1;B.1;B.2;B.3;B.5;B.6;B.7;B.133;B.81;B.136
# Final SNP: B.136 found/depth: 11/11
# Query file: /data/project/genomes/GCA_000807905.1.fasta sha256:57f4df3e9755607532bcf142ebdd9c2cbb4b06ddbe296fe19473887531f76849
# Typing settings: database sha256:f3fc3d4822ae378ebb5922afcb3c41a311c1ac9d0102dc9395f508a02fa104a2 min_required_hits=3 strictness=0.7 mask=0 aligner=minimap2
```

Excerpt of `samples/GCA_000807905.1_not_called.txt`, with an ancestral canSNP, an ancestral
canSNP flagged at an alignment edge, and a canSNP without alignment data:

```tsv
Name	Reference	Sequence	Pos	QuerySeq	QueryPosition	Ancestral base	Derived base	Target base	Edge flag
B.187	FSC200	CP003862.1	50301	JUJV01000028.1	4126	C	T	C	
...
B.85	OSU18	CP000437.1	946694	JUJV01000118.1	11780	G	T	G	edge
...
B.36	OSU18	CP000437.1	557870	N/A	N/A	G	T	N/A	
# Query file: /data/project/genomes/GCA_000807905.1.fasta sha256:57f4df3e9755607532bcf142ebdd9c2cbb4b06ddbe296fe19473887531f76849
# Typing settings: database sha256:f3fc3d4822ae378ebb5922afcb3c41a311c1ac9d0102dc9395f508a02fa104a2 min_required_hits=3 strictness=0.7 mask=0 aligner=minimap2
```

### Typing settings line

The `# Typing settings:` line holds `database sha256:<sha256>`, `min_required_hits`, `strictness`
and `mask`. For assemblies it adds `aligner`; for read samples it adds `platform`, `min_depth`, `min_allele_fraction`, `min_mapq`,
`min_baseq` and `max_reads` (the values used, whether given or taken from the platform defaults;
`max_reads=0` means that all reads were mapped), and with
`--local-references` also `local_references=on` and `local_rounds`.

```text
# Typing settings: database sha256:f3fc3d48...a104a2 min_required_hits=3 strictness=0.7 mask=0 platform=illumina min_depth=5 min_allele_fraction=0.8 min_mapq=20 min_baseq=20 max_reads=300000
```

### Reads line

For read samples, `<sample>_snps.txt` has a `# Reads:` line before the `# Query file:` line:

| Field | Content |
|---|---|
| `platform` | Sequencing platform of the sample. |
| `files` | Number of read files (1 or 2). |
| `sites_with_reads` | canSNP sites with at least one counted read, out of all canSNP sites. |
| `median_depth` | Median depth over the sites with reads. |
| `mixed_sites` | Number of sites at which no base reached `--min-allele-fraction`. |

With `--local-references` the line adds `local_reference` (the genome used, or `none`),
`local_node`, `screen_share`, `screens`, `local_rounds` and `conflicts`, and with `--fetch-local`
also `fetched`, the number of genomes downloaded for the sample. [method.md](method.md) describes
local references.

```text
# Reads: platform=illumina files=2 sites_with_reads=324/328 median_depth=20 mixed_sites=0
```

## Site counts file (reads)

For read samples, `<sample>_site_counts.tsv` gives the base counts at every canSNP of the
database. It is a tab-separated table with a header row and no `#` lines.

| Column | Content |
|---|---|
| `Name`, `Reference`, `Sequence`, `Pos` | The canSNP and its position, as in the canSNP files. |
| `Ancestral base`, `Derived base` | Bases recorded by the scheme. |
| `A`, `C`, `G`, `T` | Number of reads with each base, both strands, after the mapping and base quality filters. |
| `depth` | Sum of the four counts. |
| `fraction_derived` | Fraction of the counted reads with the derived base, two decimals. |
| `state` | `derived`, `ancestral`, `other_base`, `mixed` (no base reached the allele fraction), `low_depth` (fewer reads than `--min-depth`), `no_reads`, and with `--local-references` also `not_lifted` and `conflict`. All states other than the first three count as no data. |
| `source` | Genome the reads were mapped to for these counts: the reference genome, or a local reference; two genomes, comma-separated, for a conflict. |
| `lifted_position` | With a local reference, the position the canSNP was lifted to, as `<sequence>:<position><strand>`; otherwise empty. |

```tsv
Name	Reference	Sequence	Pos	Ancestral base	Derived base	A	C	G	T	depth	fraction_derived	state	source	lifted_position
B.187	FSC200	CP003862.1	50301	C	T	0	0	0	0	0	0.00	no_reads	FSC200	
B.116	FSC200	CP003862.1	63364	T	G	0	0	0	29	29	0.00	ancestral	FSC200	
M.1	SCHUS4.2	AJ749949.2	75124	C	T	0	0	0	19	19	1.00	derived	SCHUS4.2	
```

The thresholds behind these states are described in [method.md](method.md#site-states-from-reads).

## summary.html

`summary.html` is written in `--outdir` at the end of every run that had at least one sample; a
run that stops with an error writes none. It is a single HTML file that opens in a web browser without a network connection. It shows:

- the run information: CanSNPer2 version, date, database and its sha256, scheme, reference
  genomes, aligner and typing settings (the same lines as in `snps_summary.txt`);
- a table with one row per sample: `Sample`, `Final SNP`, `Support`, `No data` (canSNPs on the
  SNP path without alignment data) and `Status`, with the samples that need attention first;
- the scheme tree, painted with the site states of the sample selected in the table, and the SNP
  path of that sample with the bases at each canSNP;
- links to the per-sample files of the selected sample.

`Status` is the call status: `clean`, `review` or `not called`, and in addition `failed` (the
alignment or processing of the sample failed) and `reused` (the result was read from an earlier
run). The rules are given in [method.md](method.md#call-status). Each `review` and `not called`
row lists its reasons.

The button `Download table (TSV)` saves the table as a tab-separated file with the run
information as `#` lines and the columns `sample`, `final_snp`, `support`, `no_data`, `status`,
`snp_path` and `notes`. The button `Save SVG` saves the tree as it is drawn.

## snps_summary.txt (--summary)

With `--summary`, `snps_summary.txt` lists every sample of the run, including samples that were
not called. It starts with `# ` lines on the run, followed by a tab-separated table:

| Column | Content |
|---|---|
| `sample` | Sample name. |
| `final_snp` | Final SNP, or `NA`. |
| `support` | `found/expected`: derived canSNPs on the SNP path / canSNPs on the SNP path; `NA` when not called. |
| `snp_path` | The SNP path, `;`-separated, as in `<sample>_snps.txt`. |

```text
# CanSNPer2 version: 2.6.0
# Date: 2026-11-25 13:52
# Database: francisella_tularensis.db sha256:f3fc3d4822ae378ebb5922afcb3c41a311c1ac9d0102dc9395f508a02fa104a2
# Scheme: francisella_tularensis (verified), FOI-Bioinformatics/CanSNPer2-data master@d7493fd, built 2026-11-25 with CanSNPer2 2.6.0
# Reference genomes: OSU18 (GCA_000014605.1), FSC200 (GCA_000168775.2), FTNF002-00 (GCA_000017785.1), LVS (GCA_000009245.1), SCHUS4.2 (GCA_000008985.1)
# Aligner: minimap2
# Typing settings: min_required_hits=3 strictness=0.7 mask=0
sample	final_snp	support	snp_path
GCA_000018925.1	M.1	4/4	T/N.1;T.1;A/M.1;M.1
GCA_000746455.1	A.I.14	9/9	T/N.1;T.1;A/M.1;A.1;A.I.1;A.I.7;A.I.12;A.I.2;A.I.14
GCA_000807905.1	B.136	11/11	T/N.1;T.1;B.1;B.2;B.3;B.5;B.6;B.7;B.133;B.81;B.136
```

A `# Note:` line follows when the scheme is not verified or the database was changed after it
was built. For read runs the `# Aligner:` line reads `minimap2-reads`; the read thresholds are
recorded in the `# Typing settings:` line of each `<sample>_snps.txt`, not here.

`--summary` also writes `summary_tree.html`, the scheme tree with the final SNPs of all samples
marked.

## Tree pages (--save-tree)

With `--save-tree`, `samples/<sample>_tree.html` shows the scheme tree of one sample, painted
with its site states, with its final SNP, support and SNP path. `summary.html` shows the same
tree for every sample, so these pages are optional. Like `summary.html`, they are single HTML
files that need no network connection.

## Log files

The log files are written to `--log-dir` (default `logs`, relative to `--workdir`).

- `cansnper2-type-<Mon-DD-YYYY>.log` is the log of a typing run, for example
  `cansnper2-type-Nov-25-2026.log`. When that name is taken, the time is added as `HHMM`
  (`cansnper2-type-Nov-25-2026-1404.log`), and then a number (`...-1404.2.log`). By default it holds
  warnings and errors; with `-v` also progress information, and with `--debug` debug messages.
  With `--quiet` it holds errors only.
- `<reference>_<sample>.<tag>.log` holds the output of the aligner for one sample and one
  reference genome. The tag is `mauve`, `minimap2`, `mummer4`, or `reads` for read mapping.
  With `--local-references`, `<genome>_<sample>.local.log` holds the mapping to a local reference.

The other verbs name their log after the verb in the same way, in their own `--log-dir`
(default `logs` in the current folder): `cansnper2-db-create-<date>.log` (and `db-import`,
`db-export`, `db-info`, `db-convert`, `db-modify`, `db-remove`), `cansnper2-scheme-get-<date>.log`
and `cansnper2-ref-get-<date>.log`. `-v`, `--debug` and `--quiet` set their level in the same way.
`cansnper2 scheme list` writes no log, and `cansnper2 selftest` none either.

Error messages are printed to standard error; results, listings and progress messages to
standard output.
`scheme list` and `selftest` write no log file.

## Temporary files

Each run aligns in its own folder `run_<random>` below a `CanSNPer2` folder. Without
`--tmp-dir` this is the system temporary folder (which follows the `TMPDIR` environment
variable), or `<outdir>/.cansnper2_tmp` when that cannot be written. With `--tmp-dir DIR` it is
`DIR/CanSNPer2/run_<random>`.

The run folder is removed at the end of the run, and also when the run stops with an error.
With `--keep-temp` it is kept and its path is printed:

```text
Alignments and temporary files are kept in /data/project/tmp/CanSNPer2/run_t9t5s7qb
```

It holds one alignment file per reference genome and sample, `<reference>_<sample>.<suffix>`,
with the suffix `xmfa` (progressiveMauve), `paf` (minimap2), `delta` (MUMmer4) or `pileup`
(reads, with the BAM file and the site lists next to it), and `<sample>.query_sha256`, the
sha256 of the query the alignments were made from. progressiveMauve alignments have a
`.sequences` table next to them. A kept folder can be read again without aligning:

```sh
cansnper2 type GCA_000018925.1.fasta -d francisella_tularensis.db --skip-alignment \
    --tmp-dir tmp/CanSNPer2/run_t9t5s7qb -o results_again
```

Give the same `--aligner` as the run that made the alignments (the default is minimap2), since
the aligner determines the file names that are read. With `--skip-alignment`, CanSNPer2 refuses alignments that were made
from another query file with the same sample name.

## Results from an earlier run

When `<sample>_snps.txt` already exists in the output folder, the sample is not typed again: its
call is read from the file, printed with `(from an earlier run, use --rerun to type it again)`,
and shown as `reused` in `summary.html`. The earlier file is used only when its `# Query file:`
sha256 and its `# Typing settings:` line match the current run; otherwise the run stops and asks
for `--rerun` or another `--outdir`. `--rerun` types all samples again and replaces their files.

For assemblies the aligner is part of the `# Typing settings:` line, so a run with another
`--aligner` into the same output folder stops instead of reusing the result. A result file
written before 2.6.0, which does not record the aligner, stops the run in the same way.
To compare aligners, type into a separate `--outdir` for each.
