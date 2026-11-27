# How CanSNPer2 makes a call

This page describes, from the code of CanSNPer2 2.6.0, how a sample is typed: how the base at
each canSNP position is read from an alignment or from mapped reads, and how the tree walk turns
those site states into a final SNP and a call status. It is written for users who need to know
the exact rules, for example to judge a result or to build a scheme. For reading results in
practice, see [interpreting-results.md](interpreting-results.md).

## Terms

| Term | Definition |
|---|---|
| scheme | The published or user-made set of text files: tree, canSNPs, reference genomes and an optional genomes file. |
| database | The SQLite `.db` file built from a scheme; the file that `cansnper2 type -d` reads. |
| canSNP | Canonical SNP: a SNP whose derived allele defines a branch of the scheme tree. |
| node | A position in the scheme tree, named after its canSNP (for example `B.Br.001`). A node can have several canSNPs, ranked from 1 within the node; exactly one of them is marked `typed`, and CanSNPer2 calls only that one. |
| ancestral / derived | The two allele states of a canSNP that the scheme records. |
| other base | A base that is neither the ancestral nor the derived allele. |
| no data | No aligned base (or too few reads) at a canSNP position. |
| site state | The per-canSNP result for a sample: derived, ancestral, other base or no data (for reads, also mixed, which counts as no data). |
| SNP path | The derived canSNPs from the root to the final node. |
| final SNP | The result of typing: the deepest node whose path is confirmed, as named in the output files. |
| support | `found/expected`: derived canSNPs on the path / canSNPs on the path. |
| call status | `clean`, `review` or `not called`, as shown in `summary.html` (see [Call status](#call-status)). |
| query | An input file given to `cansnper2 type`. |
| sample | The unit that is typed and reported: one assembly, or one read set (single file or pair). |
| reference genome | A genome of the scheme on which canSNP positions are defined (`<genome>.fna` in `--ref-dir`). |
| representative genome | A genome of the genomes file, used by read typing with `--local-references`. |

## Overview

A database holds the scheme tree, the typed canSNP of each node, and for each canSNP the
reference genome, the reference sequence and the 1-based position on that sequence where it is
defined, with its ancestral and derived base. A node can carry several canSNPs, ranked from 1
within the node (`brucella_genus`, for example, has up to three per node). Exactly one canSNP
of each node is marked `typed`; CanSNPer2 reads and calls only the typed canSNPs. The others
are alternatives, which are kept in the database and written back by `cansnper2 db export`. A
database in which no canSNP is marked `typed` is refused. For each sample, CanSNPer2 aligns the assembly
(or maps the reads) to every reference genome of the database that carries canSNPs, reads the
query base (or the base counts) at each canSNP position of that reference, and classifies each
canSNP as derived, ancestral, other base or no data. Each typed canSNP is read from exactly one
reference genome, the one the database assigns it to. The tree walk then takes the derived
canSNPs, from the deepest node upwards, and reports as final SNP the deepest node whose path
from the root is confirmed by the rules described below. A sample without a confirmed node is
reported as `NA`.

## Site states from assemblies

### Preparing the inputs

Before alignment, each query and reference FASTA file is checked and, where needed, a cleaned
copy is written to the run directory (the original file is not changed):

- Every sequence character that is not A, C, G, T (upper or lower case) or N is replaced with
  N. This covers gap characters (`-`, `.`, `*`) and IUPAC ambiguity codes; an ambiguity code at
  a canSNP position therefore gives no data, not an other base. Positions are not shifted.
- Blank lines and whitespace inside sequence lines are removed. A byte order mark is removed.
- A gzip-compressed query is refused, as are an empty query, a query without a FASTA header, a
  header without sequence, and two sequences with the same name (the first word of the header).

Each reference file is checked once per run. It is refused when a sequence name is empty or
used twice, when a sequence that the database places canSNPs on is missing, when a sequence is
shorter than its largest canSNP position, or when the reference base is neither the ancestral
nor the derived base at one fifth or more of its canSNPs. Fewer mismatches give a warning.
Only reference genomes of the database that carry typed canSNPs are aligned against; other
`.fna` files in `--ref-dir` are skipped with a warning. A reference genome of the database
without a file stops the run unless `--allow-missing-references` is given; its canSNPs then
have no data.

### Running the aligner

The query is aligned against each reference genome in a separate process; all alignments of
one query run at the same time. A round of alignments that takes longer than
`--alignment-timeout` (default 60 minutes; 0 for no limit) is stopped. The aligner is chosen
with `--aligner`; minimap2 is the default from 2.6.0 (progressiveMauve was the default up to
2.4.0). The command lines are defined in
`CanSNPer2/modules/Aligners.py`.

| `--aligner` | Command | Output |
|---|---|---|
| `minimap2` (default) | `minimap2 -c --cs=long -x asm20 -r 1k,10k -N 50 -p 0.1 -o <out> <reference> <query>` | PAF |
| `progressiveMauve` | `progressiveMauve --output <out> <reference> <query>` | xmfa |
| `mummer4` | `nucmer -t 1 --delta=<raw> <reference> <query> && delta-filter -1 <raw> > <out>` | delta |

The settings, and the reasons given for them in the code:

- minimap2 uses the assembly preset `asm20`, which allows for up to about 5 % sequence
  divergence. The long-join bandwidth is lowered from 100 kb to 10 kb (`-r 1k,10k`), and
  secondary chains are kept down to 10 % of the best score (`-N 50 -p 0.1`, instead of 5 and
  0.8). Both changes were made so that short stretches between IS elements and rearrangements,
  which separate Francisella genomes of different subspecies, are aligned; with the preset
  values more canSNP positions had no data (see [validation.md](validation.md)).
- progressiveMauve runs with its default settings. It writes a progress file (`.sslist`) next
  to each input, so both genomes are linked into the run directory before alignment. When
  progressiveMauve exits with code 11 or is aborted (signal 6), the alignment is retried once
  with every `-` in the query replaced by N.
- MUMmer4 `nucmer` runs with its default settings except one thread (`-t 1`): with the default
  of two threads, repeated runs on the same genomes gave different alignments after the
  filter at a few canSNP positions, while one thread gave the same output each time.
  `delta-filter -1` keeps, for every reference and query position, only the best alignment,
  allowing rearrangements. `nucmer` and `delta-filter` run as one shell command, so a failure
  of either fails the alignment.

### Reading the alignment

For every canSNP position that an aligned block of the reference covers, the reader of the
aligner reports the query base on the forward strand of the reference:

- minimap2 (PAF): only primary alignments (`tp:A:P`) with a mapping quality of at least 20 are
  read. A secondary alignment, or a primary one with a lower mapping quality, can carry the
  base of another copy of a repeated region. The query base is taken from the `cs` tag, which
  is on the reference strand.
- progressiveMauve (xmfa): the reference is the first sequence of each block and the query the
  second. Blocks that hold only one of the two genomes carry no comparison and are skipped. On
  a block where the reference runs on the reverse strand, the query base is complemented.
- MUMmer4 (delta): the delta file holds coordinates and indel positions but no bases, so the
  query bases are read from the query FASTA file named in the delta file; that file must still
  exist when the alignment is read again with `--skip-alignment`. Bases of reverse alignments
  are complemented.

A query gap at a canSNP position is recorded as `-`. When a position is covered by more than
one aligned block with different query bases (possible with minimap2 and MUMmer4 in repeats or
duplicated contigs, not with progressiveMauve), it is recorded as N and a warning names it.

### Site states

The recorded base is compared with the bases of the canSNP:

| Query base | Site state | Value in the code |
|---|---|---|
| the derived base | derived | 1 |
| the ancestral base | ancestral | 2 |
| another of A, C, G, T | other base | 3 |
| `-` (gap) or N | no data | 0 |
| no aligned block covers the position | no data (the canSNP is absent from the states) | - |

canSNPs that no alignment covered are listed in `<sample>_not_called.txt` with the target base
`N/A`, so that every canSNP of the tree appears in one of the two per-sample result files (see
[output-formats.md](output-formats.md)).

### Positions near alignment edges and `--mask`

A canSNP closer than 25 bases to either end of the aligned block it lies in is flagged `edge`
in the result files and listed in a warning. The flag does not change the site state; an edge
flag on a canSNP of the SNP path sets the call status to `review`. Ends of contigs and of
alignment blocks are where alignment errors are most likely.

`--mask N` (default 0, off) leaves out canSNPs within N bases of either end of an aligned
block: they are not evaluated, so in the tree walk they count as no data. The number of masked
positions per reference is logged. `--mask` applies to assemblies only.

### References with more than one sequence

Since CanSNPer2 2.3.0, a canSNP position is stored as a reference sequence (the first word of
its FASTA header, a GenBank accession in the published schemes) and a 1-based position on that
sequence. minimap2 and MUMmer4 report positions per sequence, so their alignments are read
directly. progressiveMauve counts positions along the reference sequences joined end to end;
before alignment CanSNPer2 writes the name and length of each reference sequence, in file
order, to `<alignment>.sequences`, and the xmfa reader converts the positions back with this
table. The order of the sequences in the reference file and sequences without canSNPs (such as
plasmids) therefore do not affect the result. The sequence names in the reference file must be
the ones the database uses; a reference with other names (for example RefSeq instead of
GenBank accessions) is refused by the reference check. Databases with positions along joined
sequences, from before 2.3.0, are described in [upgrading.md](upgrading.md).

## Site states from reads

When the queries are FASTQ files, each sample is mapped to each reference genome instead of
aligned. A run types either reads or assemblies, not both.

### Mapping and pileup

For each reference genome the following runs as one `bash -o pipefail` command, so that a
failure of any step fails the sample:

```text
minimap2 -a -x <preset> -t <threads> <reference> <(head -n <4 x max-reads> <reads>)...
  | samtools sort -@ <threads> -o <out>.bam -
samtools index <out>.bam
samtools mpileup -l <out>.sites.bed -a -B -A -d 0 -q <min-mapq> -Q <min-baseq> -o <out> <out>.bam
```

The BED file lists the canSNP positions of that reference, so the pileup covers only those
positions. The `mpileup` options are chosen as follows: `-B` turns off base alignment quality,
which lowers base qualities next to indels in a way that depends on the reference; `-A` counts
pairs that are not properly paired, which are frequent against a distant reference; `-d 0`
removes the depth limit; `-a` writes a line for every position, also without reads. Reads below
`--min-mapq` and bases below `--min-baseq` are not counted. Secondary, QC-failed and
duplicate-flagged reads are excluded by the `mpileup` defaults, and a base covered by both reads
of a pair is counted once.

Only the first `--max-reads` reads of each read file are mapped (default 300000; for paired-end
reads, 300000 read pairs). Each file is cut with `head` in a process substitution, after
`gzip -cd` when it is compressed, so no shortened copy is written to disk. The first reads of a
file written by a sequencer are distributed over the whole genome. A file sorted by position,
such as one converted from a sorted BAM file, is not; for such a file the first reads cover only
part of the genome, and `--max-reads 0` maps all reads. For 2 x 150 base reads, 300000 read pairs
correspond to about 90 Mb, or 45x on a 2 Mb genome. For long reads the same number of reads is
usually the whole run.

`--threads` (default 4) is the number of threads for one sample. A sample is mapped to all
reference genomes of the scheme at the same time, and the threads are divided between these
mappings, at least one each: with 4 threads and the five reference genomes of
`francisella_tularensis`, each mapping uses one thread for minimap2 and for `samtools sort`. A
mapping to a local reference runs alone and uses all threads.

### Platform defaults

`--platform` sets the minimap2 preset and the default thresholds (`ParsePileup.PLATFORMS`).
Options given on the command line replace the defaults.

| `--platform` | minimap2 preset | `--min-mapq` | `--min-baseq` | `--min-depth` | `--min-allele-fraction` | `--max-reads` |
|---|---|---|---|---|---|---|
| `illumina` (default) | `sr` | 20 | 20 | 5 | 0.8 | 300000 |
| `ont` | `map-ont` | 20 | 10 | 5 | 0.7 | 300000 |
| `hifi` | `map-hifi` | 20 | 20 | 5 | 0.8 | 300000 |

### Site state of a canSNP from its counts

The A, C, G and T counts at a canSNP position (bases that passed the quality filters) are
classified in this order:

| Counts | Site state | Label in `<sample>_site_counts.tsv` |
|---|---|---|
| no reads | no data | `no_reads` |
| fewer than `--min-depth` bases | no data | `low_depth` |
| no base reaches `--min-allele-fraction` of the depth | no data | `mixed` |
| the derived base reaches the fraction | derived | `derived` |
| the ancestral base reaches the fraction | ancestral | `ancestral` |
| another base reaches the fraction | other base | `other_base` |

The base that reaches the fraction is the most frequent one. With local references (below) two
further labels occur, `not_lifted` and `conflict`; both count as no data.

Example, Illumina defaults (`--min-depth 5`, `--min-allele-fraction 0.8`), ancestral C,
derived T:

| Counts | Depth | Most frequent base share | Site state |
|---|---|---|---|
| T 3 | 3 | - | no data (`low_depth`) |
| T 18, C 2 | 20 | 0.90 | derived |
| T 12, C 8 | 20 | 0.60 | no data (`mixed`) |
| A 9, C 1 | 10 | 0.90 | other base |

### Mixed sites and the possibly-mixed flag

CanSNPer2 types isolates. A mixed site gives no evidence to the tree walk, like a site without
data, but it is counted. When a sample has at least `--mixed-sites` mixed canSNP sites (default
3; the count is over all canSNPs of all reference genomes, not only those on the path), the
sample is flagged as possibly mixed: a warning is logged, the call status is set to `review`
with the mixed canSNPs named, and the count is given in the `Reads:` line of
`<sample>_snps.txt`. The flag never changes the final SNP. The strains in a mixture are not
resolved. A region that the sample has twice while the reference has it once gives the same
pattern, because reads of both copies map to the single copy of the reference.

Edge flags and `--mask` are not used for reads: read ends are filtered by base quality.

## Tree walk and the final SNP

The tree walk is implemented in `NewickTree._final_call` and `NewickTree._confirm_path`
(`CanSNPer2/modules/NewickTree.py`). Its inputs are the site state of every canSNP and the
list of derived canSNPs.

### Label nodes

Each node is represented in the tree walk by its typed canSNP; alternative canSNPs of a node
are not read. A node of the database tree without a typed canSNP is a label node; it is shown as `-name-`
and is not a SNP. Label nodes are skipped when a path is checked and do not count in the
support. In the verified schemes of CanSNPer2-data every node, the root included, has a typed
canSNP, so label nodes occur only in user-made schemes.

### Rules for confirming one node

For a node with a derived canSNP, `_confirm_path` walks its actual ancestor chain, from the
node itself up to the root. Each canSNP on the chain (the node included) is one of:

- called: derived in the sample; counts as a hit;
- no data: no aligned base, a gap or N, a masked position, or for reads no reads, low depth,
  mixed, not lifted or conflict; counts against the fraction but does not break the chain;
- against: ancestral or other base; the chain is broken and the node cannot be confirmed.

The node is confirmed when all of the following hold:

1. no canSNP on the chain is ancestral or other base;
2. the number of called canSNPs on the chain is at least `--min-required-hits` (default 3);
3. the called fraction, called / canSNPs on the chain, is at least `--strictness` (default
   0.7);
4. the canSNP directly above the node (the nearest one, skipping label nodes) is called. This
   rule prevents a single derived canSNP on a neighboring branch, whose parent had no data,
   from being reported instead of the lineage supported by the other canSNPs.

Rule 4 applies only when the chain has more than one canSNP.

### Choosing the final SNP

1. If the sample has fewer derived canSNPs in the whole tree than `--min-required-hits`, no
   node is called (message: "The number of SNPs found is too small").
2. The derived canSNPs are sorted by their depth in the tree, deepest first, and each is
   checked with the rules above.
3. The first confirmed node sets the depth of the call. Further confirmed nodes at the same
   depth are collected; the check stops at the first confirmed node that is shallower.
4. If exactly one node is confirmed at that depth, it is the final SNP. Its support is
   called / canSNPs on its chain.
5. If two or more nodes are confirmed at the same depth, none is chosen: the sample is `NA`, the
   reason (in the log and in `summary.html`) names the competing paths, the result files give each
   candidate with its support (`# Competing paths:`, `competing_paths` in the JSON), and the call
   status is `not called`.
6. If no node is confirmed, the sample is `NA`.

The SNP path of a called sample is the chain of the final SNP restricted to its derived
canSNPs. Derived canSNPs outside that chain are not part of the path; they are kept in the
result table and set the call status to `review`.

### Worked example

Take a branch of six canSNPs, S1 (root) to S6, with a sibling X6 of S6 below S5. Default
settings: `--min-required-hits 3`, `--strictness 0.7`. These outcomes were checked by running
the tree walk on a test database.

| Site states | Final SNP | Support | Reason |
|---|---|---|---|
| S1-S6 derived except S3 no data; X6 ancestral | S6 | 5/6 | 5 called, 5/6 = 0.83 >= 0.7, S5 called |
| S1-S6 derived except S5 no data | S4 | 4/4 | S6 fails rule 4 (S5 has no data); S4 is the deepest confirmed node |
| S2 and S3 no data, others derived | `NA` | - | S6: 4/6 = 0.67 < 0.7; S5: 3/5 = 0.6; S4: S3 above it has no data; S1: 1 hit < 3 |
| S3 ancestral, others derived | `NA` | - | S4, S5 and S6 are broken at S3; S2 has 2 hits < 3 |
| S6 ancestral, others derived | S5 | 5/5 | S6 is not derived and is not checked |
| S1-S6 and X6 derived | `NA` | - | S6 and X6 are confirmed at the same depth |

In the first row, the call leans on the strictness tolerance for S3, which sets the call status
to `review`.

## Call status

`summary.html` gives each sample one status. For a typed sample it starts as `clean` and is
set as follows (`CanSNPer2/modules/CanSNPer2.py`); the reasons are listed for the sample in
`summary.html`.

| Status | Condition |
|---|---|
| `not called` | No final SNP (`NA`): fewer derived canSNPs than `--min-required-hits`, no derived canSNP, no node confirmed, or two confirmed nodes at the same depth. |
| `review` | A final SNP was called, and at least one of: a canSNP on the path of the final SNP had no data (the call depends on `--strictness`); a canSNP on the SNP path is flagged `edge`; one or more derived canSNPs lie outside the SNP path, on another branch; for reads, the sample has at least `--mixed-sites` mixed sites. |
| `clean` | A final SNP was called and none of the `review` conditions holds. |

A `not called` sample keeps that status; the reasons of the other conditions may be listed
with it (for example the possibly-mixed flag). Two further statuses describe samples that were
not typed in the run:

| Status | Condition |
|---|---|
| `failed` | The alignment or mapping produced no output (with `--keep-going`), or the processing of the sample failed. |
| `reused` | The result file of an earlier run in the same `--outdir` was read instead of typing the sample again. A result file made from another input file (by sha256) or with another database or settings stops the run instead; a result file from an older version that records neither is reused as it is. For assemblies the aligner is one of the settings (see [output-formats.md](output-formats.md#results-from-an-earlier-run)). `--rerun` types the sample again. |

Sources of no data at a canSNP, as named in the result files and the log:

- no aligned block covered the position (listed in `<sample>_not_called.txt` with target base
  `N/A`); this includes canSNPs of a reference genome left out with `--allow-missing-references`;
- the query has a gap or N there, or two aligned blocks gave different bases (target base `-`
  or `N`);
- the position was masked with `--mask` (counted in the log);
- for reads, the label in `<sample>_site_counts.tsv`: `no_reads`, `low_depth`, `mixed`,
  `not_lifted` or `conflict`.

`cansnper2 type` exits with status 1 when no sample of the run was typed (all `NA`, failed or
skipped) and no earlier result was reused.

## Local references

Reads map to a reference genome only where they are close to it. For a scheme that spans a
genus with one reference genome, canSNPs of lineages far from that reference may get no reads,
and the call stops above them. `--local-references` adds rounds in which the reads are mapped
to a representative genome of the scheme that is closer to the sample, and the canSNPs are
read at their positions in that genome (`CanSNPer2/modules/LocalReferences.py`,
`CanSNPer2/modules/LiftOver.py`). It applies to reads only.

### Representative genomes

The representatives are the genomes with role `representative` in the genomes file of the
scheme, each placed at a node. Their files are read from `<ref-dir>/local/<genome>.fna`.
`cansnper2 ref get --local-references` downloads them:

```sh
cansnper2 ref get -d scheme.db -o references/ --local-references
cansnper2 type -d scheme.db --ref-dir references/ --local-references s1_R1.fastq.gz s1_R2.fastq.gz
```

Representatives without a file are not used, and a warning says how many are missing. With
`--fetch-local`, the representatives that the choice below screens or chooses, and that have no
file yet, are downloaded from NCBI during typing into the same folder, where they remain for
later runs; a genome that cannot be downloaded is passed over. Without `--fetch-local`, typing
does not download.

```sh
cansnper2 type -d scheme.db --ref-dir references/ --local-references --fetch-local s1_R1.fastq.gz s1_R2.fastq.gz
```

### When a round is made

Before each round, the tree walk is run on the current site states (the provisional call). A
round is made only when it can change the call: when there is no call, or when a canSNP on the
path from the root to the call, or the canSNP of a child of the call, has no data (no reads,
low depth or not lifted). A mixed site counts as data here. When every such canSNP has data,
the call is settled and no round is made. At most `--local-rounds` rounds are made (default 2).

### Choosing the local reference

The local reference is chosen by walking down the tree from the provisional call (from the root
when there is no call):

1. At each node, the child subtrees that hold representatives are the options; a subtree whose
   own canSNP is ancestral is left out.
2. If the derived canSNPs lie in exactly one option, or only one option exists, the walk goes
   there without a screen.
3. Otherwise a screen is made: the first `--screen-reads` reads of each read file (default
   100,000) are mapped at once, with the platform preset and `--secondary=no`, to up to three
   representatives of each option (taken from different nodes). Each genome is scored by the
   matched bases of the primary alignments on it, and the walk follows the subtree of the best
   genome. If the best genome matches less than 10 % of the screened bases, the walk stops
   there; if this happens at the first step, no local reference is chosen.
4. The walk ends at a node without further options; a representative of that node (or, when
   it has none, of its subtree) is the local reference.

The rounds also end when the same genome is chosen again, when no local reference is found, or
when the chosen genome does not fit its place in the tree (below).

### Lifting the canSNP positions

Each reference genome of the scheme is aligned to the local genome with the minimap2 assembly
command used for assembly typing (`asm20`, same options and mapping-quality rule). A canSNP is
lifted only when the local genome has its ancestral or derived base at the aligned position; a
gap, another base, or a position covered by two alignments with different bases is not lifted.
The lift is cached in `<ref-dir>/local/lift/<local>__<reference>.lift.tsv` (or in the run
directory when that folder cannot be written) and made again when either genome file, the
canSNPs, the minimap2 version or the options change.

The local genome is used only if it carries the derived base at more of the lifted canSNPs on
the path from the root to its own node than the ancestral base; otherwise a warning is logged
and the round is not made. This protects against a genome placed at the wrong node of the
scheme.

### Mapping and merging the counts

The reads are mapped to the local genome with the same commands and thresholds as before, and
counted at the lifted positions (complemented for canSNPs lifted to the reverse strand). Per
canSNP, the counts are merged with the earlier ones:

- the new counts are used when the canSNP was lifted and they reach `--min-depth`;
- otherwise the earlier counts are kept;
- when both reach `--min-depth` and give different decided states (derived, ancestral or other
  base), the canSNP is labeled `conflict` and counts as no data.

The `source` and `lifted_position` columns of `<sample>_site_counts.tsv` give the genome each
count came from and the lifted position. The `Reads:` line of `<sample>_snps.txt` names the
local reference, its node, its share of the screened bases at the last screen, the number of
screens, the number of rounds and the number of conflicts (with `--fetch-local` also the number
of genomes downloaded for the sample).

### Close screens

After the rounds, each screen step in which the second subtree had at least 20 % of the matched
bases and at least half the share of the first is checked against the final counts. If the
canSNPs of the two subtrees are not one derived and the other ancestral, a warning says that
the sample may hold more than one strain or lie between the two subtrees. Sibling genomes share
most of their sequence, so a close screen alone is expected between them and is not reported
when the counts separate them. This warning does not change the call or the call status.
