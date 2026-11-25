# Input formats

This page describes every file that `cansnper2 type` and the `cansnper2 db` verbs read: query
files, list files and samplesheets, the node-labels file, the scheme files, the modification file
of `db modify`, and the reference genome files. It is a reference for users who prepare inputs by
hand or in a pipeline. The options themselves are listed in [cli.md](cli.md).

## Query files

A query is a file given to `cansnper2 type`. CanSNPer2 tells the kind of a query from its content,
not from its file name. It reads the first non-blank characters of the file (after a UTF-8 byte
order mark, and after decompression when the file is gzip-compressed):

| First character | Kind |
|---|---|
| `>` | FASTA assembly |
| `@` | FASTQ reads |
| anything else | a list file (see [List files](#list-files)) |

Assemblies and reads cannot be typed in one run, because they are read with different settings.
A run that is given both stops with the message `reads and assemblies cannot be typed in one run`
and lists the files that are not reads.

### FASTA assemblies

An assembly is one FASTA file per sample, with one or more sequences (complete genomes, or
contigs of a draft assembly). The order of the sequences does not matter.

- The file must not be gzip-compressed; decompress it first.
- Each sequence is identified by the first word of its header line. These names must be set
  and unique within the file.
- The file must hold at least one header and one sequence line.
- Characters other than `A`, `C`, `G`, `T` (gaps, IUPAC ambiguity codes) are replaced by `N`, and
  blank lines and whitespace inside sequence lines are removed, on a temporary copy before
  alignment. [method.md](method.md#preparing-the-inputs) describes this step.

### FASTQ reads

Reads are FASTQ files, plain or gzip-compressed. Read typing needs minimap2 and samtools in
`PATH`. The sequencing platform is given with `--platform` (`illumina`, `ont` or `hifi`; default
`illumina`) or per sample in a samplesheet.

A sample is one read file (single-end or long reads) or two (paired-end reads). Without a
samplesheet, CanSNPer2 groups the files into samples by their names:

1. The endings `.fastq`, `.fq` and `.gz` are removed (in any letter case).
2. For `--platform illumina`, a name that ends in `_R1` or `_R2`, optionally followed by
   `_001`, or in `_1` or `_2`, is one mate of a pair. The two mates with the same remaining
   name form one sample of that name.
3. Any other file is a single-end sample named after the file without its endings.

Long reads (`--platform ont` or `hifi`) are never paired: each file is one sample, and a mate
suffix stays part of its name.

| Files | Samples |
|---|---|
| `S1_R1_001.fastq.gz`, `S1_R2_001.fastq.gz` | `S1` (paired) |
| `S2_1.fq.gz`, `S2_2.fq.gz` | `S2` (paired) |
| `S3.fastq.gz` | `S3` (single-end) |
| `S4_R1.fastq.gz` alone | refused |
| `S5_R1.fastq.gz`, `S5_1.fastq.gz` | refused |

A file named as one mate whose other mate is not given is refused, so that half of a pair is
never typed as single-end reads unnoticed. To type such a file as single-end reads, list it in a
samplesheet. Two files that give the same sample name and are not the two mates of one pair are
also refused; rename one of them or use a samplesheet.

### Samplesheet

`--samplesheet FILE` gives the read samples explicitly and replaces the query arguments; a run
cannot be given both.

- Tab-separated text. The first non-blank line is the header, which must begin with `sample` and
  `reads1`. The columns `reads2` and `platform` are optional, and other columns are not read.
- One sample per row. `sample` is the sample name; `reads1` is the read file of a single-end
  sample or the first mate of a pair; `reads2` the second mate.
- `platform` is `illumina`, `ont` or `hifi`. An empty or missing value takes `--platform`.
  A samplesheet may hold samples of several platforms.
- Relative paths are read from the folder of the samplesheet, not from the current folder.
- Blank lines are skipped and spaces around the values are removed.
- A sample name may be listed only once, may contain only letters, digits, `.`, `_` and `-`,
  and may not start with `.`. Every read file must exist. These checks are made before anything
  is typed.

```tsv
sample	reads1	reads2	platform
isolate7	reads/isolate7_R1.fastq.gz	reads/isolate7_R2.fastq.gz	illumina
isolate8	reads/isolate8.fastq.gz		ont
```

```sh
cansnper2 type --samplesheet samples.tsv -d francisella_tularensis.db
```

### List files

A list file names the queries of a run, one file path per line. It is recognized by its content:
a file that exists, is not empty, is not gzip-compressed, and whose first non-blank character is
neither `>` nor `@`.

- A list file must be the only query argument.
- Blank lines are skipped; every other line is a path. Spaces at the start and end of a line are
  removed.
- Relative paths are read from the current folder, not from the folder of the list file.
- The paths may name assemblies or read files; read files are grouped into samples by the rules
  above.
- A FASTA file without its header line would be taken as a list file. When the first lines look
  like sequence and none of them names an existing file, the run stops with a message that asks
  for a header line.

```sh
cansnper2 type queries.txt -d francisella_tularensis.db
```

### Sample names

Every result file, alignment file and row of the summary is named after the sample.

- An assembly is named after its file name without the last extension:
  `GCA_000018925.1.fasta` gives `GCA_000018925.1`, and `isolate.contigs.fna` gives
  `isolate.contigs`.
- A read sample is named as described above, or by the `sample` column of a samplesheet.

Two different files that give the same sample name (for example `a/x.fasta` and `b/x.fasta`)
stop the run before anything is aligned; the same file given twice is typed once. A file name
that gives an empty sample name is refused. Names taken from file names are not otherwise
restricted; names given in a samplesheet follow the rule above. Since the name becomes part of
file names, letters, digits, `.`, `_` and `-` are the safe choice for file names as well.

## Node-labels file

`--node-labels FILE` adds labels to the tree drawings in `summary.html`, `<sample>_tree.html`
and `summary_tree.html` (see [output-formats.md](output-formats.md)). It does not change the
call.

- Tab-separated text, one label per line: `<node>`, `<kind>`, `<rank>`, `<label>`, and an
  optional `<color>`.
- `node` is a node name of the database, such as `B.1`. A label whose node is not in the
  database is skipped with one warning that names it.
- `kind` is `node`, which draws a badge under the node, or `clade`, which draws a dashed box
  around the subtree below the node.
- `rank` is free text (for example `subspecies` or `marker`). It is carried by the label and is
  not taken from the depth of the node in the tree. Labels without a color are colored by rank.
- `color` is a color for the drawing, such as `#7b3294`.
- Blank lines and lines starting with `#` are skipped. A first line whose first field is `snp`,
  `snp_id`, `node` or `name` is read as a header and skipped.

```tsv
node	kind	rank	label	color
B.1	clade	subspecies	holarctica
M.1	clade	subspecies	mediasiatica	#7b3294
A.I.14	node	clade	A.I.14 group
```

## Scheme files

A scheme is a set of tab-separated text files. `cansnper2 db create` and `db import` build a
database from them, and `cansnper2 scheme get` fetches them from CanSNPer2-data and builds the
database itself. The full format of the tree, canSNP and references files is defined in the
[`schemes/README.md`](https://github.com/FOI-Bioinformatics/CanSNPer2-data/blob/master/schemes/README.md)
of CanSNPer2-data; this section summarizes them, and specifies the genomes file and the checks
that CanSNPer2 applies. Building and maintaining databases is described in
[schemes-and-databases.md](schemes-and-databases.md).

| File | Option | Content |
|---|---|---|
| `<name>_tree.txt` | `--tree` | One line per path from the root to a node, with the node names in tab-separated columns. |
| `<name>_canSNPs.txt` | `--cansnps` | Every canSNP of every node: node, rank within the node, publication, reference genome, reference sequence and 1-based position, derived and ancestral base, score, and whether the scheme types it (`typed` is `yes` for exactly one canSNP per node). |
| `<name>_references.txt` | `--references` | One row per reference genome: short genome name, strain, GenBank and RefSeq assembly accession, assembly name. |
| `<name>_genomes.txt` | `--genomes` | Optional: genomes placed at nodes, with their role (see below). |

A node can have several canSNPs, ranked from 1 within the node. Exactly one of them is marked
`typed` `yes`, and CanSNPer2 calls only that one; the others are alternatives that the database
keeps and that `cansnper2 db export` writes back.

A canSNP file of CanSNPer2 2.3.0 (one canSNP per node, with an `snp_id` column) is still read:
each row becomes the typed canSNP of rank 1 of its node. Files from before 2.3.0, with a single
`position` column along the joined reference sequences, are refused; convert them with
`cansnper2 db convert` (see [schemes-and-databases.md](schemes-and-databases.md)).

### Genomes file

The genomes file places genomes at nodes of the tree. CanSNPer2 stores it in the database, and
read typing with `--local-references` takes its candidate genomes from the rows with the role
`representative`. Databases built before 2.6.0 hold no genomes.

- Tab-separated text. The first line is the header and must begin with `genome`, `node` and
  `role`; a leading `#` is allowed. Further columns (CanSNPer2-data uses `basis` and `source`) are
  allowed and are not stored.
- One genome per row; blank lines are skipped.

| Column | Content |
|---|---|
| `genome` | The genome. Representative genomes are downloaded by this name as an NCBI assembly accession (for example `GCF_000008985.1`) and stored as `<ref-dir>/local/<genome>.fna`. |
| `node` | A node of the tree. Empty for an outgroup, and only for an outgroup. |
| `role` | `representative` (aligned when the scheme was made; used as a local reference), `member` (placed at the node through a representative) or `outgroup` (aligned, but placed at no node). |

The file is refused, with every problem listed, when a row has no genome or no role, a role is
not one of the three, a genome is placed at a node that is not in the tree, a genome other than
an outgroup has no node, an outgroup has a node, or a genome is listed twice.

```tsv
genome	node	role
GCF_000008985.1	A.I.13	representative
GCA_009813935.1	B.53	member
```

`cansnper2 db export` writes the genomes of a database back as `<database>_genomes.txt`.

### Checks before a database is built

`db create`, `db import` and `scheme get` check the tree, canSNP and references files before
anything is written. All problems are reported together, each with its file and line, under
`The source files cannot be built into a database, N problem(s):`, and no database is written.
`db modify` applies the same row checks to the files it adds (see
[Modification file](#modification-file-db-modify)).

| File | Refused when |
|---|---|
| any | The text is not valid UTF-8. (A UTF-8 byte order mark is removed, a file ending in `.gz` is read compressed, and lines ending in CR LF are read.) |
| tree | The file is empty; a line has an empty column; a line starts at another root than the first line; the root is given a parent; a node is given more than one parent. Columns may be separated by tabs or `;`. |
| canSNPs | The file is empty or has a header but no rows; a column is named twice; a required column is missing; the file has a `position` column instead of `sequence_position`; a row has too few or too many columns. |
| canSNPs, per row | `rank` is not a positive integer; the ranks of a node are not 1, 2, 3 ... in file order; `score` is neither empty nor a number; `typed` is not `yes` or `no`; a node of the file has no canSNP marked `typed` `yes`, or more than one; `reference` is empty; `node` is not in the tree; `genome` is not in the references file; `sequence` is empty; `sequence_position` is not a positive integer or is larger than 2147483647; a base is not `A`, `C`, `G` or `T`; derived and ancestral base are equal; one position on one sequence of one genome is annotated twice. |
| references | The file is empty; a column is named twice or missing; a row does not have as many columns as the header; a genome is listed twice; `genome` contains characters other than letters, digits, `.`, `_` and `-`, or starts with `.`; `genbank_id` is not an assembly accession of the form `GCA_` or `GCF_`, nine digits and a version (such as `GCA_000008985.1`); two genomes share a `genbank_id`; a genome name or accession repeats the name or accession of another genome. |

In the canSNP file, lines that start with `#` are skipped, and trailing columns named `comment`,
`comments`, `note`, `notes`, `free text` or `free_text` may be left out of a row. CanSNPer2-data
applies further rules to the files it publishes (LF line endings, no underscore in genome names,
a README per scheme); its `schemes/README.md` lists them.

## Modification file (db modify)

`cansnper2 db modify --mod-file FILE --parent NODE` adds nodes to the tree of a database or
replaces part of it. `--parent` is required with `--mod-file` and must name a node of the
database.

- Tab-separated text. The first line is the header: `parent` and `child`, in this order or as
  `child` and `parent` (the case of the first name does not matter). A third column, a level
  name, is optional.
- One link per line: a parent node and a child node, and optionally a level. Blank lines are
  skipped; a line with fewer than two or more than three columns is refused.
- Nodes that are already in the database are linked by their name; other names become new
  nodes.
- Without `--replace`, the links are added below `--parent`. With `--replace`, the nodes below
  `--parent` that the file does not give again are removed; a file without rows is then refused,
  as it would remove the whole branch.

```tsv
parent	child
B.136	B.136.x1
B.136.x1	B.136.x2
```

```sh
cansnper2 db modify -d mod.db --mod-file mod.txt --parent B.136
```

New nodes have no canSNPs. Give their canSNPs with `--cansnps` (and new reference genomes with
`--references`), in the same command or a later one; these files have the format of the scheme
files above and are checked against the tree and the reference genomes of the database. The
change is made on a copy of the database, which replaces it only when the result is still a valid
tree and every check passed; otherwise the database is left unchanged. A modified database
records the status `unknown` and the date of the change (see
[reproducibility.md](reproducibility.md)).

## Reference genome files

`cansnper2 type` aligns each query against the reference genomes of the database. They are read
from the folder given with `--ref-dir` (default `references/`, relative to `--workdir`).

- Each reference genome is one uncompressed FASTA file named `<genome>.fna`, where `<genome>` is
  the `genome` column of the references file. `cansnper2 ref get` and `cansnper2 scheme get` write
  these files (as links into `<ref-dir>/source/`).
- A reference genome may hold several sequences (chromosomes, plasmids). Each sequence is
  identified by the first word of its header line, and each canSNP by that name and its 1-based
  position on the sequence, so the order of the sequences does not matter and sequences without
  canSNPs are allowed. See [method.md](method.md#references-with-more-than-one-sequence).
- The sequence names must be the names the canSNP file gives in its `sequence` column. The
  published schemes use GenBank sequence accessions (for example `AJ749949.2`). The RefSeq
  assembly of the same genome names its sequences differently (`NC_`, `NZ_`); `cansnper2 ref get
  -s refseq` renames them to the GenBank accessions after download. A file with other names is
  refused.
- `.fna` files in the folder that are not reference genomes of the database are skipped, with one
  warning that names them.
- Representative genomes for `--local-references` are read from `<ref-dir>/local/<genome>.fna`.

Before the first sample is aligned, each reference file is checked against the canSNPs that the
database places on it. The run stops when a sequence name is empty or used twice, when a sequence
named by the database is missing, when a sequence is shorter than its largest canSNP position, or
when at a fifth or more of its canSNPs the reference base is neither the ancestral nor the derived
base. Fewer such mismatches give a warning. A reference genome of the database that carries
canSNPs and has no file stops the run, unless `--allow-missing-references` is given.
