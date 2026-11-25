# Schemes and databases

This page is for users who fetch published typing schemes and for users who build, change or
convert their own databases. It describes the `cansnper2 scheme`, `ref` and `db` commands. The
formats of the scheme files are described in [input-formats.md](input-formats.md), and every
option is listed in [cli.md](cli.md).

A scheme is a set of text files: a tree, a canSNP file, a reference genomes file and, optionally,
a genomes file. `cansnper2 type` does not read these files directly. It reads a database, an
SQLite file built from them.

## Published schemes

### List the schemes

The published schemes are kept in
[CanSNPer2-data](https://github.com/FOI-Bioinformatics/CanSNPer2-data). `scheme list` prints
the schemes of a branch with their status:

```sh
cansnper2 scheme list
cansnper2 scheme list --branch dev
```

```text
SNP schemes in master of FOI-Bioinformatics/CanSNPer2-data:
  bacillus_anthracis             verified
  brucella_classical             verified
  coxiella_burnetii              verified
  francisella_tularensis         verified
  yersinia_pestis                verified
```

The branch `master` holds verified schemes only. The branch `dev` holds the same schemes and, in
addition, draft schemes that are being evaluated (on 25 November 2026: `brucella_2017`,
`brucella_anthropi`, `brucella_genus` and `brucella_suis`). The status is read from the Status
row of the README of each scheme; a scheme without one is listed as `status not given`.
`scheme list` writes no log file.

### Fetch a scheme and build its database

```sh
cansnper2 scheme get francisella_tularensis
cansnper2 scheme get brucella_suis --branch dev
```

`scheme get` does the following:

1. It resolves the branch to a commit through the GitHub API and reads every file at that
   commit, so that a change pushed during the download cannot mix two versions of the scheme.
2. It downloads the canSNP, tree and reference genomes files, and the genomes file when the
   scheme has one.
3. It compares the SHA-256 of each file with the checksum recorded in the manifest of the
   repository (`database/<name>.manifest.json`). A file that differs stops the command.
4. It checks the files and builds the database with the installed version of CanSNPer2 (the
   checks are those of `db create`, see [Build a database](#build-a-database)). The database is
   written under a temporary name and renamed when it is complete, so a failed build leaves no
   file behind.
5. It writes the provenance file `<name>.source.json` next to the database, where `<name>` is
   the database file name without `.db` (for `francisella_tularensis.db`,
   `francisella_tularensis.source.json`).
6. It downloads the reference genomes from NCBI GenBank into `./references`, as `ref get` does
   (see [Reference genomes](#reference-genomes)).

The database is written to `<NAME>.db` unless `-d` names another file. An existing database is
not replaced unless `--overwrite` is given. Further options:

| Option | Effect |
| --- | --- |
| `--branch REF` | read the scheme from another branch, a tag or a commit; `--branch dev` is needed for draft schemes |
| `--no-references` | build the database without downloading the reference genomes |
| `--local-references` | also download the representative genomes of the scheme into `references/local`, for read typing with `cansnper2 type --local-references` |
| `--allow-unverified` | build a verified scheme although its files cannot be checked against a manifest |

The reference genomes are always downloaded from GenBank into `./references`. To download them
from RefSeq or into another folder, give `--no-references` and run `ref get` afterwards.

### Verified and draft schemes

| Situation | Verified scheme | Draft scheme |
| --- | --- | --- |
| a file differs from the checksum in the manifest | refused | refused |
| no manifest, or the manifest has no checksum for a file | refused unless `--allow-unverified` is given | built with a warning |
| the scheme is built | no warning | a warning that the scheme is a draft |

When the `schemes.toml` file of the repository lists known issues for a scheme (canSNPs whose
calls may be unreliable), `scheme get` prints a warning that names them, and the database records
them.

### Provenance file

`<name>.source.json` records where the database came from:

| Key | Content |
| --- | --- |
| `scheme`, `status` | the scheme name and its status from the scheme README |
| `repository`, `ref`, `commit` | the repository URL (or the local directory), the branch asked for and the commit read |
| `sources` | the path of each scheme file in the repository and its SHA-256 |
| `source_type` | the format of the source files (`CanSNPer`) |
| `built`, `cansnper2` | the date of the build and the CanSNPer2 version that built it |
| `uncommitted_sources` | with `--data-dir` in a git checkout only: scheme files with uncommitted changes |

Because a commit can be given with `--branch`, a database can be built again from the same
files with the `commit` value of this file.

### Other repositories and local directories

Schemes can be read from any GitHub repository, or local directory, that is laid out as
CanSNPer2-data:

```sh
cansnper2 scheme list --data-dir my-schemes
cansnper2 scheme get my_scheme --data-dir my-schemes
cansnper2 scheme get my_scheme --repo owner/my-schemes
```

Only `schemes/<name>/<name>_canSNPs.txt`, `<name>_tree.txt` and `<name>_references.txt` are
required. A scheme that has a SNP file of 2.3.0 (`<name>_snp.txt`) instead of a canSNP file is
built from it with a warning. `<name>_genomes.txt`, the scheme `README.md` (with its Status row), `schemes.toml` and
`database/<name>.manifest.json` are used when present.

- `--repo OWNER/NAME` reads the default branch of that repository unless `--branch` is given.
- `--data-dir DIR` reads the directory as it is on disk and cannot be combined with `--repo` or
  `--branch`. A file that differs from its manifest, or a missing manifest, gives a warning
  instead of an error. When the directory is in a git checkout, the commit is recorded, and
  scheme files with uncommitted changes are listed in `source.json` and marked with
  `+uncommitted` after the commit in the database.
- A scheme whose canSNP file has a `position` column (positions along the joined reference
  sequences, the format before 2.3.0) is refused. A local file in that format can be converted
  with `db convert` (see [Convert databases from earlier versions](#convert-databases-from-earlier-versions)).

### GitHub API limits

Listing the schemes and resolving a branch to a commit use the GitHub API, which allows 60
anonymous requests per hour from one address. On shared machines, such as compute clusters or
continuous integration runners, this limit can be reached. Set the environment variable
`GITHUB_TOKEN` (or `GH_TOKEN`) to a GitHub token to raise it; the token is sent to
`api.github.com` only. When a branch cannot be resolved to a commit, `scheme get` reads the files
from the branch itself, prints a warning, and records the commit as `unknown` in the database.

## Reference genomes

### Download the reference genomes of a database

```sh
cansnper2 ref get -d francisella_tularensis.db
cansnper2 ref get -d francisella_tularensis.db -s refseq -o refseq_references/
```

`ref get` downloads, for each reference genome of the database, the assembly
`<accession>_<assembly name>_genomic.fna.gz` from the NCBI FTP site, using the `genbank_id`
(default, `-s genbank`) or the `refseq_id` (`-s refseq`) of the reference genomes file. The
files are unpacked into `<outdir>/source/`, and `<outdir>/<genome>.fna` is created as a link to
the unpacked file. `<genome>` is the genome name of the scheme, the name `cansnper2 type` looks
for in `--ref-dir`. The default folder is `references`.

```text
references/
  FSC200.fna -> source/GCA_000168775.2_ASM16877v2_genomic.fna
  LVS.fna -> source/GCA_000009245.1_ASM924v1_genomic.fna
  ...
  source/
```

Each downloaded genome is checked against the canSNPs the database places on it: every sequence
named by a canSNP must be in the file, each sequence must reach the largest canSNP position on
it, and the base at each canSNP position must be the ancestral or the derived base. A genome
that lacks this base at a fifth of its canSNPs or more is refused; fewer mismatches give a
warning. The command ends with exit status 1 and lists the genomes that could not be fetched or
did not pass the check.

A genome that is already unpacked in `<outdir>/source/` and passes the check is not downloaded
again, so `ref get` can be run again after an interrupted download. A file `<genome>.fna` that
is not a link (for example a genome copied there by hand) is checked in the same way and, when it
passes, used instead of a download; when it does not pass, `ref get` keeps it, reports the
genome as not fetched and names the problem.

### Sequence names, GenBank and RefSeq

A canSNP is placed on a sequence by the first word of the FASTA header of the reference genome,
which in the published schemes is the GenBank accession of the sequence (for example
`AJ749949.2`). The order of the sequences in the file does not matter, and sequences without
canSNPs, such as plasmids, are allowed.

GenBank assemblies carry these names. In a RefSeq assembly the sequences are named by RefSeq
accessions (`NC_`, `NZ_`), so `ref get -s refseq` downloads the NCBI assembly report of each
genome and renames each sequence to its GenBank accession. A sequence that the report marks as
different from its GenBank sequence (Relationship `<>`) keeps its RefSeq name, with a warning;
if a canSNP lies on such a sequence, the reference check refuses the genome.

### Representative genomes for read typing

```sh
cansnper2 ref get -d francisella_tularensis.db --local-references
```

`--local-references` also downloads the representative genomes listed in the genomes file of the
scheme into `<outdir>/local/<accession>.fna`. A `GCF_` accession is downloaded from RefSeq and a
`GCA_` accession from GenBank, whatever `-s` says. These genomes are used by
`cansnper2 type --local-references`, which is described in [method.md](method.md). A database
without representative genomes gives an error.

## Build a database

### Create a database from scheme files

```sh
cansnper2 db create -d my_scheme.db --tree my_scheme_tree.txt \
    --cansnps my_scheme_canSNPs.txt --references my_scheme_references.txt
cansnper2 ref get -d my_scheme.db
cansnper2 type sample1.fna sample2.fna -d my_scheme.db
```

`--genomes my_scheme_genomes.txt` adds the genomes file, which is needed only for read typing
with `--local-references`. The formats of the four files are described in
[input-formats.md](input-formats.md).

Before anything is written, the files are checked with the rules of the `validate.py` script of
CanSNPer2-data. Among other things, the tree must have one root and each node one parent; each
canSNP must lie on a node of the tree and on a genome of the reference genomes file, have a
positive position and two different bases out of A, C, G and T, and the canSNPs of a node must be
ranked 1, 2, 3 ... with exactly one marked `typed yes`; each reference genome needs a GenBank or
RefSeq assembly accession with its version (such as `GCA_000008985.1`). Every problem is listed
with its file and line, and no database is written when there is one.

`db create` writes a new file only: an existing database is refused, and so is a database path
in a folder that does not exist. A build that fails part way removes the file it had started.
`db create` also accepts only some of the files, for example the tree alone, and the database can
then be completed with `db modify` (see [Change a database](#change-a-database)); `db import`
needs all three.

### Scheme information

The options `--scheme-name`, `--scheme-status` (`verified`, `draft` or `unknown`) and
`--scheme-source` record which scheme the database holds. They are accepted by `db create` and
`db import` only.

```sh
cansnper2 db create -d my_scheme.db --tree my_scheme_tree.txt \
    --cansnps my_scheme_canSNPs.txt --references my_scheme_references.txt \
    --scheme-name my_scheme --scheme-status draft --scheme-source "my-schemes@1a2b3c4"
```

The information is kept in the `scheme_info` table of the database. `db info` prints it,
followed by the size of the database:

```sh
cansnper2 db info -d francisella_tularensis.db
```

```text
scheme	francisella_tularensis
status	verified
source	FOI-Bioinformatics/CanSNPer2-data
ref	master
commit	d7493fd26acfe08ff57beeed1e00ebcce0607cc3
known_issues	none
format	2.4
built_with	2.6.0
built	2026-11-25

nodes	328
canSNPs	328
typed canSNPs	328
canSNPs per node	1
reference genomes	5
```

| Key | Content |
| --- | --- |
| `scheme` | the name of the scheme |
| `status` | `verified`, `draft` or `unknown` |
| `source` | the repository (`owner/name`) or local directory the scheme files were read from, or the text of `--scheme-source` |
| `ref` | the branch, tag or commit asked for (`scheme get` from GitHub only) |
| `commit` | the commit the files were read at; `unknown` when it could not be resolved; followed by `+uncommitted` (a local checkout with uncommitted scheme files) or `+modified` (changed with `db modify` or `db remove`, see [Change a database](#change-a-database)); `modified` alone when no commit was recorded before the change |
| `known_issues` | canSNPs that `schemes.toml` lists as known issues, comma-separated; `none` when there are none |
| `format` | the database format (`2.4`) |
| `built_with` | the CanSNPer2 version that built the database, or that made the last change with `db modify` or `db remove` |
| `built` | the date of the build; a later change does not alter it and is recorded in `modified` |
| `modified` | the date of the last change with `db modify` or `db remove` |

`format`, `built_with` and `built` are always written. For a database built with `db create` or
`db import`, `scheme`, `status` and `source` are the values given with `--scheme-name`,
`--scheme-status` and `--scheme-source`; record them when you build a scheme of your own. Keys
without a value are not printed, and an empty value is printed as `none`. A database
without the table (built before 2.3.0, or by hand) is reported with the line
`the database records no scheme information (built with CanSNPer2 before 2.3.0 or by hand); check where it came from`.
The counts are the nodes of the tree, the canSNPs and how many of them are typed, the largest
number of canSNPs on one node, and the reference genomes.

When `cansnper2 type`, `ref get` or `db export` open a database, they report its scheme once. A
draft or unknown status, known issues, a changed database or a database without scheme
information give a warning; a verified scheme is reported only with `--verbose`. The scheme
information is never a reason to stop a command.

## Export and import

### Export a database

```sh
cansnper2 db export -d francisella_tularensis.db -o export/
```

The default format, `source`, writes the files that `db create` and `db import` read, named
after the database:

| File | Content |
| --- | --- |
| `<db>_tree.txt` | one line per node, the tab-separated path from the root to the node |
| `<db>_canSNPs.txt` | every canSNP of every node, ranked within its node, with the header of CanSNPer2-data |
| `<db>_references.txt` | the reference genomes |
| `<db>_genomes.txt` | the genomes of the scheme, only when the database has them |

The command prints the `db import` command line that rebuilds the database from these files.
The date column of the canSNPs is not exported; a database built from the files records its own
build date. The scheme information is not exported either, so it is given again with the
`--scheme-*` options of `db import`. Existing files are not overwritten: the export stops and
writes nothing. A database whose tree has more than one root, a node with more than one parent,
or canSNPs on a genome that is not among its reference genomes cannot be written in this format
and is refused. The format name `tab` is accepted as another name for `source`.

`--format newick` writes the tree only, to `<outdir>/<db>_tree.nwk` with `-o`, or to standard
output without it:

```sh
cansnper2 db export -d francisella_tularensis.db --format newick -o export/
cansnper2 db export -d francisella_tularensis.db --format newick > francisella_tularensis.nwk
```

### Import the files again

`db import` builds a new database from the exported files. It runs the same checks as
`db create`, needs all three of `--tree`, `--cansnps` and `--references`, and in addition checks
that the result forms a single rooted tree before it keeps the database. A scheme can therefore
be exported, edited as text and imported:

```sh
cansnper2 db export -d my_scheme.db -o scheme/
cansnper2 db import -d my_scheme_edited.db --tree scheme/my_scheme_tree.txt \
    --cansnps scheme/my_scheme_canSNPs.txt --references scheme/my_scheme_references.txt \
    --genomes scheme/my_scheme_genomes.txt --scheme-name my_scheme --scheme-status unknown
```

The rebuilt database has the same tree, canSNPs and reference genomes as the exported one. The
destination must not exist.

## Change a database

`db modify` and `db remove` change a database in place. Each change is made on a copy, which
replaces the database only when every step succeeded and the result is still a valid tree; a
failed change leaves the database unchanged. A successful change marks the scheme information:
the status becomes `unknown`, `+modified` is added to the commit, and the date and CanSNPer2
version are recorded in `modified` and `built_with`. `cansnper2 type` then warns that the
database no longer matches the scheme files.

### Add nodes to the tree

`--mod-file` gives tree links to add below the node named by `--parent`. The file has a header
row `parent` and `child` (an optional third column names a level) and one tab-separated link per
row:

```text
parent	child
B.103	B.103.1
B.103.1	B.103.1.1
```

```sh
cansnper2 db modify -d my_scheme.db --mod-file new_branch.txt --parent B.103
```

With `--replace`, the nodes below `--parent` that the file does not list are removed, together
with their canSNPs:

```sh
cansnper2 db modify -d my_scheme.db --mod-file subtree.txt --parent A.II.2 --replace
```

`--parent` is required with `--mod-file`. A change that would give a node a second parent, create
a cycle or link a node to the root is refused.

### Add canSNPs and reference genomes

```sh
cansnper2 db modify -d my_scheme.db --cansnps added_canSNPs.txt
cansnper2 db modify -d my_scheme.db --references added_references.txt
cansnper2 db modify -d my_scheme.db --mod-file new_branch.txt --parent B.103 \
    --cansnps added_canSNPs.txt --references added_references.txt
```

The files have the formats of the scheme files. They add rows; they do not change rows the
database already holds. A canSNP on a node that is neither in the database nor in the
modification file, on a genome that is not listed, at a rank its node already has, a second typed
canSNP of a node, a position the database already holds, or a genome or accession that is
already in the database is refused, and all problems are listed. To change existing canSNPs or
reference genomes, export the database, edit the files and import them.

### Remove a branch

```sh
cansnper2 db remove -d my_scheme.db --parent B.103
```

`db remove` removes the node given with `--parent`, every node below it and their canSNPs. The
root of the tree cannot be removed.

## Convert databases from earlier versions

Databases and canSNP files written by CanSNPer2 before 2.4.0 are refused by the current version
with a message that names `cansnper2 db convert`. `db convert` writes a converted copy into the
folder given with `-o`; the input is not changed, and an existing output file is not
overwritten.

```sh
cansnper2 db convert -d old.db -o converted/
cansnper2 db convert -d old.db --ref-dir references/ -o converted/
cansnper2 db convert --cansnps my_scheme_snp.txt -o converted/
```

| Input | What changes | Reference files |
| --- | --- | --- |
| database or SNP file of 2.3.0 (one canSNP per node) | each canSNP becomes rank 1 of its node, typed, without a score | not needed |
| database or SNP file from before 2.3.0 (positions along the joined reference sequences) | positions are moved to the sequences of the reference genomes, then as above | needed, `--ref-dir` |

For formats before 2.3.0, each canSNP is checked against its reference genome
`<ref-dir>/<genome>.fna`, and nothing is written when a base matches neither allele. The
sequences of the reference file must be in the order in which the old positions were counted
(the order of the NCBI assembly file). For a database, missing reference files are downloaded
from GenBank with the accessions the database records; for a SNP file, every reference file must
already be in `--ref-dir`. A converted SNP file is written as the canSNP file of its scheme
(`<name>_snp.txt` becomes `<name>_canSNPs.txt`), and a converted database records the status
`unknown` and the converted file as its source.

The formats are described in more detail in [upgrading.md](upgrading.md).
