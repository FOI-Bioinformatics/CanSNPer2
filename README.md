# CanSNPer2
CanSNPer2: A toolkit for SNP-typing NGS data.

CanSNPer2 types assemblies, or sequencing reads directly (see [Read typing](#read-typing)).

Databases supplied can be found at https://github.com/FOI-Bioinformatics/CanSNPer2-data

* Francisella tularensis
* Bacillus anthracis
* Yersinia pestis
 
#### Upcoming databases
* Brucella
* Coxiella


## Installation

bioconda (recommended, conda works too but mamba is faster):

```
mamba create -n cansnper2 -c conda-forge -c bioconda cansnper2
mamba activate cansnper2
```

The package installs one command: `cansnper2` (see [The cansnper2
command](#the-cansnper2-command-since-260)). Bioconda builds before 2.6.0 installed
the separate commands `CanSNPer2`, `CanSNPer2-database`, `CanSNPer2-download` and
`CanSNPer2-test`; those are replaced and no longer installed. The recipe for the
package is in [conda/meta.yaml](conda/meta.yaml).

From source, in any Python 3.11+ environment:

```
python -m pip install .        # pipx install . works too
```

## Requirements (for manual install conda will install all dependencies)

* progressiveMauve (the default aligner, see [Aligners](#aligners))
* minimap2 (optional, `--aligner minimap2`)

Read typing
* minimap2 and samtools (see [Read typing](#read-typing))

## The cansnper2 command (since 2.6.0)

A single command with a verb for every job. Since 2.6.0 it is the only command that is
installed: `CanSNPer2`, `CanSNPer2-database`, `CanSNPer2-download` and `CanSNPer2-test`
are replaced by the verbs and their options. `cansnper2 --version` prints the version,
`--help` documents any verb. If the old scripts of a previous installation are still in
your PATH, remove them.

```sh
cansnper2 type query.fna -d scheme.db          # type assemblies (was CanSNPer2)
cansnper2 selftest                             # check the installation (was CanSNPer2-test)
cansnper2 db create -d new.db --tree tree.txt --cansnps snps.txt --references refs.txt
cansnper2 db import -d new.db --tree tree.txt --cansnps snps.txt --references refs.txt
cansnper2 db export -d scheme.db --format source --outdir export/    # or newick, tab
cansnper2 db info -d scheme.db
cansnper2 db convert -d old.db --ref-dir references/ --outdir converted/
cansnper2 db modify -d scheme.db --mod-file changes.txt --parent NODE --replace
cansnper2 db remove -d scheme.db --parent NODE
cansnper2 scheme list                          # schemes of CanSNPer2-data (was --list_schemes)
cansnper2 scheme get francisella_tularensis    # build its database and fetch references
cansnper2 ref get -d scheme.db --source genbank --outdir references/  # was CanSNPer2-download
```

`db create` and `db import` also take `--scheme-name`, `--scheme-status` and
`--scheme-source` to record what scheme the database holds; `db modify` also takes
`--cansnps` and `--references` for changes to the annotations and the reference list;
`db convert` takes `--cansnps` instead of `-d` to convert a SNP source file. `scheme get`
and `scheme list` take `--repo owner/name`, `--branch ref` and `--data-dir DIR` for
another source repository, and `scheme get` takes `--overwrite`, `--no-references` and
`--allow-unverified`.

Every verb has the same logging options: `--log-dir DIR` for the log files, `-v` or
`--verbose` for process information, `--debug` for debug output and the traceback of an
unexpected error, and `--quiet` to suppress warnings (the old `--supress`).

The options of `type` are the arguments of CanSNPer2 with one consistent spelling:
`--database` replaces `-db`, and underscores became hyphens (`--min-required-hits`,
`--alignment-timeout`, `--aligner-path`, `--allow-missing-references`,
`--keep-going`, `--keep-temp`, `--skip-alignment`, `--ref-dir`, `--tmp-dir`,
`--log-dir`); `--supress` became `--quiet`, `--annotations` became `--node-labels`,
and `--no-snp-files` the counterpart of the old `--no_snpfiles`. The rest are spelled
as before: `--outdir`, `--save-tree`, `--summary`, `--workdir`, `--strictness`,
`--mask`, `--rerun`, `--aligner`, `--verbose`. `--read_input` (never implemented) and
the hidden `--no_save_tree` are not offered by the verb.

The sections of this guide below still show the old spelling while they are being
migrated; they map to the verbs as follows, with underscores in option names replaced
by hyphens:

| old command | cansnper2 verb |
| --- | --- |
| `CanSNPer2 ...` | `cansnper2 type ...` |
| `CanSNPer2-test` | `cansnper2 selftest` |
| `CanSNPer2-database --create` | `cansnper2 db create` |
| `CanSNPer2-database --import` | `cansnper2 db import` |
| `CanSNPer2-database --export` | `cansnper2 db export` (`--export_format` is `--format`) |
| `CanSNPer2-database --info` | `cansnper2 db info` |
| `CanSNPer2-database --convert` | `cansnper2 db convert` |
| `CanSNPer2-database --mod_file ...` | `cansnper2 db modify ...` |
| `CanSNPer2-database --remove --parent` | `cansnper2 db remove --parent` |
| `CanSNPer2-download --list_schemes` | `cansnper2 scheme list` |
| `CanSNPer2-download --scheme NAME` | `cansnper2 scheme get NAME` |
| `CanSNPer2-download --database DB` | `cansnper2 ref get -d DB` |

## User guide CanSNPer2 (for custom databases see below)

1. Build a database from a scheme in [CanSNPer2-data](https://github.com/FOI-Bioinformatics/CanSNPer2-data) and download its reference genomes

```sh
cansnper2 scheme list                       # schemes on master (verified)
cansnper2 scheme get francisella_tularensis # builds francisella_tularensis.db, downloads references/
```

The scheme files (SNPs, tree and references) are fetched and checked against the SHA-256
recorded in the manifest of CanSNPer2-data, and the database is built locally by the installed
CanSNPer2. `francisella_tularensis.source.json` beside the database records the repository
commit and the hash of each file it was built from. A verified scheme whose files cannot be
checked (no manifest, or a file the manifest has no hash for) is refused unless
`--allow-unverified` is given. `--branch dev` gives access to draft
schemes that have not yet been verified; building one prints a warning. `--branch` also takes a
commit, which makes a database reproducible from its `source.json`.

Custom schemes can be read the same way from a local directory or from another GitHub
repository, as long as it is laid out as CanSNPer2-data. Only
`schemes/<name>/<name>_{canSNPs,tree,references}.txt` is required; `schemes.toml`, the scheme
`README.md` (with its Status) and `database/<name>.manifest.json` are used when present.

```sh
cansnper2 scheme list --data-dir my-schemes
cansnper2 scheme get my_scheme --data-dir my-schemes       # a local folder or clone
cansnper2 scheme get my_scheme --repo owner/my-schemes     # default branch of that repository
```

A local directory is read as it is on disk: a file that differs from its manifest gives a
warning rather than an error, and when the directory is in a git checkout, `source.json`
records the commit and lists the scheme files with uncommitted changes.

Listing schemes and resolving a branch to a commit use the GitHub API, which allows 60
anonymous requests an hour per address. On shared machines, set `GITHUB_TOKEN` (or `GH_TOKEN`)
to a GitHub token to raise the limit; it is sent to the API only.

From 2.3.0 the database also records which scheme it holds, in a `scheme_info` table: the
scheme, its status (`verified`, `draft` or `unknown`), the repository or folder, the branch and
commit (`unknown` when the GitHub API could not resolve it), the SNPs that `schemes.toml` lists
as known issues, and the CanSNPer2 version and date of the build. `cansnper2 db info -d <db>` prints
the table. Typing, `cansnper2 ref get` and `cansnper2 db export` report
the scheme once when they open the database; a draft or unknown status, known issues, or a
database without the table (built before 2.3.0 or by hand) give a warning, and a verified scheme
is reported only with `--verbose`. The table is informative only and never stops a run.
A database changed with `db modify`, `db remove` or their `--cansnps` and `--references` files no longer
matches the scheme files it was built from, so the change sets its status to `unknown`, adds
`+modified` to the commit (or records `modified` when there was none), records the date in
`modified` and the CanSNPer2 version in `built_with`; typing then gives a warning. A database
converted with `db convert` records status `unknown` and the converted file as its source.
`db info` prints `none` for an empty value, such as `known_issues` of a scheme without any.

```
$ cansnper2 db info -d brucella_suis.db
scheme	brucella_suis
status	draft
source	FOI-Bioinformatics/CanSNPer2-data
ref	dev
commit	abc1234000000000000000000000000000000000
known_issues	B.1
format	2.4
built_with	2.4.0
built	2026-11-21
```

2. Alternatively, with a database you already have, download its references

```sh
cansnper2 ref get -d my_database.db
```

3. Run genomes
```sh
cansnper2 type fastadir/*.fasta -d francisella_tularensis.db --summary
```

The summary parameter gives `snps_summary.txt`, a table of the final SNP, its support and the SNP path of every sample, and a tree (self-contained html, opens in any browser) with the final SNPs of the run colored. Every run also writes `summary.html`, the results overview; see [Output files](#output-files).

For more options `cansnper2 type --help`

## Quick start custom databases
Create database, download references annotated in the database and type genomes
```sh
cansnper2 db create -d francisella_tularensis.db --cansnps canSNPs.txt --tree tree.txt --references references.txt
cansnper2 ref get -d francisella_tularensis.db
cansnper2 type sample1.fasta sample2.fasta -d francisella_tularensis.db
```

The source files are checked before a database is built, with the rules of
[CanSNPer2-data](https://github.com/FOI-Bioinformatics/CanSNPer2-data)'s `validate.py`. Files
that do not form one scheme are refused, every problem is listed with its file and line, and no
database is written. The same checks apply to `db import` and `scheme get`.
`db create` builds a new database: an existing database file is refused (a database from
before 2.3.0 with the hint to convert it), and a build that fails part way leaves no file behind.

| File | Refused |
|---|---|
| tree | a line that does not start at the root of the first line, a node with more than one parent, an empty column |
| canSNP | a sequence_position that is not a positive integer, the old header with a `position` column, a base other than A, C, G or T, equal derived and ancestral bases, a genome missing from the references file, a node missing from the tree, a (sequence, sequence_position) on one genome listed twice, a row with too few columns, a rank that is not a positive integer, ranks of a node that are not 1, 2, 3 ... in order, a score that is not a number, a `typed` that is not yes or no, a node with no typed canSNP or with more than one |
| references | a genome listed twice, a row with too few columns |

`#` comment lines, blank lines, a trailing comment column in the canSNP file, `;` separated
trees and nodes without a SNP are accepted.

### Database verbs

```
cansnper2 db create  -d new.db --tree tree.txt --cansnps snps.txt --references refs.txt
cansnper2 db import  -d new.db --tree tree.txt --cansnps snps.txt --references refs.txt
cansnper2 db export  -d db [--format source | newick | tab] [-o outdir]
cansnper2 db info    -d db
cansnper2 db convert  -d db            (or: --cansnps old_snps.txt)
                    --ref-dir references/ -o outdir
cansnper2 db modify  -d db --mod-file changes.txt [--parent NODE] [--replace]
                    [--cansnps snps.txt] [--references refs.txt]
cansnper2 db remove  -d db --parent NODE
```

`db create` and `db import` also record the scheme when you give `--scheme-name`,
`--scheme-status` (verified, draft or unknown) and `--scheme-source`. Every db verb takes
`--log-dir DIR`, `-v/--verbose`, `--debug` and `--quiet`; `cansnper2 db <action> --help`
lists all options of that action.

#### Export

`db export` writes the three source files that `db create` reads, in the format of the
[CanSNPer2-data](https://github.com/FOI-Bioinformatics/CanSNPer2-data) source folder:

| File | Content |
|---|---|
| `<db>_tree.txt` | one line per node, the tab separated path from the root to the node |
| `<db>_canSNPs.txt` | `node rank reference genome sequence sequence_position derived_base ancestral_base score typed` |
| `<db>_references.txt` | `genome strain genbank_id refseq_id assembly_name` |

The canSNP file has the header CanSNPer2-data requires. It holds every canSNP of every node,
ranked within the node, with `typed` marking the ones the scheme types; `score` is empty for a
canSNP taken from a publication. The date column is not exported, and a
database built from the files gets the build date. A database created from these files has the
same tree, SNPs and reference genomes, so a scheme can be exported, edited as text and rebuilt:

```sh
cansnper2 db export -d my_scheme.db -o scheme/
cansnper2 db import -d rebuilt.db --tree scheme/my_scheme_tree.txt \
    --cansnps scheme/my_scheme_canSNPs.txt --references scheme/my_scheme_references.txt
```

A database whose tree has a node with more than one parent, more than one root, or SNPs that
point to a missing reference genome cannot be written in this format and is refused.
`--format newick` writes the tree only. Up to version 2.2.2 the default export (`tab`)
wrote node and link tables that could not rebuild the database; `tab` now gives the source files.

#### Import

`db import` is the counterpart of `db export`: it reads the three source files - the base
CanSNPer2 format - and writes the new database. Export, edit
the text, import back:

```sh
cansnper2 db import -d new_scheme.db --tree scheme/my_scheme_tree.txt \
    --cansnps scheme/my_scheme_canSNPs.txt --references scheme/my_scheme_references.txt
```

The source files carry no scheme information, so `db export` does not write the `scheme_info`
table. `db import` and `db create` record the format, the CanSNPer2 version and the build date,
and the scheme, its status and its source when they are given with `--scheme-name`,
`--scheme-status` (`verified`, `draft` or `unknown`) and `--scheme-source`.

The same loaders `db create` uses do the work, so the schema is always current, and the result
is checked with the export validator before the database is kept: source files that do not
form a clean scheme are refused and leave no half-written database behind. The destination
must not already exist. (Old databases made by releases before the root fixes cannot be fed
in directly; rebuild them from their source files, the ones shipped in
[CanSNPer2-data](https://github.com/FOI-Bioinformatics/CanSNPer2-data), and they come out
single-rooted and exportable.)

### Scheme and reference verbs

```
cansnper2 scheme list   [--repo owner/name] [--branch ref] [--data-dir DIR]
cansnper2 scheme get NAME [-d db] [--branch ref] [--repo owner/name] [--data-dir DIR]
                        [--overwrite] [--no-references] [--allow-unverified]
cansnper2 ref get -d db [--source genbank | refseq] [-o references/]
```

`scheme get` fetches the scheme files, builds the database (default `<scheme>.db`) and then
downloads the references, unless `--no-references`. `ref get` only fetches the references of
an existing database. Both take `--log-dir DIR`, `-v/--verbose`, `--debug` and `--quiet`;
the `scheme` verbs take the same source options (`--repo`, `--branch`, `--data-dir`) as
`scheme list`.

Example structure of references.txt
```
genome strain  genbank_id      refseq_id       assembly_name
OSU18   OSU18   GCA_000014605.1 GCF_000014605.1 ASM1460v1
FSC200  FSC200  GCA_000168775.2 GCF_000168775.2 ASM16877v2
FTNF002-00      FTNF002-00      GCA_000017785.1 GCF_000017785.1 ASM1778v1
LVS     LVS     GCA_000009245.1 GCF_000009245.1 ASM924v1
SCHUS4.2        SCHUS4  GCA_000008985.1 GCF_000008985.1 ASM898v1

```
Example structure of canSNPs.txt (`sequence` is the GenBank accession of the sequence the
canSNP lies on and `sequence_position` the 1-based position on that sequence; `rank` counts
from 1 within a node, `score` is what the canSNP scored when the scheme was made and is empty
for a canSNP taken from a publication, and `typed` marks the canSNPs CanSNPer2 calls, one per
node)
```
node	rank	reference	genome	sequence	sequence_position	derived_base	ancestral_base	score	typed
T/N.1	1	Svensson2009	SCHUS4.2	AJ749949.2	83976	A	G		yes
T.1	1	Svensson2009	SCHUS4.2	AJ749949.2	1165688	G	A		yes
M.1	1	Vogler2009	SCHUS4.2	AJ749949.2	75124	T	C		yes
M.1	2	Vogler2009	SCHUS4.2	AJ749949.2	75299	A	G	480	no
...
```

A node carries its canSNPs at the ranks 1, 2, 3 ... in the order of the file, and exactly one
of them is `typed yes`. CanSNPer2 calls that one; the others are the alternatives of the node,
which an assay may use, and which `db export` writes back. The SNP file of 2.3.0 and earlier,
one canSNP per node with a `snp_id` and a `strain` column, is still read by `--cansnps`: its
canSNPs become rank 1 of their node, typed and without a score.
Example structure of tree.txt
```
T/N.1
T/N.1	T.1		
T/N.1	T.1	B.1
T/N.1	T.1	A/M.1
T/N.1	T.1	A/M.1	A.1
T/N.1	T.1	A/M.1	M.1
...


```
cansnper2 type --help
```
usage: cansnper2 type [-h] [-d DB] [-o DIR] [--save-tree] [--no-snp-files]
                      [--summary] [--node-labels FILE] [--ref-dir DIR]
                      [--allow-missing-references] [--workdir DIR]
                      [--min-required-hits MIN_REQUIRED_HITS]
                      [--strictness STRICTNESS] [--mask MASK] [--keep-going]
                      [--rerun]
                      [--aligner {minimap2,mummer4,progressiveMauve}]
                      [--alignment-timeout MINUTES] [--aligner-path PATH]
                      [--skip-alignment] [--keep-temp] [--tmp-dir DIR]
                      [--log-dir DIR] [-v] [--debug] [--quiet]
                      [query ...]

type query assemblies against a CanSNP database

options:
  -h, --help            show this help message and exit

Required arguments:
  query                 file(s) to align (fasta), or one list file with a file
                        path per line (default: [])
  -d, --database DB     CanSNP database (default: None)

Output options:
  -o, --outdir DIR      output directory; the results overview summary.html is
                        written there, the per-sample files in its samples
                        folder (default: results)
  --save-tree           also write a standalone tree page per sample
                        (samples/<sample>_tree.html, the overview already
                        shows every tree) (default: False)
  --no-snp-files        do not write the per-sample snp files, the call is
                        still printed (default: True)
  --summary             output a summary file and tree with all called SNPs
                        (default: False)
  --node-labels FILE    TSV file annotating tree nodes: <snp
                        name>\t<node|clade>\t<rank>\t<label>[\t<color>]
                        (default: None)

Run options:
  --ref-dir DIR         folder with the reference genomes of the database
                        (default: references/)
  --allow-missing-references
                        type even when reference genomes have no file in
                        --ref-dir (default: False)
  --workdir DIR         working directory the other folder options are
                        relative to (default: ./)
  --min-required-hits MIN_REQUIRED_HITS
                        minimum number of derived SNPs on the path to the
                        final SNP for a call (default: 3)
  --strictness STRICTNESS
                        fraction of the SNPs on the path that must be called
                        (0-1) (default: 0.7)
  --mask MASK           ignore SNPs closer than this many bases to an
                        alignment edge (default: 0)
  --keep-going          continue with the other samples when one fails
                        (default: False)
  --rerun               type samples again although their result file exists
                        (default: False)
  --aligner {minimap2,mummer4,progressiveMauve}
                        program that aligns the query against each reference
                        (default: progressiveMauve)
  --alignment-timeout MINUTES
                        stop the alignments of a query that take longer than
                        this (0 for no limit) (default: 60)
  --aligner-path PATH   the aligner program, or the folder that holds it
                        (default: found in PATH) (default: )
  --skip-alignment      do not align, read the alignments already in --tmp-dir
                        (default: False)
  --keep-temp           keep the temporary alignment files (default: False)
  --tmp-dir DIR         folder for temporary alignments (default: CanSNPer2 in
                        the system temporary folder, see TMPDIR) (default:
                        None)

Logging options:
  --log-dir DIR         folder for the log files (default: logs/)
  -v, --verbose         print process information (default: None)
  --debug               print debug information, and the traceback of an
                        unexpected error (default: None)
  --quiet               do not print warnings (default: 30)
```

### Annotating the tree

`--node-labels FILE` draws meaning that SNP calling itself does not know about onto every
html tree view (the per-sample viewers, the summary tree and the results dashboard):

```
snp	kind	rank	label	color
A/M.1	node	type	Type A	#2471a3
B.1	node	type	Type B	#d35400
B.1	node	subspecies	subsp. anthracis
A.1	node	subspecies	subspecies A group
B.16	node	subspecies	subspecies B group
B.6	clade	marker	AMR feature (separated at the B.6 / B.12 split)	#c0392b
```

- `kind node` draws a colored badge under the SNP dot; a node may carry several badges
  (B.1 above is both a type and a subspecies level).
- `kind clade` draws a dashed rounded box around the whole subtree below the named SNP —
  for a feature that one branch of the tree separates, like an antibiotic resistance
  marker carried by the lineage under B.6.
- `rank` is a free-text level carried by the annotation itself and is never derived from
  the depth in the tree: the database tree has uneven resolution, so B.1 and B.16 can
  both be subspecies levels while standing at different depths. The legend groups
  annotations by rank.
- `color` is optional; without it ranks get a stable palette color. Annotations naming a
  SNP that is not in the database are skipped with a warning.

### Query input

Queries are FASTA files (plain text; gzip-compressed files must be decompressed first), or
FASTQ files of sequencing reads, see [Read typing](#read-typing). A
file is a query list instead when its first non-blank line does not start with `>`: the
list holds one file path per line. This is decided by the content, not by the file name, so
a FASTA file named `x.txt` is typed as a FASTA file. A query list must be the only query
argument; combining it with other arguments stops the run with a message. An empty file, a
file without a `>` header and a file that has a header but no sequence stop the run (or the
sample, with `--keep-going`) with a message naming the file.

### Sample and sequence names

A sample is named after its file name without the extension (`a/x.fasta` gives `x`), and
its results are written as `<outdir>/samples/x_snps.txt` and
`<outdir>/samples/x_not_called.txt`. Two
different files with the same name in one run stop the run before anything is aligned;
rename one or type them into separate `--outdir`. The same file given twice is typed once.

`x_snps.txt` ends with a `# Query file:` line that gives the input it was made from and the
SHA-256 of that file. A later run uses an existing result again (without `--rerun`) only
when the input has the same SHA-256; a different input with the same name is not typed and
the reason is reported. Result files from versions before this line was added are used
again as before.

Sequences are identified by the first word of their header line, in references and in
queries; these names must be set and unique, otherwise the run (reference) or the sample
(query) stops with a message naming the file.

### References of several sequences

Each SNP is stored as the reference sequence it lies on (chromosome or plasmid, named by
its GenBank accession, for example `AE008917.1`) and the 1-based position on that
sequence. A sequence is identified by the first word of its fasta header. The order of
the sequences in the reference file does not matter, and a position is the same as the one
published for that chromosome or plasmid.

Every reference is checked before it is used, at download and at the start of a run: it
must contain each sequence that the database places SNPs on, each sequence must be at least
as long as its last SNP position, and its base at each SNP position must be the ancestral
or the derived base of that SNP. A reference that lacks a named sequence stops with an
error naming the sequence. Extra sequences in the file are allowed. A reference that lacks the
expected base at a fifth of its SNPs or more is refused, fewer mismatches are reported as a
warning.

`cansnper2 ref get --source refseq` renames each sequence of the downloaded RefSeq assembly
to its GenBank accession, so that the file matches the database. A sequence that the NCBI
assembly report marks as different from its GenBank sequence (Relationship `<>`) is not
renamed, as a position on one need not hold on the other; one warning names these
sequences, and the reference check refuses the genome when a SNP lies on one of them.
`--source genbank` (the default) needs no renaming.

The tables `x_snps.txt` and `x_not_called.txt` have the columns `Name`, `Reference`,
`Sequence`, `Pos`, `QuerySeq`, `QueryPosition`, `Ancestral base`, `Derived base`,
`Target base` and `Edge flag`. `Sequence` and `Pos` are the reference sequence and the
position on it. The columns `ReferenceSeq` and `ReferenceSeqPosition` of version 2.2 are
removed, as `Sequence` and `Pos` now carry that information.

`x_snps.txt` has one table row for every called SNP. For a called sample the SNP path is
the called SNPs from the root of the tree down to the final SNP. Derived SNPs that were
called outside the SNP path, on a competing branch, are not in the path; the sample is marked
`review` in `summary.html` with a reason that names them. For a sample that is not called
the path lists all called SNPs.

### Output files

`summary.html` sits directly in `--outdir`, the one file to open; everything a run writes
per sample goes into `--outdir/samples/`, created when the first sample needs it. Each
sample gives these files there:

- `x_snps.txt`: the table of called SNPs (a header row and one row per SNP), followed by
  `#` comment lines:

  ```
  # SNP path: Root;abc_205;1522E;5754E;1344D;3740D
  # Final SNP: 3740D found/depth: 6/6
  # Query file: /data/x.fasta sha256:0930d284...
  # Typing settings: database sha256:9f5a591a... min_required_hits=3 strictness=0.7 mask=0
  ```

  The table can be read with any reader that skips comment lines, for example
  `pandas.read_csv("x_snps.txt", sep="\t", comment="#")`. The header row is written also
  when no SNP was called. Result files of earlier versions, whose trailing lines have no
  `#`, are read as before when a run uses them again.
- `x_not_called.txt`: the same columns for the SNP positions that were not called, and for
  the positions without alignment data (`Target base` `N/A`), followed by the
  `# Query file:` and `# Typing settings:` lines.
- `x_tree.html` (only with `--save-tree`): the tree of the database painted with the states of the
  sample. The header names the sample, the scheme, the final SNP and its support, and the
  run (CanSNPer2 version, aligner, settings, database and its SHA-256, reference genomes,
  date). The SNP path is listed next to the tree, from the root to the final SNP, with the
  ancestral and derived base and the base found in the sample at each SNP; clicking an entry
  moves the tree to that SNP.

Every run with queries writes `summary.html`, the results overview, and `--summary` adds
`snps_summary.txt` and `summary_tree.html`:

- `summary.html` opens in a browser without a server or a network connection. The header
  records the scheme (name, status and source from the scheme information of the
  database), the reference genomes, the database file and its SHA-256, the CanSNPer2
  version, the aligner, `min_required_hits`, `strictness`, `mask` and the date of the run. A
  scheme whose status is draft or unknown, a database without scheme information and a
  scheme with known issues get the notice that the run logs as an amber banner.
  The table lists each sample with its final SNP, the support of the call (derived SNPs
  found on the SNP path / SNPs on the SNP path, for example `6/6`), the number of SNP
  positions on the SNP path without alignment data, and a status:
  - `clean`: the SNP path is confirmed; every SNP on it has alignment data, none is near an
    alignment edge, and no derived SNP lies outside it;
  - `review`: the SNP path is confirmed, but SNPs on it have no alignment data or are
    flagged at an alignment edge, or derived SNPs were called outside it;
  - `not called`: no SNP path met `--min-required-hits` and `--strictness`, no derived SNP
    was found, or two SNP paths competed at the same depth;
  - `failed`: the alignment or the processing of the sample failed;
  - `reused`: the result was read from an earlier run (`--rerun` computes it again).

   Samples that need a look are listed first. The page opens on a support overview of the
   whole run: every sample's SNP path drawn on the tree at once, the edge to a node and the
   node itself growing with the number of paths running through it, `×N` counts at the names
   and the final SNPs with their counts in the collapsible detail - where the genomes of the
   run ended up. Two filter boxes narrow the table, by sample name and by any SNP in a
   sample's SNP path (also a SNP without alignment data, e.g. `B.33`); the Overview button
   returns to the joint view. Clicking a sample paints the tree with its
   states and shows its SNP path, the reasons for its status, links to its result files in
   `samples/` and the labels that `--node-labels` gives the final SNP; the side sections
   start collapsed and open by hand. In the tree, derived SNPs are green,
  ancestral SNPs grey, another base blue, a position without alignment data an open circle,
  and a SNP that was not tested (a reused result) a small open circle; the SNP path is drawn
  as a thick green line and an amber ring marks a SNP near an alignment edge. "Download
  table (TSV)" saves the table with the run as `#` lines, "Save SVG" saves the tree as it is
  drawn, and the page prints with the header, the table and the tree.
- `snps_summary.txt` starts with `#` lines on the run (version, date, database and its
  SHA-256, scheme, reference genomes, aligner, settings), followed by the header row
  `sample`, `final_snp`, `support`, `snp_path` and one tab separated row per sample. `support`
  is `NA` for a sample that was not called.
- `summary_tree.html`: the tree with the final SNPs of the run in green.

### Upgrading from 2.3

Version 2.4.0 changes where the canSNPs of a scheme come from. Up to 2.3.0 a scheme had a
SNP file, `<name>_snp.txt`, with one canSNP per node. From 2.4.0 it has a canSNP file,
`<name>_canSNPs.txt`, which holds every canSNP of every node, ranked within the node, with
`typed` marking the ones CanSNPer2 calls. The database stores all of them, so `db export`
writes them back, and a later release can call more than one canSNP per node without
another change of format.

A database of 2.3.0 is refused:

```
The database old.db cannot be used. Its snp_annotation table holds one canSNP per node,
without the rank, score and typed columns of the current format. Convert it with cansnper2 db convert,
or fetch the current scheme with cansnper2 scheme get.
```

`db convert` upgrades it, and a `<name>_snp.txt` as well. Neither needs reference files,
because no position moves:

```sh
cansnper2 db convert -d my_scheme.db -o upgraded/
cansnper2 db convert --cansnps my_scheme_snp.txt -o upgraded/
```

`--cansnps` also reads a `<name>_snp.txt` as it is, as rank 1 of each node,
typed and without a score, so a scheme in that format still builds.

### Upgrading from 2.2

Version 2.3.0 changed how a SNP position is stored. Up to 2.2 it was one number along all
the sequences of the reference joined in the order of the reference fasta file, which made
a result depend on that order. From 2.3.0 the database and the source files hold the
sequence and the position on it (`sequence`, `sequence_position`; the `position` column
is gone), and the result tables hold the columns listed above.

A database or SNP source file in the old format is refused:

```
The database old.db cannot be used. It stores SNP positions along the joined reference
sequences (the format before 2.3.0). Convert it with cansnper2 db convert, or fetch
the current scheme with cansnper2 scheme get.
```

A database made by 2.3.0 cannot be read by 2.2.x, which stops with `no such column: position`.

`db convert` brings a SNP source file or a database of an older format to the current one, in
one step, and writes it to the folder given with `-o`. The input is not changed.

A file or database from before 2.3.0 has its positions located on the sequence of its reference
genome, which is read from `--ref-dir` (`<genome>.fna`; missing files are downloaded). A SNP whose
base in the reference is neither its derived nor its ancestral base is reported, and nothing is
written. A file or database of 2.3.0 only gains the rank, score and typed columns of the canSNP
format, so it needs no reference files: every canSNP becomes rank 1 of its node, typed, without
a score. A converted SNP file is written as the canSNP file of its scheme
(`<name>_snp.txt` becomes `<name>_canSNPs.txt`).

```sh
# a SNP source file from before 2.3.0; --references is optional and names the genomes of the file
cansnper2 db convert --cansnps old_snp.txt --references old_references.txt \
    --ref-dir references/ -o converted/

# a database from before 2.3.0
cansnper2 db convert -d old.db --ref-dir references/ -o converted/

# a SNP file or database of 2.3.0: no reference files are needed
cansnper2 db convert --cansnps my_scheme_snp.txt -o converted/
cansnper2 db convert -d my_scheme.db -o converted/
```

The reference files must have their sequences in the order in which the old positions were
counted (the order of the NCBI assembly file). The sequence names written are the first
words of their fasta headers; missing reference files are downloaded from GenBank, so the
names are GenBank accessions. progressiveMauve alignments made by 2.2 must be made again,
as `--skip-alignment` needs the sequence table (`<alignment>.sequences`) that 2.3.0 writes
next to an xmfa alignment. minimap2 and MUMmer4 alignments made by 2.2 are read as they
are, because their reference coordinates were always positions on each sequence.

### When something is wrong

Every command checks its inputs before it starts to write or to align, and a problem ends
the command with one message that names the file or the database and says what to do. The
exit status tells the kind of stop:

- 0: the command did what it was asked to do (warnings may have been printed);
- 1: an input could not be used; the message says which and why;
- 2: an unexpected error, printed as `CanSNPer2 stopped: <type>: <message>; run with --debug
  for details`. Run the command again with `--debug` to see the traceback, and report it
  together with the command line.

A database is never left half changed: `db create`, `db import` and `scheme get` remove
a database they had started, and a modification (the files of `db modify` and `db remove`)
is made on a copy that replaces the database
only when every step succeeded.

| Message (beginning) | Command | What to do |
|---|---|---|
| `The source files cannot be built into a database, N problem(s):` followed by `<file>:<line>: ...` | `db create`, `db import`, `scheme get` | Correct the listed rows; all problems are listed together. |
| `<file>:<line>: genbank_id '...' is not an assembly accession` | as above, and `--references` | Give the GenBank (GCA_) or RefSeq (GCF_) assembly accession with its version, such as GCA_000008985.1. |
| `<file>:<line>: genbank_id '...' of genome '...' is already used by genome '...'` | as above, and `--references` on an existing database | Each genome needs its own assembly accession, as one reference file is written per accession. |
| `The files cannot be added to <db>, N problem(s):` | `db modify` with `--cansnps` or `--references` | A canSNP on a node that is not in the tree of the database or the modification file, on a genome that is not listed, a short row, or a genome, a rank of its node, a second typed canSNP of its node, or a position the database already holds. Correct the files; the database was not changed. |
| `The selected parent node (X) could not be found`, `The node X could not be found`, `X is a root of the tree` | `db modify`, `db remove` | Give a node of the database with `--parent`; the root cannot be removed. |
| `The SNP X could not be written to <db>`, `The node X could not be written to <db>` | `db modify`, `db remove` | The disk is full or the database cannot be written to; free space or check the permissions, then run again. |
| `scheme X on owner/name@ref uses positions on the joined reference (before 2.3.0)` | `scheme get` | The branch holds the format before 2.3.0. Use CanSNPer2 2.2 for it, or a branch converted for 2.3.0. |
| `... manifest.json was not found ...` or `... records no hash for ...`, `the files of the verified scheme X cannot be checked` | `scheme get` | The files cannot be checked against the manifest of the repository. Report it to the repository, or give `--allow-unverified` when the source is trusted. |
| `The database X could not be built from the scheme Y`, `The database X could not be written` | `scheme get` | An sqlite error or a file system error; free space or check the permissions of the folder. No partial file is left behind. |
| `X lists no reference genomes, there is nothing to download` | `ref get` | Build the database again with its references file. |
| `REF: the database has no genbank_id for this genome` | `ref get` | Give the assembly accession of the genome in the references file and build the database again. |
| `The SNP positions of X could not be read` | `ref get`, `type` | The database is damaged; export it with `db export`, correct the files and build it again. |
| `N reference genome(s) of X have no file <genome>.fna in DIR: ...` | `cansnper2 type` | Fetch them with `cansnper2 ref get -d X -o DIR`. Without them the SNPs on these genomes would not be compared. To type without them on purpose, for example to test a scheme with a reference genome left out, give `--allow-missing-references`: the run then continues with one warning that names the missing genomes. It still stops when no genome of the database has a file. |
| `The database X has no canSNPs`, `The database X lists no reference genomes` | `cansnper2 type` | Build the database again with its canSNP file and its references file. |
| `The database X has N canSNP(s) but none of them is marked typed` | `cansnper2 type` | Every node needs one canSNP marked `typed yes`; build the database again from a canSNP file that marks them. |
| `The database X has N SNP(s) on node ids that are not in the nodes or tree table` | `cansnper2 type` | Export the database, correct the tree or canSNP file and build it again. Nodes without canSNPs are allowed. |

These cases give a warning and keep exit status 0: a draft scheme or a scheme of unknown
status; a reference whose base matches neither the ancestral nor the derived base at fewer
than a fifth of its SNPs; reference files in `--ref-dir` that are not reference genomes with
SNPs in the database; reference genomes without a file in `--ref-dir` when
`--allow-missing-references` is given; a database changed after it was built; `--no-references`; and a draft scheme, or a scheme read with
`--data-dir`, whose files cannot be checked against a manifest.

## Aligners

CanSNPer2 aligns each query genome against every reference genome of the database
and reads the query base at each SNP position from the alignment. Three aligners are
available:

* `progressiveMauve` (default), which writes xmfa alignments. The xmfa counts reference
  positions along the reference sequences joined end to end, so the names and lengths of
  those sequences are written next to each alignment (`<alignment>.sequences`), and
  `--skip-alignment` reads them again. This table is written for progressiveMauve only.
* `minimap2`, run as `minimap2 -c --cs=long -x asm20 -r 1k,10k -N 50 -p 0.1` and
  written as PAF. Only primary alignments with a mapping quality of at least 20 are
  read: a secondary alignment, or a primary one with a lower mapping quality, can carry
  the base of another copy of a repeated region. The long-join bandwidth (`-r`) and the
  secondary chain settings (`-N`, `-p`) differ from the preset so that the short
  stretches between IS elements and rearrangements, which separate Francisella genomes
  of different subspecies, are aligned as well. PAF gives positions on each reference
  sequence, so the alignment is read again with `--skip-alignment` without any other file.
* `mummer4`, MUMmer4 `nucmer` with its default settings and one thread (`-t 1`) followed
  by `delta-filter -1`, written as delta. With the default of two threads, repeated runs on
  the same genomes could give different alignments after the filter, and so a different
  base or no data at a few canSNP positions; one thread gives the same result every time. The 1-to-1 filter keeps, for every reference and query position,
  only the best alignment, allowing rearrangements. The delta file names the reference
  and query fasta files it was made from. The reference coordinates are positions on each
  reference sequence; the reader takes only the query bases from the query file, so that
  file must still be in place when the alignment is read again with `--skip-alignment`,
  and the reference file is not needed. `nucmer` and `delta-filter` must
  both be installed (bioconda package `mummer4`).

On 937 F. tularensis genomes typed with the Francisella scheme, both minimap2 and
MUMmer4 gave the same final call as progressiveMauve for every genome, and neither
had a SNP position with a different base:

| | minimap2 | MUMmer4 |
|---|---|---|
| positions with the same base as progressiveMauve | 305,513 | 305,693 |
| positions with data only from progressiveMauve | 291 | 111 |
| positions with data only from this aligner | 169 | 177 |
| time for the 937 genomes, three parallel runs | 12.9 min | 15.0 min |

On 268 genomes of other Francisellaceae species and genera, neither gave a final call.
progressiveMauve, run on a subset of 77 of them, called 5 genomes of other genera
F. tularensis. minimap2 does not align genomes of other genera; MUMmer4 does, at about
43 s per genome against under 4 s for minimap2.

Each aligner is a class in `CanSNPer2/modules/Aligners.py` that gives the command line,
the alignment file ending, how its exit codes are judged and the reader of its output.
A reader is a subclass of `AlignmentReader` (`CanSNPer2/modules/AlignmentReader.py`)
that walks the alignment file and reports, for each SNP position it covers, the query
base on the forward strand of the reference (`observe()`); the classification into
derived, ancestral, other base and no data, `--mask`, the edge flag and the result
rows are shared by all readers. SNP positions are the reference coordinates of the
database: the reference sequence (the first word of its fasta header) and the position on
it. progressiveMauve counts along the sequences joined end to end and its reader converts
the positions. A new aligner is registered in `ALIGNERS` and can then be selected with `--aligner`.

Additional information about running CanSNPer2 could be found in the [wiki](https://github.com/FOI-Bioinformatics/CanSNPer2/wiki).

About this software
===================
Copyright (C) 2019 David Sundell @ FOI bioinformatics group  

CanSNPer2 is implemented as a Python 3 package. It is open source software made available
under the [GPL-3.0 license](LICENSE).

If you experience any difficulties with this software, or you have suggestions, or want
to contribute directly, you have the following options:

- submit a bug report or feature request to the
  [issue tracker](https://github.com/FOI-Bioinformatics/CanSNPer2/issues)
- contribute directly to the source code through the
  [github](https://github.com/FOI-Bioinformatics/CanSNPer2) repository. 'Pull requests' are
  especially welcome.

## Read typing

CanSNPer2 types sequencing reads without an assembly. Read files are recognised by their
content (FASTQ, plain or gzip-compressed); a run types either reads or assemblies, not both.

```
CanSNPer2 --database brucella_genus.db --refdir references/ sample_R1_001.fastq.gz sample_R2_001.fastq.gz
CanSNPer2 --database brucella_genus.db --refdir references/ --platform ont barcode05.fastq.gz
CanSNPer2 --database brucella_genus.db --refdir references/ --samplesheet samples.tsv
```

**Samples.** Paired-end files are found by their names: `_R1`/`_R2` (with or without the
Illumina `_001`) or `_1`/`_2`. The sample is named after the file without the mate part and
the endings (`B_iv_01_M_S1_R1_001.fastq.gz` gives `B_iv_01_M_S1`). A file named as one mate
whose other mate was not given stops the run. Reads of `--platform ont` or `hifi` are never
paired. `--samplesheet` gives the samples explicitly instead of query arguments: a tab
separated file with the header `sample`, `reads1` and optionally `reads2` and `platform`;
relative paths are read from the folder of the sheet.

**Mapping and counting.** For each reference of the database, the reads are mapped with
`minimap2 -a -x <preset>`, sorted and indexed with samtools, and `samtools mpileup -a -B -A
-d 0 -q <min_mapq> -Q <min_baseq>` lists the read bases at the canSNP positions only. Reads
below the mapping quality and bases below the base quality are not counted; duplicate,
secondary and QC-failed reads are left out (mpileup defaults), and a base covered by both
reads of a pair is counted once.

**The call of a canSNP** from its counts of A, C, G and T:

| counts | state |
|---|---|
| no reads | no data |
| fewer than `--min_depth` | no data (low depth) |
| derived base at least `--min_allele_fraction` | derived |
| ancestral base at least `--min_allele_fraction` | ancestral |
| another base at least `--min_allele_fraction` | other base |
| no base reaches `--min_allele_fraction` | no data (mixed) |

The tree walk and `--min_required_hits` and `--strictness` are those of assembly typing.

**Platform defaults**, which the options override:

| `--platform` | minimap2 preset | `--min_mapq` | `--min_baseq` | `--min_depth` | `--min_allele_fraction` |
|---|---|---|---|---|---|
| illumina (default) | sr | 20 | 20 | 5 | 0.8 |
| ont | map-ont | 20 | 10 | 5 | 0.7 |
| hifi | map-hifi | 20 | 20 | 5 | 0.8 |

On reads simulated from 30 *Brucella* and 20 *F. tularensis* genomes (wgsim for Illumina,
badread for ONT, 15 kb reads with 0.1 % errors for HiFi) and typed against their assembly
calls, no read call was on another branch than the assembly call, at any depth or platform.
With fewer reads than needed the call stops higher up the tree or is NA:

| | 5x | 10x | 20x | 30x | 50x |
|---|---|---|---|---|---|
| Illumina, same call as the assembly | 0 of 50 | 37 of 50 | 40 of 50 | | 41 of 50 |
| ONT | | 42 of 50 | | 50 of 50 | |
| HiFi | | 46 of 50 | | 50 of 50 | |

The Illumina misses from 20x on are *Brucella* species far from the reference *B. melitensis*
16M, whose reads map at too few canSNPs (see [Local references](#local-references)). The ONT
reads of badread have errors close to random; systematic errors of real ONT reads, as at
homopolymers, were not part of the test.

**Mixtures.** CanSNPer2 types isolates. A sample with at least `--mixed_sites` (default 3)
mixed canSNP sites is reported as possibly mixed, in the log, in the results overview and
in the `Reads:` line of its result file; the strains in it are not resolved. The flag only
asks for a review; it never changes the call. A canSNP in a region the sample has twice while
the reference has it once is mixed too, as reads from the second copy map to the single copy
of the reference: in the evaluation above, 1 of the 50 genomes (an *F. tularensis* draft with
90-96 % copies of three canSNP regions on short contigs) was flagged with 3 to 5 mixed sites,
the other 49 had at most 2. Mixtures of two genomes on different branches had 8 to 17 mixed
sites at 50/50, 3 to 11 at 80/20 and 0 to 3 at 95/5; no mixture was called on another branch
than its majority genome, the call stopped higher up or was that of the majority genome.

**Output.** The result files have the columns of assembly typing, with `QuerySeq` `reads`
and `QueryPosition` the read depth. `<sample>_snps.txt` has the extra comment line
`Reads: platform=... files=... sites_with_reads=... median_depth=... mixed_sites=...` (the median
depth is over the canSNP sites that have reads), and the typing settings
line names the read thresholds. `<sample>_site_counts.tsv` lists every canSNP with its counts
of A, C, G and T, the depth, the fraction of the derived base and its state. `--keep_temp`
keeps the sorted reads (`<reference>_<sample>.pileup.bam`) and the pileups, and
`--skip_alignment --tmpdir` types from them again. `--mask` does not apply to reads.

### Local references

Reads map to a reference only where they are close to it: short reads need roughly 95 %
identity over their length, while an assembly aligns with far less. For a scheme that spans
a genus and has one reference genome, many canSNPs of the species far from that reference
get no reads, and the call stops above them. Two *B. anthropi* read sets typed with
`brucella_genus` (reference *B. melitensis* 16M) had reads at 169 and 171 of its 283 canSNPs,
and stopped one node above `BR.anthropi`; one of the draft assemblies of the same samples was
called at `BR.anthropi.1.1.1.1`.

`--local_references` adds rounds against a genome of the scheme closer to the sample:

1. A round is made only when it can change the call: when there is no call, or when a canSNP
   on the path to the call or of a child of the call has no data (no reads, too few, or not
   lifted). A terminal call with data on its whole path, or a call whose children all have
   data, is settled, and no round is made. This is checked before every round.
2. The local reference is one of the representative genomes of the terminal nodes, from the
   genomes file of the scheme (`CanSNPer2-database --genomes`). It is found by walking down
   the tree from the call of the first round (from the root without a call). At each node,
   the child subtrees that hold representatives are the options, less those whose own canSNP
   is ancestral. When the derived canSNPs of the first round lie in one of them, the walk
   goes there.
3. Otherwise a subsample of the reads (`--screen_reads`, default 100,000 per file) is mapped
   at once to up to three representatives of each option, from different terminal nodes, and
   the walk follows the subtree of the genome with the most matched bases. When that genome
   matches less than 10 % of the screened bases, the reads are not from the genomes below,
   and the walk stops there. The walk ends at a terminal node, whose representative is the
   local reference. When two subtrees get similar shares of the screened bases (the second at
   least 20 % and at least half the first), a warning says so, unless the final counts
   separate the two: the canSNP of one subtree derived and that of the other ancestral.
   Sibling genomes share most of their sequence, so a close screen alone is expected between
   them.
4. The canSNP positions are lifted from the scheme reference to the local reference by
   aligning the two genomes with minimap2 (asm20, as assembly typing does); a canSNP is lifted
   only when the local genome has its ancestral or derived base there. The lift is cached in
   `<refdir>/local/lift/`. A local genome that carries the ancestral base at more of the
   lifted canSNPs on the path to its own node than the derived base is not used.
5. The reads are mapped to the local reference and counted at the lifted positions. Per
   canSNP, these counts are used when they reach `--min_depth`, otherwise those of the
   first round. When both reach it and give different bases, the canSNP is a `conflict` and
   counts as no data.
6. This repeats from the new call while the chosen genome changes, at most `--local_rounds`
   times (default 2).

The local genomes are read from `<refdir>/local/<accession>.fna`;
`CanSNPer2-download --database <db> --outdir <refdir> --local_references` fetches the
representatives of a scheme (221 genomes, about 860 MB, for `brucella_genus`). Alternatively,
`--fetch_local` downloads from NCBI during typing only the genomes that step 2 or 3 uses and that
are not in the folder yet: the chosen genome when the counts decide, the screened genomes
otherwise. A downloaded genome stays in the folder for later runs; a genome that cannot be
fetched is left out, as a missing file is. Typing then depends on network access to NCBI, so the
option is off by default; on compute nodes without it, fetch the genomes beforehand.

The site count file gives for each canSNP the genome its counts
came from (`source`) and the lifted position (`lifted_position`, contig:position and strand),
and the `Reads:` line names the local reference, its node, its share of the screened bases at
the last screen, the number of screens and the number of conflicts; with `--fetch_local` it also
gives the number of genomes fetched for the sample (`fetched=`).

**Evaluation.** The nine *Brucella* genomes of the evaluation above whose 50x Illumina reads
missed their assembly call against 16M (all outside the classical cluster) were typed again
with `--local_references`; when a genome was itself a representative (five of them), it was
taken out of the local references, so that the sample was not mapped to its own genome. The
other four are members of a representative at 99 % average nucleotide identity or closer.

| genome | assembly call | 16M only | with local references |
|---|---|---|---|
| GCA_039674445.1 | `BR.chromiisoli` | NA | `BR.chromiisoli` |
| GCF_026256505.1 | `BR.intermedia.1.4` | NA | `BR.intermedia.1.4` |
| GCA_031999685.1 | `BR.intermedia.1.1.1` | `BR.1.1.1.1.1.1` | `BR.intermedia.1.1.1` |
| GCF_000182645.1 (held out) | `BR.intermedia.1.3.1` | NA | `BR.intermedia.1.3.1` |
| GCF_008932555.1 | `BR.intermedia` | NA | `BR.intermedia` |
| GCA_008932375.1 (held out) | `BR.haematophila` | `BR.1.1.1.1.1.2` | `BR.haematophila` |
| GCF_003664555.1 (held out) | `BR.teleogrylli.1.1` | `BR.1.1.1.1.1.2` | `BR.teleogrylli` |
| GCA_018798825.1 (held out) | `BR.sp018798825` | NA | `BR.sp018798825` |
| GCF_003049685.1 (held out) | `BR.pituitosa.1.1.2` | NA | `BR.pituitosa.1.1.2` |

Eight reached the assembly call and one stopped a node above it; none was called on another
branch. The canSNPs with reads rose from 60-158 to 175-259 of 283. The counts of the first
round chose the local reference in eight samples without a screen; one needed three screens.
One sample (*B. pituitosa*) went first to a genome of another species and reached its own in
a second round. The two refbio *B. anthropi* read sets (B_iv_02, B_iv_04) went from one node
above `BR.anthropi` to `BR.anthropi.1.1.1.1`, the call of the draft assembly of B_iv_02, with
one local round each.

