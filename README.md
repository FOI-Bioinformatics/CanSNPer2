# CanSNPer2
CanSNPer2: A toolkit for SNP-typing NGS data.

Future planned implementations
* Read input - Call SNPs directly from reads instead of an assembly

Databases supplied can be found at https://github.com/FOI-Bioinformatics/CanSNPer2-data

* Francisella tularensis
* Bacillus anthracis
* Yersinia pestis
 
#### Upcoming databases
* Brucella
* Coxiella


## Installation
Installing using bioconda using mamba (recommended instead of conda)

```
mamba create -n cansnper2 cansnper2
mamba activate cansnper2
```

Installing from the repository
```
python setup.py install
```

## Requirements (for manual install conda will install all dependencies)

* progressiveMauve (the default aligner, see [Aligners](#aligners))
* minimap2 (optional, `--aligner minimap2`)

Read input (not yet implemeted)
* BWA
* Minimap2

## User guide CanSNPer2 (for custom databases see below)
1. Build a database from a scheme in [CanSNPer2-data](https://github.com/FOI-Bioinformatics/CanSNPer2-data) and download its reference genomes

```sh
CanSNPer2-download --list_schemes                       # schemes on master (verified)
CanSNPer2-download --scheme francisella_tularensis      # builds francisella_tularensis.db, downloads references/
```

The scheme files (SNPs, tree and references) are fetched and checked against the SHA-256
recorded in the manifest of CanSNPer2-data, and the database is built locally by the installed
CanSNPer2. `francisella_tularensis.source.json` beside the database records the repository
commit and the hash of each file it was built from. A verified scheme whose files cannot be
checked (no manifest, or a file the manifest has no hash for) is refused unless
`--allow_unverified` is given. `--branch dev` gives access to draft
schemes that have not yet been verified; building one prints a warning. `--branch` also takes a
commit, which makes a database reproducible from its `source.json`.

Custom schemes can be read the same way from a local directory or from another GitHub
repository, as long as it is laid out as CanSNPer2-data. Only
`schemes/<name>/<name>_{snp,tree,references}.txt` is required; `schemes.toml`, the scheme
`README.md` (with its Status) and `database/<name>.manifest.json` are used when present.

```sh
CanSNPer2-download --data_dir my-schemes --list_schemes
CanSNPer2-download --data_dir my-schemes --scheme my_scheme       # a local folder or clone
CanSNPer2-download --repo owner/my-schemes --scheme my_scheme     # default branch of that repository
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
as known issues, and the CanSNPer2 version and date of the build. `CanSNPer2-database --database
<db> --info` prints the table. Typing, `CanSNPer2-download --database` and `--export` report
the scheme once when they open the database; a draft or unknown status, known issues, or a
database without the table (built before 2.3.0 or by hand) give a warning, and a verified scheme
is reported only with `--verbose`. The table is informative only and never stops a run.
A database changed with `--mod_file`, `--annotation`, `--references` or `--remove` no longer
matches the scheme files it was built from, so the change sets its status to `unknown`, adds
`+modified` to the commit (or records `modified` when there was none), records the date in
`modified` and the CanSNPer2 version in `built_with`; typing then gives a warning. A database
converted with `--convert` records status `unknown` and the converted file as its source.
`--info` prints `none` for an empty value, such as `known_issues` of a scheme without any.

```
$ CanSNPer2-database --database brucella_suis.db --info
scheme	brucella_suis
status	draft
source	FOI-Bioinformatics/CanSNPer2-data
ref	dev
commit	abc1234000000000000000000000000000000000
known_issues	B.1
format	2.3
built_with	2.3.0
built	2026-11-21
```

2. Alternatively, with a database you already have, download its references

```sh
CanSNPer2-download --database my_database.db
```

3. Run genomes
```sh
CanSNPer2 --database francisella_tularensis.db fastadir/*.fasta --summary
```

The summary parameter gives `snps_summary.txt`, a table of the final SNP, its support and the SNP path of every sample, and a tree (self-contained html, opens in any browser) with the final SNPs of the run colored. Every run also writes `summary.html`, the results overview; see [Output files](#output-files).

For more options CanSNPer2 --help

## Quick start custom databases
Create database, download references annotated in the database and run CanSNPer2
```sh
CanSNPer2-database --database francisella_tularensis.db --annotation canSNPs.txt --tree tree.txt --reference references.txt --source_type CanSNPer --create
CanSNPer2-download --database francisella_tularensis.db
CanSNPer2 sample1.fasta sample2.fasta --database francisella_tularensis.db --save_tree
```

The source files are checked before a database is built, with the rules of
[CanSNPer2-data](https://github.com/FOI-Bioinformatics/CanSNPer2-data)'s `validate.py`. Files
that do not form one scheme are refused, every problem is listed with its file and line, and no
database is written. The same checks apply to `--import` and `CanSNPer2-download --scheme`.
`--create` builds a new database: an existing `--database` file is refused (a database from
before 2.3.0 with the hint to convert it), and a build that fails part way leaves no file behind.

| File | Refused |
|---|---|
| tree | a line that does not start at the root of the first line, a node with more than one parent, an empty column |
| canSNP | a sequence_position that is not a positive integer, the old header with a `position` column, a base other than A, C, G or T, equal derived and ancestral bases, a genome missing from the references file, a node missing from the tree, a (sequence, sequence_position) on one genome listed twice, a row with too few columns, a rank that is not a positive integer, ranks of a node that are not 1, 2, 3 ... in order, a score that is not a number, a `typed` that is not yes or no, a node with no typed canSNP or with more than one |
| references | a genome listed twice, a row with too few columns |

`#` comment lines, blank lines, a trailing comment column in the canSNP file, `;` separated
trees and nodes without a SNP are accepted.

### CanSNPer2-database

```
usage: CanSNPer2-database [-h] [-db] [--tree] [--annotation] [--references]
                          [--source_type] [--create] [--import]
                          [--scheme_name] [--scheme_status] [--scheme_source]
                          [--info] [--convert] [--refdir] [--mod_file]
                          [--parent] [--remove] [--replace] [--export]
                          [--export_format] [-o] [--tmpdir] [--logs]
                          [--verbose] [--debug] [--supress]

CanSNPer2-database - create, modify and export CanSNP databases

options:
  -h, --help         show this help message and exit

Required:
  -db , --database   CanSNPer2 database name

Load database:
  --tree             CanSNPer tree source file
  --annotation       CanSNPer snp source file
  --references       File containing all reference genomes listed
  --source_type      Select source file type
  --create           Create new database!
  --import           build a database from the source files written by
                     --export (--tree --annotation --references), written to
                     --database
  --scheme_name      with --create or --import: the name of the scheme,
                     recorded in the database
  --scheme_status    with --create or --import: the status of the scheme,
                     verified, draft or unknown, recorded in the database
  --scheme_source    with --create or --import: where the scheme files came
                     from (a repository or folder), recorded in the database

Database information:
  --info             print the scheme information the database records
                     (scheme, status, source, CanSNPer2 version) and exit

Convert from the format before 2.3.0:
  --convert          convert a SNP source file (--annotation) or a database
                     (--database) from positions along the joined reference
                     sequences to positions on each reference sequence,
                     written to -o
  --refdir           folder with the reference files <genome>.fna used by
                     --convert; missing files are downloaded

Database modifications:
  --mod_file         File with modifications/update to the tree
  --parent           Node (or nodes matching tree file) from which to
                     update/replace/remove
  --remove           Remove the whole branch rooted at the --parent node, with
                     its nodes, links and SNP annotations
  --replace          replace node

Export database:
  --export           Export the database to text files
  --export_format    source: the tree, SNP and reference files that --create
                     reads; newick: the tree; tab is the old name of source
                     (default source)
  -o , --outdir      outdir for database export!

Logging and debug options:
  --tmpdir           Specify tmp directory default (/tmp)
  --logs             Specify log directory
  --verbose          print process info, default no output
  --debug            print debug info
  --supress          supress warnings

available CanSNPer2 commands (installed together):
  CanSNPer2            type query assemblies or reads against a CanSNP database
  CanSNPer2-database   create, modify and export CanSNP databases
  CanSNPer2-download   fetch the reference genomes a database needs from NCBI
  CanSNPer2-test       check the installation: imports and aligner programs in PATH
Run any command with --help for its arguments; examples in README.md.
```

#### Export

`--export` writes the three source files that `--create` reads, in the format of the
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
CanSNPer2-database --database my_scheme.db --export -o scheme/
CanSNPer2-database --database rebuilt.db --import --tree scheme/my_scheme_tree.txt \
    --annotation scheme/my_scheme_canSNPs.txt --references scheme/my_scheme_references.txt
```

A database whose tree has a node with more than one parent, more than one root, or SNPs that
point to a missing reference genome cannot be written in this format and is refused.
`--export_format newick` writes the tree only. Up to version 2.2.2 the default export (`tab`)
wrote node and link tables that could not rebuild the database; `tab` now gives the source files.

#### Import

`--import` is the counterpart of `--export`: it reads the three source files - the base
CanSNPer2 format - and writes a database in the current format to `--database`. Export, edit
the text, import back:

```sh
CanSNPer2-database --import --database new_scheme.db --tree scheme/my_scheme_tree.txt \
    --annotation scheme/my_scheme_canSNPs.txt --references scheme/my_scheme_references.txt
```

The source files carry no scheme information, so `--export` does not write the `scheme_info`
table. `--import` and `--create` record the format, the CanSNPer2 version and the build date,
and the scheme, its status and its source when they are given with `--scheme_name`,
`--scheme_status` (`verified`, `draft` or `unknown`) and `--scheme_source`.

The same loaders `--create` uses do the work, so the schema is always current, and the result
is checked with the export validator before the database is kept: source files that do not
form a clean scheme are refused and leave no half-written database behind. The destination
must not already exist. (Old databases made by releases before the root fixes cannot be fed
in directly; rebuild them from their source files, the ones shipped in
[CanSNPer2-data](https://github.com/FOI-Bioinformatics/CanSNPer2-data), and they come out
single-rooted and exportable.)

### CanSNPer2-download

```
usage: CanSNPer2-download [-h] [-db] [--scheme] [--branch] [--repo]
                          [--data_dir] [--list_schemes] [--overwrite]
                          [--no_references] [--allow_unverified] [-s] [-o]
                          [--logs] [--verbose] [--debug]

CanSNPer2-download - fetch the reference genomes a database needs from NCBI

options:
  -h, --help          show this help message and exit

Database:
  -db , --database    CanSNP database (with --scheme: the database to build,
                      default <scheme>.db)

CanSNPer2-data scheme options:
  --scheme            Fetch this scheme from CanSNPer2-data and build its
                      database
  --branch            Branch, tag or commit of the repository; in
                      CanSNPer2-data master holds verified schemes, dev also
                      draft schemes (default master, with --repo the default
                      branch of that repository)
  --repo              GitHub repository laid out as CanSNPer2-data, as
                      owner/name (default FOI-Bioinformatics/CanSNPer2-data)
  --data_dir          Local directory laid out as CanSNPer2-data (a clone, or
                      a folder of custom schemes) to read schemes from instead
                      of GitHub
  --list_schemes      List the schemes of the repository or --data_dir and
                      exit
  --overwrite         Replace an existing database when building from --scheme
  --no_references     Build the database from --scheme without downloading its
                      references
  --allow_unverified  Build a verified scheme from GitHub although its files
                      cannot be checked against a manifest (no manifest, or no
                      hash for a file)

Download options:
  -s , --source       Source for download (genbank/refseq)
  -o , --outdir       reference genomes folder

Logging and debug options:
  --logs              Specify log directory
  --verbose           Verbose logging
  --debug             Debug logging, and the traceback of an unexpected error

Give a database with --database, or name a scheme with --scheme to fetch it from CanSNPer2-data (or --repo, or --data_dir), build its database and download its references in one step.
available CanSNPer2 commands (installed together):
  CanSNPer2            type query assemblies or reads against a CanSNP database
  CanSNPer2-database   create, modify and export CanSNP databases
  CanSNPer2-download   fetch the reference genomes a database needs from NCBI
  CanSNPer2-test       check the installation: imports and aligner programs in PATH
Run any command with --help for its arguments; examples in README.md.
```
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
which an assay may use, and which `--export` writes back. The SNP file of 2.3.0 and earlier,
one canSNP per node with a `snp_id` and a `strain` column, is still read by `--annotation`: its
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
CanSNPer2 help
```
usage: CanSNPer2 [-h] [-db] [-o DIR] [--save_tree] [--no_snpfiles]
                       [--summary] [--annotations FILE] [--refdir]
                       [--allow_missing_references] [--workdir] [--read_input]
                       [--min_required_hits MIN_REQUIRED_HITS]
                       [--strictness STRICTNESS] [--mask MASK] [--keep_going]
                       [--rerun] [--aligner] [--alignment_timeout MINUTES]
                       [--aligner_path] [--skip_alignment] [--keep_temp]
                       [--tmpdir] [--logdir] [--verbose] [--debug] [--supress]
                       [query ...]

CanSNPer2 - type query assemblies (or reads) against a CanSNP database

options:
  -h, --help            show this help message and exit

Required arguments:
  query                 File(s) to align (fasta), or one list file with a file
                        path per line (a file that does not start with '>'; a
                        list must be the only query argument)
  -db , --database      CanSNP database

Output options:
  -o DIR, --outdir DIR  Output directory
  --save_tree           Save the tree of each sample as <sample>_tree.html
                        (default False)
  --no_snpfiles         Do not write the result files (<sample>_snps.txt,
                        <sample>_not_called.txt); the call is still printed
  --summary             Output a summary file and tree with all called SNPs
                        not affected by no_snpfiles
  --annotations FILE    TSV file annotating tree nodes: <snp name>
                        <node|clade> <rank> <label>[ <color>]. 'node' draws a
                        badge under the SNP dot, 'clade' a dashed box around
                        its whole subtree; rank is a level carried by the
                        annotation (type, subspecies, marker...), never a tree
                        depth

Run options:
  --refdir              Specify reference directory
  --allow_missing_references
                        Type even when reference genomes of the database have
                        no <genome>.fna in --refdir; one warning names them,
                        and their SNPs are not compared (default: stop with an
                        error)
  --workdir             Change workdir default (./)
  --read_input          Select if input is reads not fasta
  --min_required_hits MIN_REQUIRED_HITS
                        Minimum number of derived SNPs on the path from the
                        root to the final SNP for a call (at least 1, default
                        3)
  --strictness STRICTNESS
                        Fraction of the SNPs on the path to the final SNP that
                        must be called; SNPs without alignment data count
                        against it, an ancestral SNP stops the path (0-1,
                        default 0.7)
  --mask MASK           Ignore SNPs closer than this many bases to the edge of
                        an alignment (default 0, off)
  --keep_going          If Error occurs, continue with the rest of samples
  --rerun               Rerun already processed files (else skip if result
                        file exists)
  --aligner             Program that aligns the query against each reference:
                        minimap2, mummer4, progressiveMauve (default
                        progressiveMauve)
  --alignment_timeout MINUTES
                        Stop the alignments of a query that take longer than
                        this, in minutes (0 for no limit, default 60)
  --aligner_path        The aligner program, or the folder that holds it
                        (default: found in PATH)
  --skip_alignment, --skip_mauve
                        Do not align, read the alignments already in --tmpdir
                        (--skip_mauve is the old name)
  --keep_temp           keep temporary files

Logging and debug options:
  --tmpdir              Specify tmp directory (default: CanSNPer2 in the
                        system temporary folder, see TMPDIR; a folder in
                        --outdir is used when that cannot be written to)
  --logdir              Specify log directory
  --verbose             Verbose output
  --debug               Debug output
  --supress             Supress warnings
```

### Annotating the tree

`--annotations FILE` draws meaning that SNP calling itself does not know about onto every
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

Queries are FASTA files (plain text; gzip-compressed files must be decompressed first). A
file is a query list instead when its first non-blank line does not start with `>`: the
list holds one file path per line. This is decided by the content, not by the file name, so
a FASTA file named `x.txt` is typed as a FASTA file. A query list must be the only query
argument; combining it with other arguments stops the run with a message. An empty file, a
file without a `>` header and a file that has a header but no sequence stop the run (or the
sample, with `--keep_going`) with a message naming the file.

### Sample and sequence names

A sample is named after its file name without the extension (`a/x.fasta` gives `x`), and
its results are written as `<outdir>/x_snps.txt` and `<outdir>/x_not_called.txt`. Two
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

`CanSNPer2-download --source refseq` renames each sequence of the downloaded RefSeq assembly
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

Each sample gives these files in `--outdir`:

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
- `x_tree.html` (with `--save_tree`): the tree of the database painted with the states of the
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
  - `not called`: no SNP path met `--min_required_hits` and `--strictness`, no derived SNP
    was found, or two SNP paths competed at the same depth;
  - `failed`: the alignment or the processing of the sample failed;
  - `reused`: the result was read from an earlier run (`--rerun` computes it again).

  Samples that need a look are listed first. Selecting a sample paints the tree with its
  states and shows its SNP path, the reasons for its status, links to its result files and
  the labels that `--annotations` gives the final SNP. In the tree, derived SNPs are green,
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

### Upgrading from 2.2

Version 2.3.0 changes how a SNP position is stored. Up to 2.2 it was one number along all
the sequences of the reference joined in the order of the reference fasta file, which made
a result depend on that order. From 2.3.0 the database and the source files hold the
sequence and the position on it (`sequence`, `sequence_position`; the `position` column
is gone), and the result tables hold the columns listed above.

A database or SNP source file in the old format is refused:

```
The database old.db cannot be used. It stores SNP positions along the joined reference
sequences (the format before 2.3.0). Convert it with CanSNPer2-database --convert, or
download the current database with CanSNPer2-download --scheme.
```

A database made by 2.3.0 cannot be read by 2.2.x, which stops with `no such column: position`.

`--convert` brings a SNP source file or a database of an older format to the current one, in
one step, and writes it to the folder given with `-o`. The input is not changed.

A file or database from before 2.3.0 has its positions located on the sequence of its reference
genome, which is read from `--refdir` (`<genome>.fna`; missing files are downloaded). A SNP whose
base in the reference is neither its derived nor its ancestral base is reported, and nothing is
written. A file or database of 2.3.0 only gains the rank, score and typed columns of the canSNP
format, so it needs no reference files: every canSNP becomes rank 1 of its node, typed, without
a score. A converted SNP file is written as the canSNP file of its scheme
(`<name>_snp.txt` becomes `<name>_canSNPs.txt`).

```sh
# a SNP source file from before 2.3.0; --references is optional and names the genomes of the file
CanSNPer2-database --convert --annotation old_snp.txt --references old_references.txt \
    --refdir references/ -o converted/

# a database from before 2.3.0
CanSNPer2-database --convert --database old.db --refdir references/ -o converted/

# a SNP file or database of 2.3.0: no reference files are needed
CanSNPer2-database --convert --annotation my_scheme_snp.txt -o converted/
CanSNPer2-database --convert --database my_scheme.db -o converted/
```

The reference files must have their sequences in the order in which the old positions were
counted (the order of the NCBI assembly file). The sequence names written are the first
words of their fasta headers; missing reference files are downloaded from GenBank, so the
names are GenBank accessions. progressiveMauve alignments made by 2.2 must be made again,
as `--skip_alignment` needs the sequence table (`<alignment>.sequences`) that 2.3.0 writes
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

A database is never left half changed: `--create`, `--import` and `CanSNPer2-download
--scheme` remove a database they had started, and a modification (`--mod_file`,
`--annotation`, `--references`, `--remove`) is made on a copy that replaces the database
only when every step succeeded.

| Message (beginning) | Command | What to do |
|---|---|---|
| `The source files cannot be built into a database, N problem(s):` followed by `<file>:<line>: ...` | `CanSNPer2-database --create`, `--import`, `CanSNPer2-download --scheme` | Correct the listed rows; all problems are listed together. |
| `<file>:<line>: genbank_id '...' is not an assembly accession` | as above, and `--references` | Give the GenBank (GCA_) or RefSeq (GCF_) assembly accession with its version, such as GCA_000008985.1. |
| `<file>:<line>: genbank_id '...' of genome '...' is already used by genome '...'` | as above, and `--references` on an existing database | Each genome needs its own assembly accession, as one reference file is written per accession. |
| `The files cannot be added to <db>, N problem(s):` | `CanSNPer2-database --mod_file` with `--annotation` or `--references` | A canSNP on a node that is not in the tree of the database or the modification file, on a genome that is not listed, a short row, or a genome, a rank of its node, a second typed canSNP of its node, or a position the database already holds. Correct the files; the database was not changed. |
| `The selected parent node (X) could not be found`, `The node X could not be found`, `X is a root of the tree` | `CanSNPer2-database --mod_file`, `--remove` | Give a node of the database with `--parent`; the root cannot be removed. |
| `The SNP X could not be written to <db>`, `The node X could not be written to <db>` | `CanSNPer2-database` | The disk is full or the database cannot be written to; free space or check the permissions, then run again. |
| `scheme X on owner/name@ref uses positions on the joined reference (before 2.3.0)` | `CanSNPer2-download --scheme` | The branch holds the format before 2.3.0. Use CanSNPer2 2.2 for it, or a branch converted for 2.3.0. |
| `... manifest.json was not found ...` or `... records no hash for ...`, `the files of the verified scheme X cannot be checked` | `CanSNPer2-download --scheme` | The files cannot be checked against the manifest of the repository. Report it to the repository, or give `--allow_unverified` when the source is trusted. |
| `The database X could not be built from the scheme Y`, `The database X could not be written` | `CanSNPer2-download --scheme` | An sqlite error or a file system error; free space or check the permissions of the folder. No partial file is left behind. |
| `X lists no reference genomes, there is nothing to download` | `CanSNPer2-download --database` | Build the database again with its references file. |
| `REF: the database has no genbank_id for this genome` | `CanSNPer2-download --database` | Give the assembly accession of the genome in the references file and build the database again. |
| `The SNP positions of X could not be read` | `CanSNPer2-download --database`, `CanSNPer2` | The database is damaged; export it with `CanSNPer2-database --export`, correct the files and build it again. |
| `N reference genome(s) of X have no file <genome>.fna in DIR: ...` | `CanSNPer2` | Fetch them with `CanSNPer2-download --database X --outdir DIR`. Without them the SNPs on these genomes would not be compared. To type without them on purpose, for example to test a scheme with a reference genome left out, give `--allow_missing_references`: the run then continues with one warning that names the missing genomes. It still stops when no genome of the database has a file. |
| `The database X has no canSNPs`, `The database X lists no reference genomes` | `CanSNPer2` | Build the database again with its canSNP file and its references file. |
| `The database X has N canSNP(s) but none of them is marked typed` | `CanSNPer2` | Every node needs one canSNP marked `typed yes`; build the database again from a canSNP file that marks them. |
| `The database X has N SNP(s) on node ids that are not in the nodes or tree table` | `CanSNPer2` | Export the database, correct the tree or canSNP file and build it again. Nodes without canSNPs are allowed. |

These cases give a warning and keep exit status 0: a draft scheme or a scheme of unknown
status; a reference whose base matches neither the ancestral nor the derived base at fewer
than a fifth of its SNPs; reference files in `--refdir` that are not reference genomes with
SNPs in the database; reference genomes without a file in `--refdir` when
`--allow_missing_references` is given; a database changed after it was built; `--no_references`; and a draft scheme, or a scheme read with
`--data_dir`, whose files cannot be checked against a manifest.

## Aligners

CanSNPer2 aligns each query genome against every reference genome of the database
and reads the query base at each SNP position from the alignment. Three aligners are
available:

* `progressiveMauve` (default), which writes xmfa alignments. The xmfa counts reference
  positions along the reference sequences joined end to end, so the names and lengths of
  those sequences are written next to each alignment (`<alignment>.sequences`), and
  `--skip_alignment` reads them again. This table is written for progressiveMauve only.
* `minimap2`, run as `minimap2 -c --cs=long -x asm20 -r 1k,10k -N 50 -p 0.1` and
  written as PAF. Only primary alignments with a mapping quality of at least 20 are
  read: a secondary alignment, or a primary one with a lower mapping quality, can carry
  the base of another copy of a repeated region. The long-join bandwidth (`-r`) and the
  secondary chain settings (`-N`, `-p`) differ from the preset so that the short
  stretches between IS elements and rearrangements, which separate Francisella genomes
  of different subspecies, are aligned as well. PAF gives positions on each reference
  sequence, so the alignment is read again with `--skip_alignment` without any other file.
* `mummer4`, MUMmer4 `nucmer` with its default settings and one thread (`-t 1`) followed
  by `delta-filter -1`, written as delta. With the default of two threads, repeated runs on
  the same genomes could give different alignments after the filter, and so a different
  base or no data at a few canSNP positions; one thread gives the same result every time. The 1-to-1 filter keeps, for every reference and query position,
  only the best alignment, allowing rearrangements. The delta file names the reference
  and query fasta files it was made from. The reference coordinates are positions on each
  reference sequence; the reader takes only the query bases from the query file, so that
  file must still be in place when the alignment is read again with `--skip_alignment`,
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
