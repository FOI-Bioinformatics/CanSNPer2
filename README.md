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
commit and the hash of each file it was built from. `--branch dev` gives access to draft
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

2. Alternatively, with a database you already have, download its references

```sh
CanSNPer2-download --database my_database.db
```

3. Run genomes
```sh
CanSNPer2 --database francisella_tularensis.db fastadir/*.fasta --summary
```

The summary parameter will give a final result file with all SNPs that could be confirmed as well as a final CanSNP tree (self-contained html, opens in any browser) with all SNPs from set colored.

For more options CanSNPer2 --help

## Quick start custom databases
Create database, download references annotated in the database and run CanSNPer2
```sh
CanSNPer2-database --database francisella_tularensis.db --annotation snps.txt --tree tree.txt --reference references.txt --source_type CanSNPer --create
CanSNPer2-download --database francisella_tularensis.db
CanSNPer2 sample1.fasta sample2.fasta --database francisella_tularensis.db --save_tree
```

### CanSNPer2-database

```
usage: CanSNPer2-database [-h] [-db] [--tree] [--annotation] [--references]
                      [--source_type] [--create] [--import] [--import_db]
                      [--mod_file] [--parent]
                      [--remove] [--replace] [--export] [--export_format] [-o]
                      [--tmpdir] [--logs] [--verbose] [--debug] [--supress]

CanSNPer2-database

optional arguments:
  -h, --help         show this help message and exit

Required:
  -db , --database   CanSNPer2 database name

Load database:
  --tree             CanSNPer tree source file
  --annotation       CanSNPer snp source file
  --references       File containing all reference genomes listed
  --source_type      Select source file type
  --create           Create new database!
  --import           Import an existing (legacy) database into the current format,
                     written to --database
  --import_db        Database to import from, with --import

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
```

#### Export

`--export` writes the three source files that `--create` reads, in the format of the
[CanSNPer2-data](https://github.com/FOI-Bioinformatics/CanSNPer2-data) source folder:

| File | Content |
|---|---|
| `<db>_tree.txt` | one line per node, the tab separated path from the root to the node |
| `<db>_snp.txt` | `snp_id reference genome position derived_base ancestral_base date` |
| `<db>_references.txt` | `genome strain genbank_id refseq_id assembly_name` |

A database created from these files has the same tree, SNPs and reference genomes, so a scheme
can be exported, edited as text and rebuilt:

```sh
CanSNPer2-database --database my_scheme.db --export -o scheme/
CanSNPer2-database --database rebuilt.db --create --tree scheme/my_scheme_tree.txt \
    --annotation scheme/my_scheme_snp.txt --references scheme/my_scheme_references.txt
```

A database whose tree has a node with more than one parent, more than one root, or SNPs that
point to a missing reference genome cannot be written in this format and is refused.
`--export_format newick` writes the tree only. Up to version 2.2.2 the default export (`tab`)
wrote node and link tables that could not rebuild the database; `tab` now gives the source files.

#### Import

`--import` reads an existing database and writes a new one in the current format, the way to
upgrade databases made by older releases. Legacy databases carry an extra orphan root left by
the old tree reader (issues #32/#34/#149): two self-linked nodes, which `--export` refuses.
`--import` drops the orphan and produces a single-rooted database that then exports and
rebuilds normally. The source is never modified.

```sh
CanSNPer2-database --import --import_db old_scheme.db --database new_scheme.db
```

A legacy shipped database (e.g. `francisella_tularensis.db`) imports to one root with all its
nodes, SNPs and reference genomes kept, and is then exportable; re-importing a database that
is already current changes nothing.

### CanSNPer2-download

```
usage: CanSNPer2-download [-h] [-db] [--scheme] [--branch] [--repo]
                          [--data_dir] [--list_schemes] [--overwrite]
                          [--no_references] [-s] [-o] [--logs] [--verbose]

CanSNPer2-download

options:
  -h, --help         show this help message and exit

Database:
  -db , --database   CanSNP database (with --scheme: the database to build,
                     default <scheme>.db)

CanSNPer2-data scheme options:
  --scheme           Fetch this scheme from CanSNPer2-data and build its
                     database
  --branch           Branch, tag or commit of the repository; in
                     CanSNPer2-data master holds verified schemes, dev also
                     draft schemes (default master, with --repo the default
                     branch of that repository)
  --repo             GitHub repository laid out as CanSNPer2-data, as
                     owner/name (default FOI-Bioinformatics/CanSNPer2-data)
  --data_dir         Local directory laid out as CanSNPer2-data (a clone, or a
                     folder of custom schemes) to read schemes from instead of
                     GitHub
  --list_schemes     List the schemes of the repository or --data_dir and exit
  --overwrite        Replace an existing database when building from --scheme
  --no_references    Build the database from --scheme without downloading its
                     references

Download options:
  -s , --source      Source for download (genbank/refseq)
  -o , --outdir      reference genomes folder

Logging and debug options:
  --logs             Specify log directory
  --verbose          Verbose logging

Give a database with --database, or name a scheme with --scheme to fetch it
from CanSNPer2-data (or --repo, or --data_dir), build its database and
download its references in one step.
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
Example structure of snps.txt (NEW headerline)
```
snp_id	strain	reference	genome	position	derived_base	ancestral_base
T/N.1	francisella	Svensson2009	SCHUS4.2	83976	A	G
T.1	francisella	Svensson2009	SCHUS4.2	1165688	G	A
M.1	francisella	Vogler2009	SCHUS4.2	75124	T	C
A/M.1	francisella	Vogler2009	SCHUS4.2	1491914	A	G
A.1	francisella	Vogler2009	SCHUS4.2	397639	T	C
B.1	francisella	Birdsell2014	OSU18	1710718	T	C
...
```
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
               [--summary] [--annotations FILE] [--refdir] [--workdir] [--read_input]
               [--min_required_hits MIN_REQUIRED_HITS]
               [--strictness STRICTNESS] [--mask MASK] [--keep_going]
               [--rerun] [--aligner] [--alignment_timeout MINUTES]
               [--aligner_path] [--skip_alignment]
               [--keep_temp] [--tmpdir]
               [--logdir] [--verbose] [--debug] [--supress]
               [query [query ...]]

CanSNPer2

optional arguments:
  -h, --help            show this help message and exit

Required arguments:
  query                 File(s) to align (fasta)
  -db , --database      CanSNP database

Output options:
  -o DIR, --outdir DIR  Output directory
  --save_tree           Save tree as PDF using ETE3 (default False)
  --no_snpfiles         Don´t save output files.
  --summary             Output a summary file and tree with all called SNPs
                        not affected by no_snpfiles
  --annotations FILE    TSV file annotating tree nodes: <snp name>  <node|clade>
                        <rank> <label> [<color>]. 'node' draws a badge under the
                        SNP dot, 'clade' a dashed box around its whole subtree;
                        rank is a level carried by the annotation (type,
                        subspecies, marker...), never a tree depth

Run options:
  --refdir              Specify reference directory
  --workdir             Change workdir default (./)
  --read_input          Select if input is reads not fasta
  --min_required_hits MIN_REQUIRED_HITS
                        Minimum sequential hits to call a SNP!
  --strictness STRICTNESS
                        Fraction of the SNPs on the path to the final SNP
                        that must be called; SNPs without alignment data count
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
  --tmpdir              Specify tmp directory default (/tmp/CanSNPer2)

Each run that aligns writes its temporary files to its own directory,
`<tmpdir>/CanSNPer2/run_<random>`, and removes only that directory when it
finishes, so several runs can share the same `--tmpdir`. With `--keep_temp` the
directory is kept and its path is logged; the alignments in it can be reused
with `--skip_alignment --tmpdir <that directory>`.
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

### Sample and sequence names

A sample is named after its file name without the extension (`a/x.fasta` gives `x`), and
its results are written as `<outdir>/x_snps.txt` and `<outdir>/x_not_called.txt`. Two
different files with the same name in one run stop the run before anything is aligned;
rename one or type them into separate `--outdir`. The same file given twice is typed once.

`x_snps.txt` ends with a `Query file:` line that gives the input it was made from and the
SHA-256 of that file. A later run uses an existing result again (without `--rerun`) only
when the input has the same SHA-256; a different input with the same name is not typed and
the reason is reported. Result files from versions before this line was added are used
again as before.

Sequences are identified by the first word of their header line, in references and in
queries; these names must be set and unique, otherwise the run (reference) or the sample
(query) stops with a message naming the file.


## Aligners

CanSNPer2 aligns each query genome against every reference genome of the database
and reads the query base at each SNP position from the alignment. Three aligners are
available:

* `progressiveMauve` (default), which writes xmfa alignments.
* `minimap2`, run as `minimap2 -c --cs=long -x asm20 -r 1k,10k -N 50 -p 0.1` and
  written as PAF. Only primary alignments with a mapping quality of at least 20 are
  read: a secondary alignment, or a primary one with a lower mapping quality, can carry
  the base of another copy of a repeated region. The long-join bandwidth (`-r`) and the
  secondary chain settings (`-N`, `-p`) differ from the preset so that the short
  stretches between IS elements and rearrangements, which separate Francisella genomes
  of different subspecies, are aligned as well. The reference sequence table
  (`<alignment>.sequences`) is written next to each alignment, so `--skip_alignment`
  can read it again.
* `mummer4`, MUMmer4 `nucmer` with its default settings followed by `delta-filter -1`,
  written as delta. The 1-to-1 filter keeps, for every reference and query position,
  only the best alignment, allowing rearrangements. The delta file names the reference
  and query fasta files it was made from; the reader takes the query bases and the
  reference sequence order from them, so those files must still be in place when the
  alignment is read again with `--skip_alignment`. `nucmer` and `delta-filter` must
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
database: for a reference with several sequences (chromosomes, plasmids) the sequences
are joined end to end in the order of the reference fasta file, as progressiveMauve
does. A new aligner is registered in `ALIGNERS` and can then be selected with `--aligner`.

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
