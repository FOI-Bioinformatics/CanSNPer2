# Command line reference

This page lists the help of every `cansnper2` command, as printed by
`cansnper2 <command> --help`. It is generated from the program by
`tools/gen_cli_docs.py`; do not edit it by hand.

Options that take a value show their default in parentheses. Paths are relative
to the current directory, except for `cansnper2 type`, where the output, log and
temporary folders are relative to `--workdir`.

| Command | Purpose |
| --- | --- |
| [`cansnper2 type`](#cansnper2-type) | Type assemblies or sequencing reads against a scheme database |
| [`cansnper2 selftest`](#cansnper2-selftest) | Check the installation: Python modules and external programs |
| [`cansnper2 scheme list`](#cansnper2-scheme-list) | List the schemes and their status |
| [`cansnper2 scheme get`](#cansnper2-scheme-get) | Fetch a scheme, build its database and download its references |
| [`cansnper2 ref get`](#cansnper2-ref-get) | Download the reference genomes of a database |
| [`cansnper2 db create`](#cansnper2-db-create) | Build a database from the files of a scheme |
| [`cansnper2 db import`](#cansnper2-db-import) | Build a database from the files written by db export |
| [`cansnper2 db export`](#cansnper2-db-export) | Export a database to text files |
| [`cansnper2 db info`](#cansnper2-db-info) | Print the scheme information and size of a database |
| [`cansnper2 db convert`](#cansnper2-db-convert) | Convert a database or canSNP file from an earlier format |
| [`cansnper2 db modify`](#cansnper2-db-modify) | Change the tree, canSNPs or reference genomes of a database |
| [`cansnper2 db remove`](#cansnper2-db-remove) | Remove the branch rooted at a node |

## cansnper2

```text
usage: cansnper2 [-h] [--version] <command> ...

CanSNPer2 assigns bacterial genomes to lineages of a canonical SNP
(canSNP) typing scheme, from assemblies or sequencing reads. It also
builds and maintains scheme databases and fetches published schemes
and their reference genomes.

positional arguments:
  <command>
    type      type assemblies or sequencing reads against a scheme database
    selftest  check the installation: Python modules and external programs
    db        create, import, export, info, convert, modify, remove
    scheme    list: the published schemes; get: build the database of one
    ref       get: download the reference genomes of a database from NCBI

options:
  -h, --help  show this help message and exit
  --version   show program's version number and exit

Run cansnper2 <command> --help for the options of a command.
Before 2.6.0 these jobs were done by the commands CanSNPer2,
CanSNPer2-database, CanSNPer2-download and CanSNPer2-test, which are
no longer installed; see docs/upgrading.md.
Documentation: https://github.com/FOI-Bioinformatics/CanSNPer2
```

## cansnper2 type

```text
usage: cansnper2 type [-h] [-d DB] [-o DIR] [--save-tree] [--no-snp-files]
                      [--summary] [--node-labels FILE] [--ref-dir DIR]
                      [--allow-missing-references] [--workdir DIR]
                      [--min-required-hits N] [--strictness FRACTION]
                      [--mask BASES] [--keep-going] [--rerun]
                      [--aligner {minimap2,mummer4,progressiveMauve}]
                      [--alignment-timeout MINUTES] [--aligner-path PATH]
                      [--skip-alignment] [--keep-temp] [--tmp-dir DIR]
                      [--platform {hifi,illumina,ont}] [--samplesheet FILE]
                      [--min-depth N] [--min-allele-fraction FRACTION]
                      [--min-mapq Q] [--min-baseq Q] [--mixed-sites N]
                      [--threads N] [--max-reads N] [--local-references]
                      [--local-rounds N] [--screen-reads N] [--fetch-local]
                      [--log-dir DIR] [-v] [--debug] [--quiet]
                      [query ...]

Type assemblies or sequencing reads against a scheme database: each query is
compared with the reference genomes of the database at the canSNP positions,
and the deepest node of the scheme tree that the derived canSNPs support is
reported.

options:
  -h, --help            show this help message and exit

Input:
  query                 assemblies (FASTA) or sequencing reads (FASTQ, may be
                        gzip-compressed) to type, or one list file with one
                        file path per line; paired read files are recognized
                        by their names (_R1/_R2, _1/_2)
  -d, --database DB     the scheme database (.db); required

Output options:
  -o, --outdir DIR      output directory; the overview summary.html is written
                        there and the per-sample files in its samples folder
                        (default: results)
  --save-tree           also write a standalone tree page per sample
                        (samples/<sample>_tree.html, the overview already
                        shows every tree)
  --no-snp-files        do not write the per-sample canSNP files; the call is
                        still printed and written to the run summary
  --summary             also write snps_summary.txt and a tree with the final
                        SNPs of all samples
  --node-labels FILE    tab-separated file of labels for tree nodes in the
                        HTML pages: node name, kind (node or clade), rank,
                        label[, color]

Run options:
  --ref-dir DIR         folder with the reference genomes of the database
                        (<genome>.fna) (default: references/)
  --allow-missing-references
                        type even when reference genomes have no file in
                        --ref-dir
  --workdir DIR         directory that relative paths of the output, log and
                        temporary folders are resolved against (default: ./)
  --min-required-hits N
                        minimum number of derived canSNPs on the path from the
                        root to a node for that node to be called (default: 3)
  --strictness FRACTION
                        minimum fraction of the canSNPs on that path that must
                        be derived; the others may lack data (0-1) (default:
                        0.7)
  --mask BASES          assemblies only: ignore canSNPs closer than this many
                        bases to the edge of an alignment block (default: 0)
  --keep-going          continue with the other queries when one cannot be
                        typed
  --rerun               type queries again although their result file exists
  --aligner {minimap2,mummer4,progressiveMauve}
                        program that aligns an assembly to each reference
                        genome (default: minimap2)
  --alignment-timeout MINUTES
                        stop an alignment that runs longer than this; 0 for no
                        limit (default: 60)
  --aligner-path PATH   the aligner program, or the folder that holds it
                        (default: searched in PATH)
  --skip-alignment      do not align; read the alignments kept in --tmp-dir by
                        an earlier run with --keep-temp
  --keep-temp           keep the temporary alignment files
  --tmp-dir DIR         folder for temporary alignments (default: CanSNPer2 in
                        the system temporary folder given by TMPDIR, or
                        <outdir>/.cansnper2_tmp when that cannot be written)

Read typing options:
  used when the queries are sequencing reads (minimap2 and samtools are
  required); thresholds that are not given take the defaults of --platform

  --platform {hifi,illumina,ont}
                        sequencing platform: illumina (short reads, paired or
                        single), ont (Oxford Nanopore) or hifi (PacBio HiFi)
                        (default: illumina)
  --samplesheet FILE    tab-separated file with the header columns sample,
                        reads1 and optionally reads2 and platform, one sample
                        per row; replaces the query arguments
  --min-depth N         minimum number of A, C, G or T bases at a canSNP
                        position, after the mapping and base quality filters,
                        for its state to be read (default: 5 for all
                        platforms)
  --min-allele-fraction FRACTION
                        minimum fraction of the reads that must carry one base
                        for the site to take that state; below it the site is
                        reported as mixed (default: illumina 0.8, ont 0.7,
                        hifi 0.8)
  --min-mapq Q          minimum mapping quality of a read (default: 20 for all
                        platforms)
  --min-baseq Q         minimum base quality (default: illumina 20, ont 10,
                        hifi 20)
  --mixed-sites N       number of mixed canSNP sites from which a sample is
                        flagged as possibly mixed (default: 3)
  --threads N           threads for minimap2 and samtools sort, per sample;
                        they are divided between the reference genomes the
                        reads are mapped to at the same time, at least one
                        each (default: 4)
  --max-reads N         number of reads mapped per sample, the first ones in
                        each read file; read pairs for paired-end reads; 0
                        maps all reads (default: 300000 for all platforms)
  --local-references    map the reads again to the representative genome of
                        the scheme closest to the sample and read the canSNPs
                        at their positions in that genome (genomes in <ref-
                        dir>/local; see ref get --local-references)
  --local-rounds N      with --local-references: maximum number of rounds of
                        choosing a closer representative genome (default: 2)
  --screen-reads N      with --local-references: number of reads per file that
                        are mapped to the candidate genomes to choose the
                        local reference (default: 100000)
  --fetch-local         with --local-references: download from NCBI during
                        typing the representative genomes that are needed and
                        not yet in <ref-dir>/local

Logging options:
  --log-dir DIR         folder for the log files (default: logs/)
  -v, --verbose         print progress information
  --debug               print debug information, and the traceback of an
                        unexpected error
  --quiet               print errors only, not warnings

examples:
  cansnper2 type sample1.fna sample2.fna -d francisella_tularensis.db
  cansnper2 type queries.txt -d scheme.db --ref-dir references/ -o results/
  cansnper2 type s1_R1.fastq.gz s1_R2.fastq.gz -d scheme.db \
      --platform illumina
  cansnper2 type --samplesheet samples.tsv -d scheme.db
```

## cansnper2 selftest

```text
usage: cansnper2 selftest [-h] [--debug]

Check that the Python modules can be imported and report which aligners and
read typing programs are found in PATH.

options:
  -h, --help  show this help message and exit
  --debug     print the traceback of an unexpected error
```

## cansnper2 scheme list

```text
usage: cansnper2 scheme list [-h] [--repo OWNER/NAME] [--branch REF]
                             [--data-dir DIR]

List the schemes of CanSNPer2-data (or of another repository or local
directory) with their status. This verb writes no log file.

options:
  -h, --help         show this help message and exit

Scheme source options:
  --repo OWNER/NAME  GitHub repository laid out as CanSNPer2-data (default:
                     FOI-Bioinformatics/CanSNPer2-data)
  --branch REF       branch, tag or commit of the repository; in
                     CanSNPer2-data, master holds the verified schemes and dev
                     also the draft schemes (default: master, or the default
                     branch of the repository given with --repo)
  --data-dir DIR     local directory laid out as CanSNPer2-data to read the
                     schemes from instead of GitHub; cannot be combined with
                     --repo or --branch

examples:
  cansnper2 scheme list
  cansnper2 scheme list --branch dev
```

## cansnper2 scheme get

```text
usage: cansnper2 scheme get [-h] [-d DB] [--repo OWNER/NAME] [--branch REF]
                            [--data-dir DIR] [--overwrite] [--no-references]
                            [--allow-unverified] [--local-references]
                            [--log-dir DIR] [-v] [--debug] [--quiet]
                            NAME

Fetch the files of a scheme, check them, build the database and download its
reference genomes from NCBI GenBank into ./references. Use ref get to download
to another folder or from RefSeq.

positional arguments:
  NAME                name of the scheme (see scheme list)

options:
  -h, --help          show this help message and exit
  -d, --database DB   the database to build (default: <NAME>.db)
  --overwrite         replace an existing database
  --no-references     build the database without downloading its references
  --allow-unverified  build a verified scheme although its files cannot be
                      checked against the checksums of its manifest
  --local-references  also download the representative genomes of the scheme
                      into references/local, for cansnper2 type --local-
                      references

Scheme source options:
  --repo OWNER/NAME   GitHub repository laid out as CanSNPer2-data (default:
                      FOI-Bioinformatics/CanSNPer2-data)
  --branch REF        branch, tag or commit of the repository; in
                      CanSNPer2-data, master holds the verified schemes and
                      dev also the draft schemes (default: master, or the
                      default branch of the repository given with --repo)
  --data-dir DIR      local directory laid out as CanSNPer2-data to read the
                      schemes from instead of GitHub; cannot be combined with
                      --repo or --branch

Logging options:
  --log-dir DIR       folder for the log files (default: logs)
  -v, --verbose       print progress information
  --debug             print debug information, and the traceback of an
                      unexpected error
  --quiet             print errors only, not warnings

examples:
  cansnper2 scheme get francisella_tularensis
  cansnper2 scheme get brucella_genus --branch dev
```

## cansnper2 ref get

```text
usage: cansnper2 ref get [-h] [-d DB] [-s {genbank,refseq}] [-o DIR]
                         [--local-references] [--log-dir DIR] [-v] [--debug]
                         [--quiet]

Download the reference genomes a database names from NCBI into a folder, as
<genome>.fna.

options:
  -h, --help            show this help message and exit
  -d, --database DB     the database whose reference genomes are downloaded;
                        required
  -s, --source {genbank,refseq}
                        NCBI assembly database to download from (default:
                        genbank)
  -o, --outdir DIR      folder for the reference genomes (default: references)
  --local-references    also download the representative genomes of the scheme
                        into <outdir>/local, for cansnper2 type --local-
                        references

Logging options:
  --log-dir DIR         folder for the log files (default: logs)
  -v, --verbose         print progress information
  --debug               print debug information, and the traceback of an
                        unexpected error
  --quiet               print errors only, not warnings

examples:
  cansnper2 ref get -d francisella_tularensis.db
  cansnper2 ref get -d scheme.db -s refseq -o references/
```

## cansnper2 db create

```text
usage: cansnper2 db create [-h] [-d DB] [--tree FILE] [--cansnps FILE]
                           [--references FILE] [--genomes FILE]
                           [--scheme-name NAME] [--scheme-status STATUS]
                           [--scheme-source SRC] [--log-dir DIR] [-v]
                           [--debug] [--quiet]

Build a database from the tree, canSNP and reference genome files of a scheme;
a database can also be started from only some of them and completed with db
modify. The files are checked first; the database is written only when they
pass.

options:
  -h, --help            show this help message and exit
  -d, --database DB     the new database to write; required
  --tree FILE           tree file of the scheme (<name>_tree.txt)
  --cansnps FILE        canSNP file of the scheme (<name>_canSNPs.txt)
  --references FILE     reference genomes file of the scheme
                        (<name>_references.txt)
  --genomes FILE        genomes file of the scheme (<name>_genomes.txt:
                        genome, node, role); optional, used by read typing
                        with --local-references

Scheme information (recorded in the database):
  --scheme-name NAME    name of the scheme, printed by db info and in the
                        typing results
  --scheme-status STATUS
                        verified, draft or unknown
  --scheme-source SRC   where the scheme files came from, e.g. a repository
                        and commit

Logging options:
  --log-dir DIR         folder for the log files (default: logs)
  -v, --verbose         print progress information
  --debug               print debug information, and the traceback of an
                        unexpected error
  --quiet               print errors only, not warnings

examples:
  cansnper2 db create -d my_scheme.db --tree my_tree.txt \
      --cansnps my_canSNPs.txt --references my_references.txt \
      --scheme-name my_scheme --scheme-status draft
```

## cansnper2 db import

```text
usage: cansnper2 db import [-h] [-d DB] [--tree FILE] [--cansnps FILE]
                           [--references FILE] [--genomes FILE]
                           [--scheme-name NAME] [--scheme-status STATUS]
                           [--scheme-source SRC] [--log-dir DIR] [-v]
                           [--debug] [--quiet]

Build a database from the files written by db export; the same checks and
options as db create.

options:
  -h, --help            show this help message and exit
  -d, --database DB     the new database to write; required
  --tree FILE           tree file of the scheme (<name>_tree.txt); required
  --cansnps FILE        canSNP file of the scheme (<name>_canSNPs.txt);
                        required
  --references FILE     reference genomes file of the scheme
                        (<name>_references.txt); required
  --genomes FILE        genomes file of the scheme (<name>_genomes.txt:
                        genome, node, role); optional, used by read typing
                        with --local-references

Scheme information (recorded in the database):
  --scheme-name NAME    name of the scheme, printed by db info and in the
                        typing results
  --scheme-status STATUS
                        verified, draft or unknown
  --scheme-source SRC   where the scheme files came from, e.g. a repository
                        and commit

Logging options:
  --log-dir DIR         folder for the log files (default: logs)
  -v, --verbose         print progress information
  --debug               print debug information, and the traceback of an
                        unexpected error
  --quiet               print errors only, not warnings

examples:
  cansnper2 db import -d new.db --tree export/scheme_tree.txt \
      --cansnps export/scheme_canSNPs.txt \
      --references export/scheme_references.txt
```

## cansnper2 db export

```text
usage: cansnper2 db export [-h] [-d DB] [--format {source,newick,tab}]
                           [-o DIR] [--log-dir DIR] [-v] [--debug] [--quiet]

Write a database as text files: the scheme files that db import reads, or the
tree in Newick format.

options:
  -h, --help            show this help message and exit
  -d, --database DB     the database; required
  --format {source,newick,tab}
                        source: the tree, canSNP, reference genome and (when
                        the database has them) genomes files that db import
                        reads; newick: the tree only, to
                        <outdir>/<database>_tree.nwk, or to standard output
                        without -o; tab: former name of source (default:
                        source)
  -o, --outdir DIR      folder for the exported files (default: the current
                        folder)

Logging options:
  --log-dir DIR         folder for the log files (default: logs)
  -v, --verbose         print progress information
  --debug               print debug information, and the traceback of an
                        unexpected error
  --quiet               print errors only, not warnings

examples:
  cansnper2 db export -d scheme.db -o export/
  cansnper2 db export -d scheme.db --format newick > scheme.nwk
```

## cansnper2 db info

```text
usage: cansnper2 db info [-h] [-d DB] [--log-dir DIR] [-v] [--debug] [--quiet]

Print the scheme information a database records (name, status, source, the
CanSNPer2 version that wrote it) and the number of nodes, canSNPs and
reference genomes.

options:
  -h, --help         show this help message and exit
  -d, --database DB  the database; required

Logging options:
  --log-dir DIR      folder for the log files (default: logs)
  -v, --verbose      print progress information
  --debug            print debug information, and the traceback of an
                     unexpected error
  --quiet            print errors only, not warnings

examples:
  cansnper2 db info -d francisella_tularensis.db
```

## cansnper2 db convert

```text
usage: cansnper2 db convert [-h] [-d DB] [--cansnps FILE] [--ref-dir DIR]
                            [-o DIR] [--log-dir DIR] [-v] [--debug] [--quiet]

Convert a database or canSNP file written by CanSNPer2 before 2.4.0 to the
current format. Files from before 2.3.0 place canSNPs on concatenated genomes;
their positions are moved to the reference sequences, which needs the
reference files (--ref-dir).

options:
  -h, --help         show this help message and exit
  -d, --database DB  the database to convert (or give --cansnps instead)
  --cansnps FILE     the canSNP source file to convert (or give --database
                     instead)
  --ref-dir DIR      folder with the reference files <genome>.fna, needed for
                     formats before 2.3.0; missing files are downloaded
  -o, --outdir DIR   folder for the converted file; required

Logging options:
  --log-dir DIR      folder for the log files (default: logs)
  -v, --verbose      print progress information
  --debug            print debug information, and the traceback of an
                     unexpected error
  --quiet            print errors only, not warnings

examples:
  cansnper2 db convert -d old.db --ref-dir references/ -o converted/
  cansnper2 db convert --cansnps old_snp.txt --ref-dir references/ \
      -o converted/
```

## cansnper2 db modify

```text
usage: cansnper2 db modify [-h] [-d DB] [--mod-file FILE] [--parent NODE]
                           [--replace] [--cansnps FILE] [--references FILE]
                           [--log-dir DIR] [-v] [--debug] [--quiet]

Change an existing database: add or replace part of the tree (--mod-file), or
add canSNPs (--cansnps) or reference genomes (--references). canSNPs and
reference genomes that the database already has are refused, not changed.

options:
  -h, --help         show this help message and exit
  -d, --database DB  the database; required
  --mod-file FILE    file with the tree modifications
  --parent NODE      node at which the tree is changed; required with --mod-
                     file
  --replace          replace the subtree below --parent with the tree of
                     --mod-file
  --cansnps FILE     canSNP file with canSNPs to add
  --references FILE  reference genomes file with genomes to add

Logging options:
  --log-dir DIR      folder for the log files (default: logs)
  -v, --verbose      print progress information
  --debug            print debug information, and the traceback of an
                     unexpected error
  --quiet            print errors only, not warnings

examples:
  cansnper2 db modify -d scheme.db --mod-file new_branch.txt --parent B.Br.001
  cansnper2 db modify -d scheme.db --mod-file subtree.txt \
      --parent B.Br.001 --replace
  cansnper2 db modify -d scheme.db --cansnps added_canSNPs.txt
```

## cansnper2 db remove

```text
usage: cansnper2 db remove [-h] [-d DB] [--parent NODE] [--log-dir DIR] [-v]
                           [--debug] [--quiet]

Remove a node and every node below it from a database.

options:
  -h, --help         show this help message and exit
  -d, --database DB  the database; required
  --parent NODE      root of the branch to remove; required

Logging options:
  --log-dir DIR      folder for the log files (default: logs)
  -v, --verbose      print progress information
  --debug            print debug information, and the traceback of an
                     unexpected error
  --quiet            print errors only, not warnings

examples:
  cansnper2 db remove -d scheme.db --parent B.Br.001
```
