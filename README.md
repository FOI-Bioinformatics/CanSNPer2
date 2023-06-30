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

* ETE3
* FlexTaxD - https://github.com/FOI-Bioinformatics/flextaxd
* progressiveMauve

Read input (not yet implemeted)
* BWA
* Minimap2

## User guide CanSNPer2 (for custom databases see below)
1. Download pre-built databases from https://github.com/FOI-Bioinformatics/CanSNPer2-data

2. Download references for database

```sh
CanSNPer2-download --database downloaded_database.db
```

3. Run genomes
```sh
CanSNPer2 --database downloaded_database.db fastadir/*.fasta --summary
```

The summary parameter will give a final result file with all SNPs that could be confirmed as well as a final CanSNP tree pdf with all SNPs from set colored.

For more options CanSNPer2 --help

## Quick start custom databases
Create database, download references annotated in the database and run CanSNPer2
```sh
CanSNPer2-database --database francisella_tularensis.db --annotation snps.txt --tree tree.txt --reference references.txt --source_type CanSNPer --create
CanSNPer2-download --database francisella_tularensis.db
CanSNPer2 sample1.fasta sample2.fasta --database francisella_tularensis.db --save_tree
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
usage: CanSNPer2 [-h] [-db] [-o DIR] [--save_tree] [--no_export]
                 [--refdir] [--workdir]
                 [--read_input] [--skip_mauve]
                 [--keep_going] [--keep_temp]
                 [--tmpdir] [--logdir] [--verbose] [--debug] [--supress]
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
  --no_export           no file output (can be used if summary only is requested)
  --summary             Output a summary file and tree with all called SNPs

Run options:
  --refdir              Specify reference directory
  --workdir             Change workdir default (./)
  --read_input          Select if input is reads not fasta
                        (not implemeted expected for version v2.1.0)
  --min_required_hits   MIN_REQUIRED_HITS
                        Minimum sequential hits to call a SNP!
  --keep_going          If Error occurs, continue with the rest of samples
  --rerun               Rerun already processed files (else skip if result file exists)

  --keep_temp           keep temporary files
  --skip_mauve          If xmfa files already exists skip step

Logging and debug options:
  --tmpdir              Specify reference directory
  --logdir              Specify log directory
  --verbose             Verbose output
  --debug               Debug output
  --supress             Supress warnings
```


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
