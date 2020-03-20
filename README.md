# CanSNPer2
The second release of CanSNPer (CanSNPer2) is exclusively written for python3. CanSNPer2 is simplified from CanSNPer1 stripped to only perform required tasks. It is also written in a modular form with separate classes for each task which allows future extentions such as a sequence read input option planned during 2020. 

## Installation
Installing using bioconda

```
conda install cansnper2
```

Installing from the repository
```
python setup.py install
```

## Requirements

* ETE3
* FlexTaxD
* progressiveMauve

## Quick start
Create database, download references and run CanSNPer2
```sh
CanSNPer2-database --database francisella_tularensis.db --annotation snps.txt --tree tree.txt --reference references.txt --source_type CanSNPer --create
CanSNPer2-download --database francisella_tularensis.db -o references
CanSNPer2 sample.fasta -db francisella_tularensis.db -o results --save_tree --refdir references --snpfile 200320
```
references.txt
```
#genome strain  genbank_id      refseq_id       assembly_name
OSU18   OSU18   GCA_000014605.1 GCF_000014605.1 ASM1460v1
FSC200  FSC200  GCA_000168775.2 GCF_000168775.2 ASM16877v2
FTNF002-00      FTNF002-00      GCA_000017785.1 GCF_000017785.1 ASM1778v1
LVS     LVS     GCA_000009245.1 GCF_000009245.1 ASM924v1
SCHUS4.2        SCHUS4  GCA_000008985.1 GCF_000008985.1 ASM898v1

```
snps.txt (NEW headerline)
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
tree.txt
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
CanSNPer2 
    --database <path_to><database_name>               ## CanSNPer2 database file
    --refdir test/references/                         ## directory where references are located if not present see CanSNPer2-download
    --export (export snpfile)                         ## output snpfile
    --snpfile snptest2.txt                            ## name of snp <tab delim> output
    -o results                                        ## output directory
    --save_tree <prefix>                              ## prefix for CanSNPer2 tree output (pdf)
    --tmpdir /tmp/                                    ## Select tempdir default global tmp
    *.fa                                              ## Any number of fasta files to run
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
