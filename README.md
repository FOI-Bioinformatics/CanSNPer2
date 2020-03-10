# CanSNPer2
CanSNPer2: A toolkit for SNP-typing NGS data. The second release of CanSNPer (CanSNPer2) is exclusively written for python3. CanSNPer2 is simplified from CanSNPer1 stripped to only perform required tasks. It is also written in a modular form with separate classes for each task which allows future extentions such as a sequence read input option planned during 2020. 

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
Run analysis using canSNPer2

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
