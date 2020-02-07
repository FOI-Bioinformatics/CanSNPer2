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

## Usage

### Run canSNPer

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

### Modify database

```
CanSNPer2-database 
    --database <path_to><database_name>
    --mod_file <path_to><modification tree (1)>
    --annotation <path_to><modification annotation (2)>
    --log database/logs 
    --parent "A/M.1" 
    --replace
```

### Create a custom database

```
CanSNPer2-database 
    --database <path_to><database_name>               ## CanSNPer2 database file 
    --tree <CanSNPer[1/2]_tree file>                  ## Format (1) or cansnper1 formatted file
    --annotation <CanSNPer[1/2]_snpfile>              ## snp annotations (2)
    --reference <CanSNPer2_genomesfile>               ## Annotation of genome references to which snps are annotated (4)
    --source_type CanSNPer(if CanSNPer1)              ## If input is formatted according to cansnper1
    --create 
    --log <path to log folder>
```


### Download references

```
CanSNPer2-download 
    --database <path_to><database_name>               ## CanSNPer2 database file (must have genome annotation loaded)
    --source genbank                                  ## Source refseq or genbank
    --outdir references                               ## reference directory (default ./reference)
```

## CanSNPer2 source file format.


### Modify database

1) Tree file (for updating the database)

The CanSNPer2 tree file must contain two columns with child and parent

| header | row |
| --- | --- |
| child               |   B.6 |
| parent              |   B.2 |
| rank                |   species    |


2) snp annotation file

The CanSNPer2 snpfile contains all four columns defining a SNP as well as two annotations
to the published paper (reference) and the genome_id specified in the genome file.
In difference to canSNPer 1 the species column is not required
as different species are expected to be placed in different databases,
(or specified with the --organism parameter might be depricated)

| header | row |
| --- | --- |
| snp_id                |   B.2 |
| strain                |  francisella |
| reference             |   Birdsell2014 |
| genome                | FSC200 |
| position              |   \<pos\> |
| ancestral_base        |   A |
| derived_base          |   T |

### Create database

3) Genome file

The CanSNPer2 genomes file contains information about the reference sequences,
The file has 7 columns. However only the first column is nessesary for running CanSNPer,
But then the automatic download script won't work and the correct reference sequences
has to be supplied manually. Observe that the first row must contain headers

| header | row |
| ------- | -------- |
| genome_id          |   OSU18 |
| strain             |   OSU18 |
| refseq_id          |   GCA_000014605.1 |
| genbank_id         |   GCF_000014605.1 |
| assembly_name      |   ASM1460v1 |
| refseq_sequence    |   NC_008369.1 |
| genbank_sequence   |   CP000437.1 |

4) genome reference file
| header | row |
| --- | --- |
| genome                |   OSU18 |
| strain                |  OSU18 |
| genbank_id            |   GCA_000014605.1 |
| refseq_id             | GCF_000014605.1 |
| assembly_name         |   ASM1460v1 |

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
