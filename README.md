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

### Create a database

```
CanSNPer2-create --database <path_to><database_name> --treefile <CanSNPer[1/2]_tree file> --snpfile <CanSNPer[1/2]_snpfile> --genomefile <CanSNPer2_genomesfile>
```

### Download references

```
CanSNPer2-download --database <path_to><database_name> --source genbank --outdir myworkdir
```

### Run canSNPer

```
CanSNPer2 --database <path_to><database_name> *.fa
```

## CanSNPer2 source file format.

1) Genome file

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


2) Tree file

The CanSNPer2 tree file must contain two columns with child and parent

| header | row |
| --- | --- |
| child               |   B.6 |
| parent              |   B.2 |


3) snp file

The CanSNPer2 snpfile contains all four columns defining a SNP as well as two annotations
to the published paper (reference) and the genome_id specified in the genome file.
In difference to canSNPer 1 the species column is not required
as different species are expected to be placed in different databases,
(or specified with the --organism parameter might be depricated)

| header | row |
| --- | --- |
| snp_id              |   B.2 |
| reference            |   Birdsell2014 |
| genome_id           | FSC200 |
| position            |   \<pos\> |
| ancestral_base      |   A |
| derived_base        |   T |


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
