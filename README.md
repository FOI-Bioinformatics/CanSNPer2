# CanSNPer2
CanSNPer2: A toolkit for SNP-typing NGS data. Copyright (C) 2019 David Sundell @ FOI bioinformatics group  VERSION 2.0.1 First release of CanSNPer2  The second release of CanSNPer (CanSNPer2) is exclusively written for python3 CanSNPer2 is simplified from CanSNPer1 stripped to only perform required tasks. It is also written in a modular form with separate classes for each task which allows future extentions such as a sequence read input option planned during 2020.  Requires progressiveMauve  This program is free software: you can redistribute and/or modify the code under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.  This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.  You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

## Installation

conda install CanSNPer2

python setup.py install

## Requirements
ETE3
FlexTaxD

## Usage

### Create a database
CanSNPer2-create --database <path_to><database_name> --treefile <CanSNPer[1/2]_tree file> --snpfile <CanSNPer[1/2]_snpfile> --genomefile <CanSNPer2_genomesfile>

### Download references
CanSNPer2-download --database <path_to><database_name> --source genbank --outdir myworkdir

### Run canSNPer
CanSNPer2 --database <path_to><database_name> *.fa


## CanSNPer2 source file format.

1)
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

</br>
2)
The CanSNPer2 tree file must contain two columns with child and parent

| header | row |
| --- | --- |
| child               |   B.6 |
| parent              |   B.2 |


3)
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
