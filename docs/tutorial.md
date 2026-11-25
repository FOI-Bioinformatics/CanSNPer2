# Tutorial

This tutorial types three public *Francisella tularensis* genome assemblies, and one public set of
Illumina reads, with the published `francisella_tularensis` scheme. It is written for new users.
It assumes CanSNPer2 is installed as
described in [installation.md](installation.md) and that the machine has network access to
GitHub and NCBI.

The output shown on this page was produced on 25 November 2026 with CanSNPer2 2.6.0, minimap2 2.30
(the default aligner) and samtools 1.22.1, with the scheme at commit `d7493fd` of CanSNPer2-data.
Other versions can give other output.

## 1. Check the installation

```sh
cansnper2 selftest
```

The output must list at least one usable aligner. The `read typing:` line says whether
sequencing reads can be typed (step 7 needs it). The meaning of each line is explained in
[installation.md](installation.md#check-the-installation).

## 2. Fetch the scheme

Work in a new, empty folder:

```sh
mkdir cansnper2-tutorial
cd cansnper2-tutorial
cansnper2 scheme get francisella_tularensis
```

```text
Built francisella_tularensis.db from the scheme francisella_tularensis in master of FOI-Bioinformatics/CanSNPer2-data
Downloaded the reference genomes of francisella_tularensis.db from genbank into references
```

The command fetched the scheme files from CanSNPer2-data, checked them against the checksums of
the repository, built the database `francisella_tularensis.db`, and downloaded the five reference
genomes of the scheme into `references/`. It also wrote `francisella_tularensis.source.json`,
which records the repository commit the database was built from, and a log file in `logs/`.
`cansnper2 db info -d francisella_tularensis.db` prints what the database holds: 328 nodes and 328 canSNPs on five reference
genomes. A node can carry several canSNPs, ranked within the node, of which exactly one is marked
`typed` and is called by CanSNPer2; in this scheme each node has one canSNP. These files are described in
[schemes-and-databases.md](schemes-and-databases.md).

`scheme get` reads the scheme through the GitHub API, which allows 60 anonymous requests per
hour from one address. On a shared network this limit can be reached, and the command stops or
warns with `HTTP Error 403: rate limit exceeded`. Set the environment variable `GITHUB_TOKEN` to
a GitHub token, or try again later (see
[schemes-and-databases.md](schemes-and-databases.md#github-api-limits)).

## 3. Download three genome assemblies

The three assemblies are complete genomes of strains from three *F. tularensis* subspecies:

| Strain | Subspecies | GenBank assembly |
| --- | --- | --- |
| SCHU S4 | *tularensis* | GCA_000008985.1 |
| WY96-3418 | *tularensis* | GCA_000016105.1 |
| FSC147 | *mediasiatica* | GCA_000018925.1 |

SCHU S4 is also one of the reference genomes of the scheme (`SCHUS4.2`). Download the assemblies
from the NCBI FTP site and unpack them; the file names become the sample names:

```sh
mkdir genomes
curl -sS https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/008/985/GCA_000008985.1_ASM898v1/GCA_000008985.1_ASM898v1_genomic.fna.gz | gunzip > genomes/SCHU_S4.fna
curl -sS https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/016/105/GCA_000016105.1_ASM1610v1/GCA_000016105.1_ASM1610v1_genomic.fna.gz | gunzip > genomes/WY96-3418.fna
curl -sS https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/018/925/GCA_000018925.1_ASM1892v1/GCA_000018925.1_ASM1892v1_genomic.fna.gz | gunzip > genomes/FSC147.fna
```

The NCBI `datasets` program (`datasets download genome accession ...`) can download the same
assemblies; it writes them into a zip archive under other file names, so rename them as above to
get the same sample names.

## 4. Type the assemblies

```sh
cansnper2 type genomes/*.fna -d francisella_tularensis.db
```

```text
Run alignments for 3 queries against the references in ./references/ using minimap2
FSC147: M.1
SCHU_S4: A.I.13
WY96-3418: A.II.6
1 warning about the results is in the log /path/to/cansnper2-tutorial/logs/cansnper2-type-Nov-25-2026.log, run with --verbose to see it
```

Each assembly was aligned to the five reference genomes, and the base of the sample was read at
every canSNP position. A canSNP is derived in a sample when the sample has the derived allele,
the allele that defines the branch of the tree, and ancestral when it has the allele of the
lineages outside that branch. The final SNP of each sample is printed: the deepest node of the
scheme tree that its derived canSNPs support, and the SNP path is the chain of derived canSNPs
from the root of the tree to it. These terms are defined in [method.md](method.md#terms). The warning in the
log says that one canSNP position (`B.184`) had no alignment data in WY96-3418; this canSNP is not
on the path of that sample and does not affect its call.

The results are in `results/`:

```text
results/
  summary.html
  samples/
    FSC147_snps.txt
    FSC147_not_called.txt
    SCHU_S4_snps.txt
    SCHU_S4_not_called.txt
    WY96-3418_snps.txt
    WY96-3418_not_called.txt
```

## 5. Open the overview

Open `results/summary.html` in a web browser (`open results/summary.html` on macOS,
`xdg-open results/summary.html` on Linux). The page works without a network connection. Its
header names the scheme, the database and the settings of the run, and its table lists the
samples:

| Sample | Final SNP | Support | No data | Status |
| --- | --- | --- | --- | --- |
| SCHU_S4 | A.I.13 | 9/9 | 0 | review |
| WY96-3418 | A.II.6 | 7/7 | 0 | review |
| FSC147 | M.1 | 4/4 | 0 | clean |

- **Support** is the number of derived canSNPs on the path from the root to the final SNP, out of
  the canSNPs on that path. All canSNPs on the three paths were derived.
- **No data** counts the canSNPs on the path without an aligned base. There were none.
- **Status** `clean` means that the path is confirmed, every canSNP on it has data, none lies
  near the edge of an alignment, and no derived canSNP lies outside the path. `review` means
  that the path is confirmed but something deserves a look.

Click a sample to see its path on the tree and the reason for its status. For SCHU_S4 the reason
is `1 derived SNP(s) called outside the SNP path, on another branch: B.83`, and for WY96-3418 it
is `... on another branch: B.97`. Each of these genomes carries the derived allele of one canSNP
that belongs to another branch of the tree, here the branch of subspecies *holarctica*. The call
itself is not changed by this, but such a canSNP is worth checking, for example in the
alignment, before the result is reported. The alignments are removed at the end of a run; type
the sample again with `--rerun --keep-temp` to keep them (see
[output-formats.md](output-formats.md#temporary-files)). The statuses and what to do about them are described
in [interpreting-results.md](interpreting-results.md).

## 6. Read a result file

`results/samples/WY96-3418_snps.txt` lists the canSNPs at which the sample carries the derived
allele, followed by four comment lines (the full paths in the `# Query file` line are shortened
here):

```tsv
Name	Reference	Sequence	Pos	QuerySeq	QueryPosition	Ancestral base	Derived base	Target base	Edge flag
T/N.1	SCHUS4.2	AJ749949.2	83976	CP000608.1	168877	G	A	A	
T.1	SCHUS4.2	AJ749949.2	1165688	CP000608.1	1126638	A	G	G	
A/M.1	SCHUS4.2	AJ749949.2	1491914	CP000608.1	441905	G	A	A	
A.1	SCHUS4.2	AJ749949.2	397639	CP000608.1	1594947	C	T	T	
A.II.1	SCHUS4.2	AJ749949.2	84150	CP000608.1	169051	C	T	T	
A.II.2	SCHUS4.2	AJ749949.2	201094	CP000608.1	1790280	C	T	T	
A.II.6	SCHUS4.2	AJ749949.2	29478	CP000608.1	114370	C	T	T	
B.97	FSC200	CP003862.1	1282065	CP000608.1	1352237	G	A	A	
# SNP path: T/N.1;T.1;A/M.1;A.1;A.II.1;A.II.2;A.II.6
# Final SNP: A.II.6 found/depth: 7/7
# Query file: .../genomes/WY96-3418.fna sha256:2978cb107c98a0d8f893abf8bd786dde5da7af8d3e2f8475857e289c4c646dd8
# Typing settings: database sha256:f3fc3d4822ae378ebb5922afcb3c41a311c1ac9d0102dc9395f508a02fa104a2 min_required_hits=3 strictness=0.7 mask=0 aligner=minimap2
```

Line by line:

- The header row names the columns. `Name` is the canSNP (and the node it defines), `Reference`
  the reference genome it is defined on, and `Sequence` and `Pos` the sequence and position in
  that genome. `QuerySeq` and `QueryPosition` give where the position aligned in the sample.
  `Ancestral base` and `Derived base` are the two alleles of the canSNP, and `Target base` is
  the base found in the sample. `Edge flag` is empty when the position is not near the edge of an
  alignment block.
- `T/N.1` to `A.II.6`: seven derived canSNPs that form one path from the root of the tree,
  through `A.1` and `A.II.1`, to `A.II.6`.
- `B.97`: a derived canSNP on another branch, defined on the reference genome FSC200. It is the
  reason for the `review` status in step 5.
- `# SNP path`: the derived canSNPs from the root to the final SNP.
- `# Final SNP: A.II.6 found/depth: 7/7`: the final SNP and its support, 7 derived canSNPs out
  of the 7 canSNPs on the path. (The second number is the depth of the node in the tree, which
  equals the number of canSNPs on the path in the published schemes; see
  [output-formats.md](output-formats.md#per-sample-cansnp-files).)
- `# Query file`: the input file and its SHA-256 checksum.
- `# Typing settings`: the SHA-256 of the database and the settings that affect the call. A
  later run into the same output folder reuses this result when the input, the database and
  these settings are the same, and stops with a message when they differ; `--rerun` types the
  sample again (see [output-formats.md](output-formats.md#results-from-an-earlier-run)).

`WY96-3418_not_called.txt` has the same columns for the other canSNPs of the scheme: those where
the sample carries the ancestral or another base, and those without data. All columns are
described in [output-formats.md](output-formats.md).

## 7. Type sequencing reads

Sequencing reads are typed with the same command. The run ERR5100339 in the European Nucleotide
Archive holds paired Illumina MiSeq reads of an *F. tularensis* subsp. *holarctica* isolate
(about 98 MB):

```sh
mkdir reads
curl -sS -o reads/ERR5100339_1.fastq.gz https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR510/009/ERR5100339/ERR5100339_1.fastq.gz
curl -sS -o reads/ERR5100339_2.fastq.gz https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR510/009/ERR5100339/ERR5100339_2.fastq.gz
cansnper2 type reads/ERR5100339_1.fastq.gz reads/ERR5100339_2.fastq.gz -d francisella_tularensis.db -o reads_results
```

```text
Run alignments for 1 query against the references in ./references/ using minimap2-reads
ERR5100339: B.103
```

CanSNPer2 recognized the files as reads, paired them by their names (`_1`, `_2`), mapped them to
each reference genome with minimap2 and counted the bases at the canSNP positions with samtools.
`reads_results/samples/ERR5100339_snps.txt` ends with these lines:

```text
# SNP path: T/N.1;T.1;B.1;B.2;B.3;B.5;B.12;B.72;B.13;B.26;B.42;B.168;B.66;B.33;B.103
# Final SNP: B.103 found/depth: 15/15
# Reads: platform=illumina files=2 sites_with_reads=328/328 median_depth=50 mixed_sites=0
```

The `# Reads` line says that all 328 canSNP positions were covered by reads, at a median depth of
50, and that no site was mixed. The read counts at each canSNP are in
`reads_results/samples/ERR5100339_site_counts.tsv`. How read sites are called is described in
[method.md](method.md).

## Next steps

- [interpreting-results.md](interpreting-results.md): support, statuses and what to report.
- [schemes-and-databases.md](schemes-and-databases.md): other schemes, and building your own.
- [cli.md](cli.md): all options of `cansnper2 type`.
