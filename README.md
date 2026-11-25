# CanSNPer2

CanSNPer2 assigns bacterial genomes to lineages of a canonical SNP (canSNP) typing scheme. A
scheme is a tree whose branches are defined by canSNPs. CanSNPer2 reads the state of each canSNP
in a genome and reports the deepest node of the tree that the derived canSNPs support. It types
assemblies, which are aligned to the reference genomes of the scheme, and sequencing reads
(Illumina, Oxford Nanopore, PacBio HiFi), which are mapped to them. It also builds, exports and
modifies scheme databases, and fetches published schemes and their reference genomes.

CanSNPer2 is intended for isolates of the organisms a scheme covers. It does not identify species,
and it does not type metagenomes or resolve mixtures of strains; see
[docs/limitations.md](docs/limitations.md).

## Published schemes

The schemes are maintained in
[CanSNPer2-data](https://github.com/FOI-Bioinformatics/CanSNPer2-data), where each scheme has a
description and the publications it is based on.

| Scheme | Organism | Status | Branch |
| --- | --- | --- | --- |
| `bacillus_anthracis` | *Bacillus anthracis* | verified | `master` |
| `brucella_classical` | classical *Brucella* species | verified | `master` |
| `coxiella_burnetii` | *Coxiella burnetii* | verified | `master` |
| `francisella_tularensis` | *Francisella tularensis* | verified | `master` |
| `yersinia_pestis` | *Yersinia pestis* | verified | `master` |
| `brucella_genus`, `brucella_2017`, `brucella_anthropi`, `brucella_suis` | *Brucella* | draft | `dev` |

`cansnper2 scheme list` prints the current list. Draft schemes are under evaluation and are
fetched with `--branch dev`; see [docs/schemes-and-databases.md](docs/schemes-and-databases.md).

## Installation

With [Bioconda](https://bioconda.github.io/), on Linux or macOS:

```sh
mamba create -n cansnper2 -c conda-forge -c bioconda 'cansnper2>=2.6.0' minimap2 samtools
mamba activate cansnper2
cansnper2 selftest
```

minimap2 aligns assemblies (the default aligner) and, with samtools, maps sequencing reads.
progressiveMauve and MUMmer4 can be installed as alternative aligners. `cansnper2 selftest`
reports which programs were found. Installation from source is described in
[docs/installation.md](docs/installation.md).

## Quick start

Fetch a scheme (this builds the database `francisella_tularensis.db` and downloads its reference
genomes into `references/`), then type assemblies (uncompressed FASTA) against it:

```sh
cansnper2 scheme get francisella_tularensis
cansnper2 type genome1.fna genome2.fna -d francisella_tularensis.db -o results
```

Paired-end Illumina reads are typed the same way:

```sh
cansnper2 type sample_R1.fastq.gz sample_R2.fastq.gz -d francisella_tularensis.db -o results
```

The final SNP of each sample is printed and written to `results/samples/`. The page
`results/summary.html` gives an overview of all samples with the support and status of each
call. [docs/tutorial.md](docs/tutorial.md) goes through a complete example, and
[docs/interpreting-results.md](docs/interpreting-results.md) explains how to read a result.

## Documentation

| Page | Contents |
| --- | --- |
| [Tutorial](docs/tutorial.md) | A worked example on public genomes |
| [Installation](docs/installation.md) | Bioconda, source installation, external programs |
| [Method](docs/method.md) | How site states are read and how the final SNP is chosen |
| [Interpreting results](docs/interpreting-results.md) | Support, call status, flags, what to report |
| [Limitations](docs/limitations.md) | Intended use and known failure modes |
| [Input formats](docs/input-formats.md) | Query files, samplesheet, scheme files |
| [Output formats](docs/output-formats.md) | Result files and their columns |
| [Schemes and databases](docs/schemes-and-databases.md) | Fetching schemes, building and changing databases |
| [Command line reference](docs/cli.md) | All commands and options |
| [Reproducibility](docs/reproducibility.md) | What is recorded, and what to report in a methods section |
| [Validation](docs/validation.md) | Evaluations of the aligners and of read typing |
| [Upgrading](docs/upgrading.md) | From the commands before 2.6.0 and from older databases |
| [Troubleshooting](docs/troubleshooting.md) | Exit codes and messages |
| [Developer notes](docs/developer.md) | Tests, continuous integration, releases |

Changes between versions are listed in [CHANGELOG.md](CHANGELOG.md).

## Citation

If you use CanSNPer2, please cite:

> Lärkeryd A, Myrtennäs K, Karlsson E, Dwibedi CK, Forsman M, Larsson P, Johansson A, Sjödin A.
> CanSNPer: a hierarchical genotype classifier of clonal pathogens. *Bioinformatics*
> 2014;30(12):1762-1764. https://doi.org/10.1093/bioinformatics/btu113

together with the CanSNPer2 version you used and the publications of the scheme, which are listed
in the description of each scheme in CanSNPer2-data (`schemes/<name>/README.md`, for example
[francisella_tularensis](https://github.com/FOI-Bioinformatics/CanSNPer2-data/blob/master/schemes/francisella_tularensis/README.md)). [CITATION.cff](CITATION.cff) holds the same
information in machine-readable form.

## License and support

CanSNPer2 is developed by the bioinformatics group of the Swedish Defence Research Agency (FOI)
and is distributed under the [GNU General Public License v3.0](LICENSE).

Questions, bug reports and suggestions are welcome in the
[issue tracker](https://github.com/FOI-Bioinformatics/CanSNPer2/issues).
[docs/troubleshooting.md](docs/troubleshooting.md) lists what to include in a report.
