# Upgrading

This page is for users of earlier CanSNPer2 versions. It maps the commands and options used before
2.6.0 to the `cansnper2` command, lists other changes that affect scripts and pipelines, and
describes how databases and SNP files of 2.3.0 and earlier are converted.

## From the commands before 2.6.0

Up to version 2.4.0, CanSNPer2 installed four commands: `CanSNPer2`, `CanSNPer2-database`,
`CanSNPer2-download` and `CanSNPer2-test`. Version 2.6.0 installs one command, `cansnper2`, with
a verb for every task. The earlier commands are no longer installed; if copies of them remain in
`PATH` from an earlier installation, remove them.

### Commands

| Before 2.6.0 | From 2.6.0 |
| --- | --- |
| `CanSNPer2 QUERY... -db DB` | `cansnper2 type QUERY... -d DB` |
| `CanSNPer2 --version` | `cansnper2 --version` |
| `CanSNPer2-test` | `cansnper2 selftest` |
| `CanSNPer2-database -db DB --create --tree T --annotation A --references R` | `cansnper2 db create -d DB --tree T --cansnps A --references R` |
| `CanSNPer2-database -db DB --import --tree T --annotation A --references R` | `cansnper2 db import -d DB --tree T --cansnps A --references R` |
| `CanSNPer2-database -db DB --export [--export_format F] [-o DIR]` | `cansnper2 db export -d DB [--format F] [-o DIR]` |
| `CanSNPer2-database -db DB --info` | `cansnper2 db info -d DB` |
| `CanSNPer2-database --convert -db DB --refdir DIR -o DIR` | `cansnper2 db convert -d DB --ref-dir DIR -o DIR` |
| `CanSNPer2-database --convert --annotation FILE --refdir DIR -o DIR` | `cansnper2 db convert --cansnps FILE --ref-dir DIR -o DIR` |
| `CanSNPer2-database -db DB --mod_file F --parent N [--replace]` | `cansnper2 db modify -d DB --mod-file F --parent N [--replace]` |
| `CanSNPer2-database -db DB --annotation A` or `--references R` (existing database) | `cansnper2 db modify -d DB --cansnps A` or `--references R` |
| `CanSNPer2-database -db DB --remove --parent N` | `cansnper2 db remove -d DB --parent N` |
| `CanSNPer2-download --list_schemes` | `cansnper2 scheme list` |
| `CanSNPer2-download --scheme NAME` | `cansnper2 scheme get NAME` |
| `CanSNPer2-download -db DB` | `cansnper2 ref get -d DB` |

`CanSNPer2-download --scheme NAME -s refseq -o DIR` built the database and downloaded the
reference genomes from RefSeq into `DIR` in one step. `cansnper2 scheme get` always downloads
from GenBank into `./references`; for other choices, build the database first and download the
references with `ref get`:

```sh
cansnper2 scheme get NAME --no-references
cansnper2 ref get -d NAME.db -s refseq -o DIR
```

### Options

As a rule, underscores in option names became hyphens (`--min_required_hits` is now
`--min-required-hits`). The table lists every option whose name changed in another way, or that
was removed. The options of the read typing and local reference functions, which were added
during the development of 2.6.0, follow the same rule.

| Command before 2.6.0 | Option before 2.6.0 | From 2.6.0 |
| --- | --- | --- |
| `CanSNPer2`, `CanSNPer2-database`, `CanSNPer2-download` | `-db` | `-d` (the long form `--database` is unchanged) |
| `CanSNPer2`, `CanSNPer2-database`, `CanSNPer2-download` | `--verbose` | `-v` or `--verbose` |
| `CanSNPer2`, `CanSNPer2-database` | `--supress` | `--quiet` |
| `CanSNPer2` | `--logdir` | `--log-dir` |
| `CanSNPer2-database`, `CanSNPer2-download` | `--logs` | `--log-dir` |
| `CanSNPer2` | `--refdir` | `--ref-dir` |
| `CanSNPer2-database` | `--refdir` | `--ref-dir` (`db convert`) |
| `CanSNPer2` | `--tmpdir` | `--tmp-dir` |
| `CanSNPer2` | `--annotations` | `--node-labels` |
| `CanSNPer2` | `--no_snpfiles` | `--no-snp-files` |
| `CanSNPer2` | `--skip_alignment`, `--skip_mauve` | `--skip-alignment` |
| `CanSNPer2` | `--read_input` | removed; sequencing reads are recognized by the content of the files |
| `CanSNPer2-database` | `--annotation` | `--cansnps` |
| `CanSNPer2-database` | `--export_format` | `--format` (`db export`) |
| `CanSNPer2-database` | `--create`, `--import`, `--export`, `--info`, `--convert`, `--remove` | the verbs `db create`, `db import`, `db export`, `db info`, `db convert`, `db remove` |
| `CanSNPer2-database` | `--source_type` | removed; the source files are always read in the CanSNPer format |
| `CanSNPer2-database` | `--tmpdir` | removed |
| `CanSNPer2-download` | `--list_schemes` | the verb `scheme list` |
| `CanSNPer2-download` | `--scheme NAME` | the verb `scheme get NAME` |
| `CanSNPer2-download` | `--no_references` | `--no-references` (`scheme get` only) |

Options that keep their spelling: `-o`/`--outdir`, `--summary`, `--workdir`, `--strictness`,
`--mask`, `--rerun`, `--aligner`, `--debug`, `--tree`, `--references`, `--parent`, `--replace`,
`--repo`, `--branch`, `--overwrite`, `-s`/`--source`.

### Other changes in 2.6.0 that affect scripts

- The default aligner for assemblies is minimap2. Up to 2.4.0 it was progressiveMauve; give
  `--aligner progressiveMauve` to keep the earlier aligner.
- The per-sample result files (`<sample>_snps.txt`, `<sample>_not_called.txt` and, with
  `--save-tree`, `<sample>_tree.html`) are written to `<outdir>/samples/`. Up to 2.4.0 they were
  written to `<outdir>` itself. `summary.html` and the files of `--summary` stay in `<outdir>`.
  The result files are described in [output-formats.md](output-formats.md).
- The log files are named after the verb, for example `logs/cansnper2-type-Nov-25-2026.log` or
  `logs/cansnper2-db-info-Nov-25-2026.log`. Up to 2.4.0 they were named
  `CanSNPer2-<date>.log`, `CanSNPer2-database-<date>.log` and `CanSNPer2-download-<date>.log`. A
  second log of the same verb on the same day has the time added as `-HHMM` (up to 2.4.0,
  `-HH:MM`).
- Databases built by 2.4.0 are read by 2.6.0 without conversion; the database format is still
  2.4.

## From 2.3 to 2.4: canSNP files

Version 2.4.0 changed where the canSNPs of a scheme come from. Up to 2.3.0 a scheme had a SNP
file, `<name>_snp.txt`, with one canSNP per node. From 2.4.0 it has a canSNP file,
`<name>_canSNPs.txt`, which holds every canSNP of every node, ranked within the node, with the
column `typed` marking the one canSNP per node that CanSNPer2 calls; the others are alternatives.
The database stores all of them, so `db export` writes them back.

A database of 2.3.0 is refused:

```text
The database old.db cannot be used. Its snp_annotation table holds one canSNP per node, without the rank, score and typed columns of the current format. Convert it with cansnper2 db convert, or fetch the current database with cansnper2 scheme get.
```

`db convert` upgrades the database, or a `<name>_snp.txt`. Neither needs reference files,
because no position changes:

```sh
cansnper2 db convert -d my_scheme.db -o upgraded/
cansnper2 db convert --cansnps my_scheme_snp.txt -o upgraded/
```

Each canSNP becomes rank 1 of its node, typed, without a score. A converted SNP file is written
as `upgraded/my_scheme_canSNPs.txt`. A `<name>_snp.txt` of 2.3.0 is also read as it is by
`db create --cansnps`, with the same result, so a scheme in that format can still be built.
Published schemes can instead be fetched again with `cansnper2 scheme get`.

## From 2.2 to 2.3: positions on each reference sequence

Version 2.3.0 changed how a canSNP position is stored. Up to 2.2 it was one number counted along
all sequences of the reference genome joined in the order of the reference FASTA file, which
made a result depend on that order. From 2.3.0 the database and the source files hold the
sequence and the position on it (the columns `sequence` and `sequence_position`; the column
`position` was removed). The result tables gained the columns `Sequence`, `QuerySeq` and
`QueryPosition`, and `Pos` is the position on `Sequence`.

A database or SNP file in the format before 2.3.0 is refused:

```text
The database old.db cannot be used. It stores SNP positions along the joined reference sequences (the format before 2.3.0). Convert it with cansnper2 db convert, or fetch the current scheme with cansnper2 scheme get.
```

`db convert` converts it to the current format in one step. The positions are located on the
sequences of the reference genomes, which are read from `--ref-dir` as `<genome>.fna`:

```sh
cansnper2 db convert -d old.db --ref-dir references/ -o converted/
cansnper2 db convert --cansnps old_snp.txt --ref-dir references/ -o converted/
```

- The sequences of each reference file must be in the order in which the old positions were
  counted, which is the order of the NCBI assembly file. The sequence names that are written
  are the first words of the FASTA headers.
- For a database, reference files missing from `--ref-dir` are downloaded from GenBank with the
  accessions the database records, so the sequence names are GenBank accessions. A SNP file
  holds no accessions, so all its reference files must be in `--ref-dir` beforehand.
  (`CanSNPer2-database --convert` could take a references file with `--references` to download
  them; `db convert` has no such option.)
- Each canSNP is checked against its reference genome. When the base at a position matches
  neither the derived nor the ancestral base, the mismatching canSNPs are listed and nothing is
  written.
- A database made by 2.3.0 or later cannot be read by 2.2.

Alignments kept from runs before 2.3.0 (`--keep-temp`) are affected as follows when they are read again
with `--skip-alignment`: progressiveMauve alignments must be made again, because the sequence
table that 2.3.0 and later write next to an xmfa alignment (`<alignment>.sequences`) is missing;
minimap2 and MUMmer4 alignments are read as they are, because their reference coordinates are
positions on each sequence.

The `db convert` command is also described in
[schemes-and-databases.md](schemes-and-databases.md#convert-databases-from-earlier-versions).
