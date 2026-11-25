# Reproducibility

This page describes what CanSNPer2 records about a typing run and a database, how to fix the
version of a scheme, and what to report in a methods section. It is written for users who need
to trace a result back to its inputs, for example in a publication, a report or an accredited
workflow.

## What CanSNPer2 records

| Where | What |
|---|---|
| `<sample>_snps.txt`, `<sample>_not_called.txt` | The `# Query file:` line (path and sha256 of the query) and the `# Typing settings:` line (sha256 of the database file and the settings that change a call). |
| `snps_summary.txt`, `summary.html`, tree pages | CanSNPer2 version, date and time of the run, database file name and sha256, the scheme line, the reference genomes with their accessions, the aligner name and the tree walk settings. |
| the database | The `scheme_info` table: scheme name, status, source, branch or commit, the CanSNPer2 version that built it, and any later change. `cansnper2 db info` prints it. |
| `<database>.source.json` | Written by `cansnper2 scheme get` next to the database: the repository, branch and commit, and the sha256 of each scheme file. |

The columns and lines of the result files are described in
[output-formats.md](output-formats.md).

CanSNPer2 does not record the versions of the external programs: the aligner (minimap2 by default,
or progressiveMauve or MUMmer4 with `--aligner`) and, for reads, minimap2 and samtools. Only the aligner name is recorded.
Note these versions yourself, for example:

```sh
minimap2 --version
samtools --version
conda list
```

## Typing settings line

The `# Typing settings:` line of each `<sample>_snps.txt` identifies the database by the sha256 of
the file, and lists `min_required_hits`, `strictness` and `mask`. For read samples it adds the
platform and the read thresholds that were used (`min_depth`, `min_allele_fraction`, `min_mapq`,
`min_baseq`), and with `--local-references` the number of local rounds.

```text
# Query file: /data/project/genomes/GCA_000807905.1.fasta sha256:57f4df3e9755607532bcf142ebdd9c2cbb4b06ddbe296fe19473887531f76849
# Typing settings: database sha256:f3fc3d4822ae378ebb5922afcb3c41a311c1ac9d0102dc9395f508a02fa104a2 min_required_hits=3 strictness=0.7 mask=0
```

CanSNPer2 compares these two lines when it finds a result file from an earlier run in the
output folder, and uses the earlier result only when both match (see
[output-formats.md](output-formats.md#results-from-an-earlier-run)).

The sha256 of the database can be computed with `sha256sum` (or `shasum -a 256` on macOS) and
compared with this line. It identifies one database file. Two databases built from the same scheme
files can have different sha256 values, because a database also records when and from which
branch or commit it was built; to compare schemes, compare the commit and the scheme file hashes
in `db info` and `<database>.source.json`.

## Scheme information in the database

`cansnper2 db info -d <database>` prints the `scheme_info` table, then the size of the database:

```text
scheme	francisella_tularensis
status	verified
source	FOI-Bioinformatics/CanSNPer2-data
ref	master
commit	d7493fd26acfe08ff57beeed1e00ebcce0607cc3
known_issues	none
format	2.4
built_with	2.6.0
built	2026-11-25

nodes	328
canSNPs	328
typed canSNPs	328
canSNPs per node	1
reference genomes	5
```

| Key | Content |
|---|---|
| `scheme` | Name of the scheme. |
| `status` | `verified`, `draft` or `unknown`. |
| `source` | The repository (`owner/name`) or local folder the scheme files were read from. |
| `ref` | The branch, tag or commit that was asked for. |
| `commit` | The commit the files were read at; `unknown` when it could not be resolved. `<commit>+uncommitted` for a local checkout whose scheme files have uncommitted changes, and `+modified` (or `modified`) once the database was changed. |
| `known_issues` | canSNPs that the repository lists as known issues, comma-separated. |
| `format` | Layout of the database (`2.4` in CanSNPer2 2.6.0). |
| `built_with` | CanSNPer2 version that built the database, or that made the last change. |
| `built` | Date the database was built. |
| `modified` | Date of the last change made with `db modify` or `db remove`. |

`format`, `built_with` and `built` are always written. For a database built with `db create` or
`db import`, the other keys are the values given with `--scheme-name`, `--scheme-status` and
`--scheme-source`; record them when you build a scheme of your own. A database without a
`scheme_info` table (built by hand) is used with a warning that it records no scheme information.

When a database is changed with `db modify` or `db remove`, its `status` becomes `unknown`, its
`commit` is marked `+modified` and the date is recorded in `modified`, because it no longer
matches the scheme files it was built from. Each typing run against such a database logs a
warning.

## The source.json file

`cansnper2 scheme get` writes `<database>.source.json` next to the database (for
`francisella_tularensis.db`, `francisella_tularensis.source.json`):

```text
{
  "built": "2026-11-25",
  "cansnper2": "2.6.0",
  "commit": "d7493fd26acfe08ff57beeed1e00ebcce0607cc3",
  "ref": "master",
  "repository": "https://github.com/FOI-Bioinformatics/CanSNPer2-data",
  "scheme": "francisella_tularensis",
  "source_type": "CanSNPer",
  "sources": {
    "schemes/francisella_tularensis/francisella_tularensis_canSNPs.txt": "d558d6da1a1f6da1ea8d478c8c4253404cb2a71718a7a41e665f96d47802c747",
    "schemes/francisella_tularensis/francisella_tularensis_genomes.txt": "59a47c1d264a3cef92dfddd773ff49ad6077c29b65a6c319f807b4fcea1468f3",
    "schemes/francisella_tularensis/francisella_tularensis_references.txt": "19e13f84fb7df0db65ed7b518881d27708bcb4665f242fba11472e75d850c7c2",
    "schemes/francisella_tularensis/francisella_tularensis_tree.txt": "5b59747e468d5fcef452b5f5a3dcdc9d10a4a6ab5ae67bec8ec4f19d4fe72887"
  },
  "status": "verified"
}
```

`sources` gives the sha256 of each scheme file as it was read. For a scheme read with
`--data-dir`, `repository` is the folder and `commit` the commit of its git checkout, if any;
scheme files with uncommitted changes are listed under `uncommitted_sources`.

When the commit cannot be resolved (for example when the GitHub API refuses anonymous requests
for a while), the files are read from the branch, the commit is recorded as `unknown` in the
database (and as `null` in `<database>.source.json`), and a warning is printed. Run the command again later to obtain a database with a known commit.

## Fixing the scheme version

A scheme on the `master` branch of CanSNPer2-data can change. To build the same database again,
fetch the scheme at the commit that the first database records:

```sh
cansnper2 db info -d francisella_tularensis.db
cansnper2 scheme get francisella_tularensis --branch d7493fd26acfe08ff57beeed1e00ebcce0607cc3 \
    -d francisella_tularensis_d7493fd.db
```

`--branch` takes a branch, tag or commit (a shortened commit is accepted as well). The scheme
files are checked against the manifest of the repository at that commit, and the commit is
recorded in the database.

Keep the reference genome files with the database. `scheme get` and `ref get` download them by
the versioned assembly accessions of the references file, and each file is checked against the
canSNP positions and bases of the database before typing.

## What to record for a methods section

For each typing run, record:

- the CanSNPer2 version (`cansnper2 --version`);
- the aligner and its version, or for reads the minimap2 and samtools versions;
- the scheme name, its status, the repository and the commit (`cansnper2 db info`);
- the sha256 of the database file (the `# Typing settings:` line);
- the settings that differ from the defaults: `--min-required-hits`, `--strictness`, `--mask`,
  and for reads `--platform` and the read thresholds;
- for reads, whether `--local-references` was used.

A template, with the placeholders in angle brackets:

> Genomes were typed with CanSNPer2 `<version>` (minimap2 `<version>`) against the `<scheme>` scheme
> (CanSNPer2-data commit `<commit>`, status `<status>`; database sha256 `<sha256>`) with the
> default settings (`--min-required-hits 3`, `--strictness 0.7`).

Give the exact CanSNPer2 version that `cansnper2 --version` prints. When another aligner was
chosen with `--aligner`, name it and its version in place of minimap2.

For read data:

> Reads were typed with CanSNPer2 `<version>` against the `<scheme>` scheme (CanSNPer2-data commit
> `<commit>`; database sha256 `<sha256>`). Reads were mapped with minimap2 `<version>` and counted
> with samtools `<version>`, with the defaults of the `<platform>` platform
> (`<settings from the Typing settings line>`).

Cite CanSNPer2 as given in `CITATION.cff` of the repository, and the publications of the scheme
as given in its README in CanSNPer2-data.
