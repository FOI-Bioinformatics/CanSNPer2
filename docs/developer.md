# Developer notes

This page is for contributors to CanSNPer2. It describes the layout of the repository, how to run
the tests, what continuous integration checks, how to add an aligner and how a release is made.
Bug reports and pull requests are welcome on the
[GitHub repository](https://github.com/FOI-Bioinformatics/CanSNPer2).

## Repository layout

| Path | Content |
| --- | --- |
| `CanSNPer2/__init__.py` | the version (`__version__`) and the command summary shown in help texts |
| `CanSNPer2/cli.py` | the `cansnper2` command: one argument parser per verb, dispatch to the run functions below |
| `CanSNPer2/CanSNPerTree.py` | the run function of `cansnper2 type`: checks of the inputs and the database, then the typing run |
| `CanSNPer2/SNPDatabase.py` | the run function of the `db` verbs: create, import, export, info, convert, modify, remove |
| `CanSNPer2/DownloadGenomes.py` | the run function of `scheme list`, `scheme get` and `ref get`; downloads of reference and representative genomes from NCBI |
| `CanSNPer2/selftest.py` | `cansnper2 selftest` |
| `CanSNPer2/modules/CanSNPer2.py` | the typing run: alignment or read mapping, site states, the call, result files |
| `CanSNPer2/modules/Aligners.py` | the aligners (progressiveMauve, minimap2, MUMmer4) and the read mapper (minimap2 and samtools) |
| `CanSNPer2/modules/AlignmentReader.py` | what all alignment readers share: classification of the query base, edge flag, `--mask`, result rows |
| `CanSNPer2/modules/ParseXMFA.py`, `ParsePAF.py`, `ParseDelta.py` | readers of progressiveMauve (xmfa), minimap2 (PAF) and MUMmer4 (delta) alignments |
| `CanSNPer2/modules/ParsePileup.py` | reader of `samtools mpileup` output for read typing; platform defaults |
| `CanSNPer2/modules/ReadInput.py` | recognition of read files, pairing of mates, the samplesheet |
| `CanSNPer2/modules/LocalReferences.py`, `LiftOver.py` | read typing with local references: choice of the genome and lift of the canSNP positions |
| `CanSNPer2/modules/NewickTree.py` | the tree walk that decides the final SNP, the Newick tree and the HTML pages |
| `CanSNPer2/modules/DatabaseConnection.py` | SQLite access, database checks, source file readers |
| `CanSNPer2/modules/InitiateDatabase.py` | creation of a database from the scheme files |
| `CanSNPer2/modules/ImportDatabase.py`, `ExportDatabase.py` | `db import` and `db export` |
| `CanSNPer2/modules/ModifyDatabase.py`, `ModifyTree.py` | `db modify` and `db remove` |
| `CanSNPer2/modules/ReadTaxonomyCanSNPer.py` | reader of the tree file |
| `CanSNPer2/modules/SourceCheck.py` | checks of the scheme files before a database is built or changed |
| `CanSNPer2/modules/ConvertPositions.py` | `db convert` |
| `CanSNPer2/modules/FetchScheme.py` | reading schemes from GitHub or a local directory, manifest checks, `source.json` |
| `CanSNPer2/modules/SchemeInfo.py`, `SchemeGenomes.py` | the `scheme_info` table and the genomes table of a database |
| `CanSNPer2/modules/ReferenceCheck.py` | checks of reference genome files against the canSNPs of a database |
| `CanSNPer2/modules/runMauve.py` | a stand-alone progressiveMauve runner (`RunMauve`) that aligns queries against a folder of references; not used by `cansnper2 type` |
| `CanSNPer2/modules/LogSetup.py`, `TopLevel.py` | log files and console output; the handler of unexpected errors and its exit status |
| `CanSNPer2/modules/FileMode.py` | permission bits of files written through temporary files |
| `tests/` | the pytest suite; `tests/helpers.py` builds small databases and scheme files |
| `tools/gen_cli_docs.py` | generates `docs/cli.md` |
| `docs/` | the documentation |
| `.github/workflows/ci.yml` | continuous integration |
| `pyproject.toml` | package metadata, the `cansnper2` command and the build settings |

## Run the tests

Install CanSNPer2 in editable mode, together with pytest (the `test` extra), and run the suite
from the repository root:

```sh
python -m pip install -e ".[test]"
python -m pytest tests
```

Tests that need an external program or data are skipped when it is missing:

| Requirement | Tests |
| --- | --- |
| minimap2 | lift-over, PAF reading with real alignments, choice of local references |
| minimap2 and samtools | read typing, local references, the read typing command line |
| `nucmer` and `delta-filter` (MUMmer4) | delta reading with real alignments |
| progressiveMauve | xmfa reading with real alignments |
| `node` (Node.js) | the scripts of the HTML pages |
| `git` | provenance of schemes read from a local git checkout |
| a CanSNPer2-data checkout | building every published scheme; found next to the repository (`../CanSNPer2-data`) or given with the environment variable `CANSNPER2_DATA` |
| Python 3.13 | the comparison of `docs/cli.md` with the generated text |

A few tests of file permissions are skipped when they run as root. `pytest -rs` lists the
skipped tests with their reasons.

`tests/test_docs.py` checks the documentation: every `cansnper2` command line in a `sh` code
block of `README.md` and `docs/*.md` must parse with the current parsers (a line that shows an
error on purpose ends with the comment `# fails`), relative links and anchors must resolve, the
commands and option spellings before 2.6.0 may appear only in `docs/upgrading.md` and
`CHANGELOG.md`, and the spelling is US English.

## Continuous integration

`.github/workflows/ci.yml` runs on every push to `master` and on every pull request, on Ubuntu
24.04 with Python 3.11, 3.12, 3.13 and 3.14. It:

1. checks out CanSNPer2 and CanSNPer2-data;
2. installs progressiveMauve, minimap2 and samtools with `apt-get`, and MUMmer4 from Bioconda;
3. installs CanSNPer2 with pip and runs `pip check`;
4. runs `cansnper2 selftest` and the help of the main verbs;
5. runs the unit tests, with `CANSNPER2_DATA` pointing to the CanSNPer2-data checkout;
6. builds the `bacillus_anthracis` database from the scheme files, exports it, rebuilds it from
   the export and compares the two exports;
7. builds the same scheme with `cansnper2 scheme get` and compares its export;
8. downloads the reference genome with `cansnper2 ref get`;
9. types that reference genome against the database with the default aligner, with
   progressiveMauve and with MUMmer4, and types reads simulated from it; the calls must agree.

## Regenerate the command line reference

`docs/cli.md` is generated from the parsers in `CanSNPer2/cli.py`. After a change to an option
or a help text, regenerate it and commit the result:

```sh
python tools/gen_cli_docs.py
```

The script must be run with Python 3.13, because argparse formats help texts slightly
differently in other versions; it refuses to run with another version. `tests/test_cli_reference.py`
checks that every option is in the file and, under Python 3.13, that the file is current.

## Add an aligner

Each aligner is a subclass of `Aligner` in `CanSNPer2/modules/Aligners.py`. A subclass sets
`name` (the value of `--aligner`), `program` (the executable), `suffix` (the ending of the
alignment file) and `log_tag`, and implements:

- `command(reference_file, query_file, output)`: the command line as a list of arguments;
- `reader(**kwargs)`: the reader of its alignment files.

It can also override `prepare()` (files the reader needs besides the alignment), `exit_status()`,
`retry_message()` and `prepare_retry()` (how a failed run is judged and retried), and
`available()` (when more than one program is needed, as for MUMmer4).

A reader is a subclass of `AlignmentReader` (`CanSNPer2/modules/AlignmentReader.py`) that
implements `read_alignment(path)`: it walks the alignment file and reports, through `observe()`,
the query base on the forward strand of the reference at each canSNP position that an aligned
block covers. Positions are those of the database: the reference sequence (the first word of its
FASTA header) and the 1-based position on it. The classification into derived, ancestral, other
base and no data, `--mask`, the edge flag and the result rows are shared by all readers.

The new class is added to `ALIGNERS` at the end of `Aligners.py`. It can then be selected with
`--aligner`, and `cansnper2 selftest` reports it. Add a reader test with a small alignment file,
and a test that runs the program, skipped when the program is not installed.

## Make a release

1. Set the version in `CanSNPer2/__init__.py` (`__version__`); `pyproject.toml` reads it from there.
2. Move the entries under `[Unreleased]` in `CHANGELOG.md` to a section for the new version with
   its date, and update `version` and `date-released` in `CITATION.cff`.
3. Regenerate `docs/cli.md` with Python 3.13 and run the test suite.
4. Merge into `master` and check that continuous integration passes.
5. Tag the release as `vX.Y.Z` (for example `v2.6.0`) and push the tag.
6. Update the Bioconda recipe in
   [bioconda-recipes](https://github.com/bioconda/bioconda-recipes) (`recipes/cansnper2/meta.yaml`):
   the version, the SHA-256 of the GitHub archive of the tag, and, when they change, the entry
   point (`cansnper2=CanSNPer2.cli:main` since 2.6.0) and the test commands.
