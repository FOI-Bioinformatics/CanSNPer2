# Changelog

This file lists the changes of each CanSNPer2 release that affect users. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/). Versions follow the version number in
`CanSNPer2/__init__.py`; releases with a Git tag are marked. Older changes are recorded in the
Git history.

## [Unreleased]

### Changed

- minimap2 is the default aligner for assemblies, instead of progressiveMauve. minimap2 is
  available for Linux and macOS and is also used for read typing. `--aligner progressiveMauve`
  gives the earlier behavior.
- Documentation reorganized into a short README and topic pages in `docs/`, including a command
  line reference generated from the program (`docs/cli.md`), a description of the typing method,
  and pages on interpreting results, limitations and reproducibility.
- Help texts and messages name options as `cansnper2` spells them (for example `--allow-unverified`
  instead of `--allow_unverified`). The help shows a default only for options that take a value.
- The conda recipe was removed from the repository; it is maintained in bioconda-recipes.

## [2.6.0] - 2026-11-25

### Added

- One command, `cansnper2`, with a verb for every task: `type`, `selftest`,
  `db create|import|export|info|convert|modify|remove`, `scheme list|get` and `ref get`.
- Typing of sequencing reads (Illumina, Oxford Nanopore, PacBio HiFi): reads are mapped to the
  reference genomes with minimap2 and the bases at the canSNP positions are counted with samtools.
  Paired files are grouped by name; a samplesheet can be given with `--samplesheet`. Per-site read
  counts are written to `<sample>_site_counts.tsv`, and samples with mixed sites are flagged.
- Local references for read typing (`--local-references`): the reads are mapped again to the
  representative genome of the scheme closest to the sample, and the canSNPs are read at their
  positions in that genome. Representative genomes are stored in the database from the genomes
  file of a scheme (`--genomes`) and downloaded with `ref get --local-references` or during typing
  with `--fetch-local`.
- An overview page of all samples (`summary.html`).

### Changed

- The commands `CanSNPer2`, `CanSNPer2-database`, `CanSNPer2-download` and `CanSNPer2-test` are
  no longer installed. Options use hyphens (`--min-required-hits`). See
  [docs/upgrading.md](docs/upgrading.md).
- The overview page `summary.html` shows the tree of every sample; a separate tree page per
  sample is written with `--save-tree`.

## [2.4.0] - 2026-11-22 (tag v2.4.0)

### Added

- Databases store every canSNP of a node, ranked, with the typed canSNP marked; the canSNP file of
  a scheme is downloaded with it.
- Databases record the scheme name, status and source (`scheme_info`), shown by `db info` and in
  the results.
- `--allow-missing-references` to type when some reference genomes are missing.

### Changed

- Result files are a table followed by `#` comment lines that record the scheme, the run settings
  and the support of each call. Result pages were redrawn.
- A result made with another database or other settings is not reused.
- Many input checks for database building, modification and typing: positions out of range,
  empty or header-only files, repeated genome names or accessions, byte order marks, non-UTF-8
  and gzip-compressed source files, and control characters in scheme options.
- The default temporary folder falls back to the output folder when the system temporary folder
  cannot be written.

## [2.3.0] - 2026-11-21

### Added

- MUMmer4 as an optional aligner (`--aligner mummer4`).
- Databases can be built from CanSNPer2-data schemes, local folders or other repositories, and
  exported to the files that database building reads.

### Changed

- canSNP positions are stored per reference sequence instead of on concatenated genomes, so
  references with several chromosomes or plasmids are handled correctly with every aligner.
  Databases and canSNP files of earlier versions are converted with `db convert` (then
  `CanSNPer2-database --convert`).
- Reference bases are checked at the canSNP positions when a database is built.
- MUMmer4 runs `nucmer` with one thread so that alignments are reproducible.

### Fixed

- Duplicate sample and sequence names are refused instead of giving wrong results.

## [2.2.2] - 2026-11-17

### Added

- minimap2 as an alternative aligner (`--aligner minimap2`).

### Fixed

- The final call follows the ancestor chain of the node: two complete paths that compete at the
  same depth give no final SNP and a message.
- canSNPs at the first or last position of an alignment block are no longer dropped.
- Typing no longer stops responding when a worker stalls; truncated reference downloads, blank
  lines in inputs and the reported version were corrected.

## [2.2.0] - 2026-11-11 (tag v2.2.0)

### Changed

- The tree is drawn as a self-contained HTML page; the ete3 dependency was removed, and Python
  3.13 and 3.14 are supported.

[Unreleased]: https://github.com/FOI-Bioinformatics/CanSNPer2/compare/v2.6.0...HEAD
[2.6.0]: https://github.com/FOI-Bioinformatics/CanSNPer2/compare/v2.4.0...v2.6.0
[2.4.0]: https://github.com/FOI-Bioinformatics/CanSNPer2/compare/v2.2.0...v2.4.0
[2.2.0]: https://github.com/FOI-Bioinformatics/CanSNPer2/releases/tag/v2.2.0
