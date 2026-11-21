# SNP positions per reference sequence

Date: 2026-11-21. Repositories: CanSNPer2, CanSNPer2-data, CanSNPCreator.

## Intent

A canSNP position is stored today as one integer along the sequences of the reference
genome (chromosomes, plasmids) joined end to end in the order of the reference fasta file.
Every result therefore depends on the order of the sequences in that file, and nothing in
the database records which sequences there are. The 2026-11-20 audit showed the effect: a
Brucella reference with its two chromosomes in another order passed the length check, and
all three aligners silently lost most calls. The base check of CanSNPer2 PR #67 now refuses
such a file, but the dependence on file order remains.

Each SNP is to be stored as a reference sequence (named by its GenBank accession, for
example AE008917.1) and a 1-based position on that sequence. The order of the sequences in
a reference file then no longer matters, positions can be compared directly with
publications that give them per chromosome, and the minimap2 and MUMmer4 readers no longer
need offsets.

## Decisions taken with the user

- Convert everything; drop joined positions (option 2). The published schemes, the
  databases and the export format change once. A database or source file with joined
  positions is refused with a message that says how to convert or rebuild it. There is no
  dual-format support.
- Sequences are named by their GenBank accession (approach A). A RefSeq download is renamed
  to the GenBank names when it is fetched. A reference file that lacks a sequence the
  database names stops the run with an error that lists the expected names.
- The result must be checked with all three aligners: progressiveMauve, minimap2, MUMmer4.
- Old CanSNPer2 versions must fail loudly on the new data, not misread it. The position
  column is therefore renamed to `sequence_position` in the source files and in the
  database. CanSNPer2 2.2.x typing with a new database stops on a database error (no
  column `position`), and `--create` on master refuses a new source file through the input
  checks of #66 (missing column `position`). Kept as `position`, an old version would read
  the chromosome II positions as joined positions and give wrong or NA calls without an
  error (2 SNPs in brucella_classical, 2 in brucella_suis, 7 in brucella_anthropi).

## Goals

1. Source files, database and results carry (sequence, position) for every SNP.
2. Typing gives the same calls and the same base at every SNP as now, with all three
   aligners, on a correct reference.
3. Typing gives the same calls whatever the order of the sequences in the reference file.
4. Every published scheme in CanSNPer2-data (master and dev) is converted and rebuilt, and
   CanSNPCreator exports the new format.

## Non-goals

- The query side: query sequence and query position are reported per contig already.
- Reading databases or source files with joined positions.
- Reference files named with RefSeq accessions supplied by the user (they must be renamed;
  the error says so).
- New SNPs on plasmids; none of the current schemes has one.

## Anti-goals

- Any change of call or base on a correct reference: the conversion must be exact.
- A silent result on a reference that lacks a named sequence: it must stop with an error.

## Approaches considered

- **A. Sequence column with GenBank names; RefSeq downloads renamed (chosen).** No new data
  file; every reader looks up (sequence, position). Cost: user-supplied RefSeq-named
  references must be renamed.
- **B. A plus a per-reference sequence table** (`<name>_sequences.txt`: GenBank name, RefSeq
  name, length) stored in the database. Accepts either naming and checks lengths without
  bases, but adds one file per scheme to keep in step with NCBI; the base check already
  catches a wrong genome. Not chosen.
- **C. Sequence index (first, second, ...) instead of a name.** Keeps the order dependence.
  Rejected.

## Design

### Source format (CanSNPer2 `--create`, CanSNPer2-data)

In `<name>_snp.txt` the `position` column is replaced by `sequence` and
`sequence_position`:

```
snp_id	strain	reference	genome	sequence	sequence_position	derived_base	ancestral_base
```

`<name>_canSNPs.txt` (CanSNPer2-data and CanSNPCreator only) likewise:

```
node	rank	genome	sequence	sequence_position	derived_base	ancestral_base	score
```

`sequence` is the GenBank accession with version of a sequence of the assembly named in
the references file; `sequence_position` is 1-based on that sequence. No two SNPs may
share (genome, sequence, sequence_position). A file with the old header is refused by
`--create` (SourceCheck) and by `validate.py` with a message naming the conversion command.

### Database

- In `snp_annotation` the column `position` is replaced by `sequence` (text) and
  `sequence_position` (integer). The new name is what makes old versions fail (see
  Decisions).
- A database whose `snp_annotation` has no `sequence_position` column is refused when it
  is opened for typing, download or export, with a message: built with joined positions,
  convert it with the conversion command, rebuild it from converted sources, or download
  the current database.
- `ExportDatabase` (source export) writes the new header; `ImportDatabase` reads it.

### Conversion command

A CanSNPer2 command converts a source SNP file with joined positions, given the reference
fasta of each genome in it, into the new format (sequence from the header's first word,
position on that sequence). It refuses a SNP whose converted base is neither the derived
nor the ancestral base, so a conversion against a reference in the wrong order cannot
pass. It also accepts a database built with joined positions and writes a converted
database, so a user who only has an old `.db` file is not stuck (my call). It is used to
convert CanSNPer2-data and is available to users with their own schemes. The exact command name and options are left to the builder; it should sit with
the other `CanSNPer2-database` operations if that fits the CLI.

### Typing

- `AlignmentReader.snplist` is keyed by (sequence, position); `snps_in` takes a sequence
  and a range. `observe` takes the sequence.
- minimap2 (`ParsePAF`): target name and position as given; the `.sequences` offset table
  is no longer needed for lookup.
- MUMmer4 (`ParseDelta`): reference name and position as given; offsets removed.
- progressiveMauve (`ParseXMFA`): the xmfa gives joined coordinates over the reference file
  that was aligned. The reader converts them to (sequence, position) with the sequence
  table of that file. The table is written next to the alignment at alignment time (the
  `prepare()` hook, as minimap2 does now), so `--skip_alignment` still works without the
  reference folder. Reverse-strand reference blocks must convert correctly.
- Reference check (run start and download): every sequence the database names for that
  reference must be present in the file (first word of the header), be long enough for
  its SNPs, and carry the ancestral or derived base at each SNP (the PR #67 rule, now per
  sequence). Extra sequences in the file are allowed; their order does not matter.
- `CanSNPer2-download --source refseq` fetches the NCBI assembly report with the genome and
  renames each sequence to its GenBank accession. `--source genbank` (default) is unchanged.

### Results

The per-sample tables become:

```
Name	Reference	Sequence	Pos	QuerySeq	QueryPosition	Ancestral base	Derived base	Target base	Edge flag
```

`Pos` is the position on `Sequence`. The `ReferenceSeq` / `ReferenceSeqPosition` columns of
PR #67, never released, are dropped in favour of this. The no-data rows (NewickTree) and the
dashboard follow the same columns.

### CanSNPer2-data

- Convert `_snp.txt` and `_canSNPs.txt` of every scheme on master and on dev with the
  conversion command; rebuild the databases and manifests with `build.py --write`.
- `validate.py`: new headers, uniqueness on (genome, sequence, sequence_position).
- `check_bases.py`: check each SNP on its named sequence; a missing sequence is a problem.
- `schemes/README.md` and the scheme READMEs that describe joined positions (for example
  the Origin section of `brucella_classical`) are updated.
- CI installs CanSNPer2 from master; see Order of work for how the conversion is tested
  against `dev_2.3.0`.

### CanSNPCreator

- `cansnp_export` writes `sequence` and `sequence_position` (it already has contig and
  position per SNP). The contig names must be the GenBank accessions of the reference.
  Schemes made from RefSeq assemblies (NC_/NZ_ names, as in the RepGenR2/GTDB workflow)
  are converted at export through the existing `--contig-aliases` option, with the alias
  table taken from the NCBI assembly report (my call).
- `cansnp_audit` and any other reader of scheme files read the new format.
- A release with a version bump.

## Verification

1. Unit tests in each repository for every changed reader, check and format.
2. The 2026-11-20 harness (52 Brucella genomes from 13 species, closed and draft, each also
   with its sequences reordered and every other one reverse-complemented), run with the new
   code and all three aligners:
   - per SNP, the same base and the same call as the old code on the same inputs;
   - with the reference sequences in another order, the same result as with the original
     reference (old code: most calls lost);
   - a reference lacking a sequence stops with an error naming it.
3. A regression on single-sequence schemes: a sample of F. tularensis, Y. pestis and
   B. anthracis genomes with all three aligners gives the same calls as before.
4. CanSNPer2-data CI (validate, check_bases, build) green on master and dev.

## Order of work

All CanSNPer2 work for this change is merged into the branch `dev_2.3.0` on
FOI-Bioinformatics/CanSNPer2 (decided by the user; created from master f6900e7), not into
master.

1. CanSNPer2 PR #67 (base check) is merged into `dev_2.3.0` (retargeted from master), and
   CanSNPer2-data PR #32 into CanSNPer2-data master (it needs no CanSNPer2 change).
2. CanSNPer2: the change on top of #67 (branch `per_sequence_positions`), PR into
   `dev_2.3.0`.
3. CanSNPer2-data: convert and rebuild on a branch of its own. Its CI installs CanSNPer2
   from `dev_2.3.0` instead of master while it is open (my call), so it is tested against
   the code that reads the new format.
4. Release: `dev_2.3.0` into CanSNPer2 master as 2.3.0, then the CanSNPer2-data
   conversion into its master with the CI install set back to master, back to back, so
   `--scheme` downloads with the released version are broken for minutes, not days.
   Then CanSNPer2-data master into dev.
5. CanSNPCreator: export and audit, release.

## Left to the builder

- Name and CLI placement of the conversion command.
- Internal data structures of the readers.
- Exact wording of error messages (modest, naming the file and the fix).

## Calls made during the design (not decided by the user)

- The conversion command also converts an existing `.db` file.
- CanSNPCreator converts RefSeq contig names at export via `--contig-aliases` and the NCBI
  assembly report.
- CanSNPer2-data CI installs CanSNPer2 from `dev_2.3.0` while the conversion is open;
  the conversion and the 2.3.0 release are merged back to back.
- Result tables keep the column name `Pos` (now the position on `Sequence`); they are not
  read back by CanSNPer2, so the rename is not needed there.
