# SNP positions per reference sequence: implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Store, read and report every canSNP as (reference sequence, position on that sequence) instead of a position on the joined reference, in CanSNPer2 and CanSNPer2-data, with identical calls for all three aligners.

**Architecture:** The database and source files get `sequence` + `sequence_position` in place of `position`. `XMFAFunctions.get_snps` and `AlignmentReader` are keyed by (sequence, position). minimap2 and MUMmer4 report per-sequence coordinates already, so their offset tables go. progressiveMauve keeps joined coordinates in its xmfa, and its reader maps them back through a sequence table written next to the alignment. A converter (`CanSNPer2-database --convert`) rewrites old source files and databases and checks each base.

**Tech Stack:** Python 3.11-3.14, sqlite3, pytest. There are no new dependencies.

**Spec:** `CanSNPer2-perseq/docs/superpowers/specs/2026-11-21-per-sequence-positions-design.md` (commits d57dae5, 21ffa80). Read it with this plan.

**Context:** On 2026-11-20 the audit showed that a Brucella reference with its two chromosomes in another order passed the length check, and every aligner then silently lost calls. PR #67 (merged into `dev_2.3.0` as 8afeef6) refuses such files, but the positions still depend on file order. The user chose a per-sequence system, with the following settled:
- convert everything;
- refuse old data;
- name sequences by GenBank accession;
- rename the column so old versions fail loudly;
- merge all CanSNPer2 work into `dev_2.3.0`.

## Global Constraints

- New SNP header: `snp_id strain reference genome sequence sequence_position derived_base ancestral_base`.
- New canSNPs header: `node rank genome sequence sequence_position derived_base ancestral_base score`.
- Database `snp_annotation` columns: `node_id, snp_id, sequence, sequence_position, ancestral_base, derived_base, reference, date, genome_i`. There is no `position` column.
- Result columns: `Name Reference Sequence Pos QuerySeq QueryPosition Ancestral base Derived base Target base Edge flag`. `Pos` is the position on `Sequence`. The PR #67 `ReferenceSeq`/`ReferenceSeqPosition` columns are removed.
- `sequence` is the GenBank accession with its version, matched against the first word of a fasta header.
- Uniqueness is on (genome, sequence, sequence_position).
- Old-format input is refused, never read. Every refusal names the conversion command.
- `MISMATCH_FRACTION = 0.2` (from PR #67) is kept, now counted over all SNPs of a reference.
- Modest scientific language in code, docs and messages. No Unicode in Nextflow files.
- Commits and PRs carry no Claude attribution; the user is the sole author.
- CanSNPer2 work happens in worktree `/Users/andreassjodin/Code/cansnp-test/CanSNPer2-perseq`, branch `per_sequence_positions`, with PRs into `dev_2.3.0`. Version `2.3.0`.
- Test command: `cd CanSNPer2-perseq && PYTHONPATH=$PWD python -m pytest -q`. Use the Python 3.11 environment that has CanSNPer2 installed in editable mode, not miniforge3's python.

## Review Focus

1. **Edge flag index:** the dashboard reads the edge flag as `row[8]` (`CanSNPer2.py:1038`). With `Sequence` inserted at index 2, edge reasons would silently disappear. The test goes in Task 2.
2. **Mauve block across a sequence boundary:** an LCB can span the end of chromosome I and the start of chromosome II in the joined reference. SNPs on both sides must be found. The test goes in Task 3.
3. **Kept 2.2.x alignments:** `--skip_alignment` on an xmfa without a `.sequences` table must stop with an error that says to realign, not type NA. The test goes in Task 3.
4. **Reference file variations:** extra sequences (a plasmid), another order, and header lines with a description after the accession must all type the same as the plain reference. The tests go in Task 4.
5. **Old databases and the converter:** an old database given to `CanSNPer2`, `CanSNPer2-download` or `--export` is refused with the conversion hint. `--convert` never overwrites its input. The tests go in Tasks 1 and 5.

---

## Part A: CanSNPer2 (worktree CanSNPer2-perseq)

### Task 1: Stored format (database, loaders, source files)

Tasks 1 and 2 together switch the format. After Task 1 only the listed test files must pass; the full suite is green again at the end of Task 2. Review and merge them as a pair.

**Files:**
- Create: `tests/helpers.py`, a shared fixture builder that replaces the near-identical `make_db` copies in about 15 test files.
- Modify:
  - `CanSNPer2/modules/InitiateDatabase.py:30-41` (schema) and `:96-117` (`add_annotation`).
  - `CanSNPer2/modules/ModifyDatabase.py:25-48` (`add_annotation`).
  - `CanSNPer2/modules/DatabaseConnection.py`: `ANNOTATION_SOURCE_COLUMNS` (:31), `check_database` (:45), `LoadCanSNPAnnotation.load_annotations` headers (:223), `XMFAFunctions.get_results`/`get_snps` (:1183-1213).
  - `CanSNPer2/modules/SourceCheck.py:113-160` (`check_snps`).
  - `CanSNPer2/modules/ExportDatabase.py:25,85-90`.
  - `CanSNPer2/modules/NewickTree.py:1180-1212` (no-data rows).
  - `CanSNPer2/SNPDatabase.py`, `CanSNPer2/CanSNPerTree.py:140`, `CanSNPer2/DownloadGenomes.py:293`: callers of `check_database`.
- Test: `tests/test_db_build_and_inputs.py`, `tests/test_create_input_checks.py`, `tests/test_export.py`, `tests/test_import.py`, `tests/test_load_annotations.py`, `tests/test_newick_tree.py`.

**Interfaces:**
- Produces, in `tests/helpers.py`:
  - `make_db(path, rows, genome="REF", joined=False)`. Each row is `(snp_id, sequence, position, ancestral, derived)`. With `joined=True` it writes the old table (a `position` column, and the sequence item of each row ignored), for refusal tests.
  - `write_sources(folder, snp_rows, ...)`, which writes the tree, SNP and references triple in the new header.
- Produces, in `DatabaseConnection`:
  - `ANNOTATION_SOURCE_COLUMNS = ("snp_id","reference","genome","sequence","sequence_position","derived_base","ancestral_base")`.
  - `JOINED_FORMAT_HINT`, the fixed sentence naming `CanSNPer2-database --convert` and the current database download.
  - `check_database(path, current_format=True)`, which raises `ConnectionError` with `JOINED_FORMAT_HINT` when `snp_annotation` has no `sequence_position` column. `current_format=False` skips that check; only `--convert` uses it.
  - `XMFAFunctions.get_snps(reference) -> (snps, positions)`. `snps` is `{(sequence, position): (ancestral, derived, snp_id)}` and `positions` is `{sequence: sorted [position]}`. A duplicate (sequence, position) raises, as now.
- Produces, in `ExportDatabase`: `SNP_COLUMNS = ("snp_id","strain","reference","genome","sequence","sequence_position","derived_base","ancestral_base")`.
- Produces, in NewickTree: no-data rows of the form `[snp_id, genome, sequence, str(position), "N/A","N/A", anc, der, "N/A", ""]`. An unannotated node gives `[name,"-","-","-","N/A","N/A","-","-","N/A",""]`.

Steps:

- [ ] **Step 1: Write `tests/helpers.py`.** Point the six test files above at it, and change their fixtures to the new rows and headers.
- [ ] **Step 2: Write the failing tests:**
  - `test_create_stores_sequence_and_sequence_position`: build with `--create` from `write_sources`. Then `SELECT sequence, sequence_position FROM snp_annotation` equals the input, and `PRAGMA table_info` has no `position`.
  - `test_old_header_is_refused_with_the_conversion_hint`: `check_sources` on a SNP file with a `position` header gives one problem containing `"--convert"`.
  - `test_same_position_on_two_sequences_is_allowed` and `test_same_sequence_and_position_twice_is_refused`, both through `check_sources`.
  - `test_empty_sequence_is_refused` (`check_sources`).
  - `test_joined_database_is_refused`: `check_database(make_db(..., joined=True))` raises `ConnectionError`, and the message contains `"--convert"`. The same call with `current_format=False` passes.
  - `test_cli_refuses_joined_database`: for `CanSNPer2`, `CanSNPer2-download` and `CanSNPer2-database --export`, run as subprocesses on a joined database, each exits non-zero and stderr contains `--convert`.
  - `test_export_import_round_trip_keeps_sequences`: export, then import, then export again; the files are byte-identical and the header equals `SNP_COLUMNS`.
  - `test_get_snps_keys_by_sequence`: two SNPs at position 50, one on `chrI` and one on `chrII`, both come back keyed apart.
  - `test_no_data_rows_carry_the_sequence`: in the NewickTree rows, `row[2] == "chrII"` and `row[3] == "50"`.
- [ ] **Step 3: Run** `pytest` on the six files. Expected: the new tests fail on the `position` column.
- [ ] **Step 4: Implement the interfaces.**
  - `SourceCheck.check_snps`: when `position` is in the header and `sequence_position` is not, report one problem: "joined positions (the format before 2.3.0); convert with CanSNPer2-database --convert". Check `sequence` for being non-empty and `sequence_position` for being a positive integer. Count uniqueness on `(genome, sequence, int(sequence_position))`.
  - `InitiateDatabase` and `ModifyDatabase.add_annotation`: write both new columns.
- [ ] **Step 5: Run** the six files. Expected: PASS.
- [ ] **Step 6: Commit** "Store SNPs as sequence and position on that sequence".

### Task 2: Readers and result tables (minimap2, MUMmer4, shared base)

**Files:**
- Modify:
  - `CanSNPer2/modules/AlignmentReader.py`.
  - `CanSNPer2/modules/ParsePAF.py`: remove the offsets, `read_sequence_table` and `write_sequence_table`.
  - `CanSNPer2/modules/ParseDelta.py`: remove the offsets; read only the query fasta.
  - `CanSNPer2/modules/Aligners.py`: remove `Minimap2.prepare`.
  - `CanSNPer2/modules/CanSNPer2.py`: `SNP_FILE_COLUMNS` (:41); remove `add_reference_sequence` (:44-51), `with_reference_sequence` (:449-461), `reference_files` and `reference_tables` (:128-129) and their call sites (:975, :1013, :1025); the edge index (:1038).
- Test: `tests/test_parse_paf.py`, `tests/test_parse_delta.py`, `tests/test_aligners.py`, `tests/test_multiproc_queue.py`, `tests/test_dashboard.py`, `tests/test_annotations.py`, `tests/test_unique_names.py`, `tests/test_download_inputs_timeout.py`.

**Interfaces:**
- Consumes: `get_snps` from Task 1.
- Produces, in `AlignmentReader`:
  - `self.snplist = {(sequence, position): (ancestral, derived, snp_id)}` and `self.snp_positions = {sequence: [sorted int]}`.
  - `snps_in(sequence, start, end) -> [int]`. It returns `[]` for an unknown sequence.
  - `observe(sequence, position, base, block_offset, block_length, query_name="-", query_pos="N/A")`.
  - `finish()`, which builds rows `[snp_id, reference, sequence, str(position), query_name, query_pos, ancestral, derived, base, edge]` in `(sequence, position)` order.
- Produces, in `AlignmentReader.py`, shared with Task 3:
  - `SEQUENCES_SUFFIX = ".sequences"`.
  - `write_sequence_table(fasta, path)`, moved from ParsePAF.
  - `read_sequence_table(path) -> [(name, length)]`, which raises `FileNotFoundError` with the realign hint.
- Produces, in CanSNPer2.py: `SNP_FILE_COLUMNS` (Global Constraints) and `EDGE_COLUMN = SNP_FILE_COLUMNS.index("Edge flag")`, used at :1038.

Steps:

- [ ] **Step 1: Change the fixtures** in the listed test files to `helpers.make_db` with sequence names that match the fixture fasta headers.
- [ ] **Step 2: Write the failing tests:**
  - `test_paf_reads_target_coordinates_without_a_table`: a PAF record on `chr2` with no `.sequences` file reads the SNP at chr2:12.
  - `test_delta_reads_reference_coordinates_without_the_reference_file`: delete the reference fasta named in the delta header. The query base is still found.
  - `test_same_position_on_two_sequences_both_called`: for both readers, chr1:50 and chr2:50 get separate rows.
  - `test_result_columns`: the header of `_snps.txt` equals `SNP_FILE_COLUMNS`, and a called row has `row[2] == sequence` and `row[3] == position`.
  - `test_dashboard_keeps_edge_reasons` (Review Focus 1): a SNP observed with `block_offset=3` appears in the dashboard row's edge list.
- [ ] **Step 3: Run** the listed files. Expected: FAIL.
- [ ] **Step 4: Implement the interfaces.**
  - `ParsePAF.read_record`: use `start = int(fields[7]) + 1`, `end = int(fields[8])`, and `wanted = self.snps_in(fields[5], start, end)`.
  - `ParseDelta.read_record`: use `self.snps_in(rname, s1, e1)`, keeping its local coordinates.
- [ ] **Step 5: Run the full suite.** It must be green except `tests/test_parse_xmfa.py`, which Task 3 fixes.
- [ ] **Step 6: Commit** "Read minimap2 and MUMmer4 alignments per reference sequence".

### Task 3: progressiveMauve reader with a sequence table

**Files:**
- Modify: `CanSNPer2/modules/ParseXMFA.py` and `CanSNPer2/modules/Aligners.py` (add `ProgressiveMauve.prepare`).
- Test: `tests/test_parse_xmfa.py`.

**Interfaces:**
- Consumes: `write_sequence_table`, `read_sequence_table` and `SEQUENCES_SUFFIX` (Task 2). `prepare()` is already called with the staged, cleaned reference that Mauve aligns (`CanSNPer2.py:497-500`).
- Produces, in ParseXMFA:
  - `self.joined`, the sorted joined positions of the SNPs.
  - `self.joined_keys = {joined: (sequence, position)}`.
  - Both are built in `read_alignment(path)` from the table at `path + SEQUENCES_SUFFIX`, with offsets running in table order. SNPs on sequences missing from the table are left out.
  - `joined_snps_in(start, end) -> [joined]`.
  - `find_snps` uses `joined_snps_in` and calls `observe(*self.joined_keys[start + i - 1], ...)`.

Steps:

- [ ] **Step 1: Write the failing tests:**
  - `test_xmfa_maps_joined_coordinates_to_sequences`: a two-sequence table (chr1 600, chr2 400) with an xmfa block at 590-620 holding SNPs at chr1:595 and chr2:15. This test covers Review Focus 2.
  - `test_xmfa_reverse_block_on_second_sequence`: a "-" reference block on chr2. The complemented base is recorded at the right chr2 position.
  - `test_xmfa_without_sequence_table_fails_with_realign_hint` (Review Focus 3): run `run()` with no table. It raises `FileNotFoundError`, and the message contains "align again".
  - `test_mauve_prepare_writes_the_table`: `ProgressiveMauve().prepare(ref, q, out)` writes `out + ".sequences"` equal to `read_sequence_table` of the reference.
  - Update the existing fixtures so each xmfa has its table.
- [ ] **Step 2: Run** `pytest tests/test_parse_xmfa.py`. Expected: FAIL.
- [ ] **Step 3: Implement the interfaces.**
- [ ] **Step 4: Run the full suite.** Expected: PASS.
- [ ] **Step 5: Commit** "Map progressiveMauve coordinates to reference sequences".

### Task 4: Reference check per sequence

**Files:**
- Modify:
  - `CanSNPer2/modules/ReferenceCheck.py`.
  - `CanSNPer2/modules/CanSNPer2.py:463-479` (`check_reference_file`).
  - `CanSNPer2/DownloadGenomes.py:66-140,151-175`.
- Test: `tests/test_reference_sequences.py`. Rewrite it; its joined-position premise is gone.

**Interfaces:**
- Produces:
  - `reference_snps(database) -> {genome: [(snp_id, sequence, position, ancestral, derived)]}`.
  - `read_sequences(path) -> [(name, sequence)]`, with upper-case sequences.
  - `check_reference(path, name, snps=None) -> str | None`. It checks, in this order:
    1. unique and non-empty names (the existing `duplicate_names`);
    2. every sequence named in `snps` is present. Otherwise it returns a reason listing the missing names and the names in the file, and saying that RefSeq-named files must be renamed to GenBank accessions or fetched with CanSNPer2-download.
    3. each sequence is at least as long as its largest SNP position;
    4. the base rule of PR #67, with `MISMATCH_FRACTION` counted over all SNPs of the reference.
- Removes: `last_snp_positions`, `joined_sequence`, `sequence_length`, `locate`, `too_short`'s joined wording, and the `last_position` parameters in `DownloadGenomes.download`.

Steps:

- [ ] **Step 1: Write the failing tests:**
  - `test_reference_in_another_order_passes_and_types_the_same`.
  - `test_reference_with_an_extra_plasmid_passes`.
  - `test_header_with_description_matches_by_first_word`. Together with the two tests above, this covers Review Focus 4.
  - `test_missing_sequence_is_refused_naming_it`: the message contains `chrII`.
  - `test_short_sequence_is_refused`.
  - `test_wrong_bases_on_one_sequence_refused_at_fraction`.
  - `test_few_wrong_bases_warn`.
  - `test_download_refuses_a_genome_without_a_named_sequence`, using `DownloadGenomes.retrieve` monkeypatched as in the current file.
- [ ] **Step 2: Run** `pytest tests/test_reference_sequences.py`. Expected: FAIL.
- [ ] **Step 3: Implement the interfaces.** Update the two callers.
- [ ] **Step 4: Run the full suite.** Expected: PASS.
- [ ] **Step 5: Commit** "Check reference files per sequence".

### Task 5: Conversion command

**Files:**
- Create: `CanSNPer2/modules/ConvertPositions.py` and `tests/test_convert_positions.py`.
- Modify: `CanSNPer2/SNPDatabase.py`, adding the `--convert` and `--refdir` options. Handle them before the database checks, as `--import` is (:118).

**Interfaces:**
- Consumes: `read_sequences`, `check_database(current_format=False)`, and `DownloadGenomes.download(genome_id, assembly, refid, source="genbank", directory=refdir, snps=None)` for genomes missing from `refdir`.
- Produces:
  - `class ConversionError(Exception)`.
  - `convert_table(source, target, genome_files) -> int`.
    - Input is any tab file with `genome`, `position`, `derived_base` and `ancestral_base` columns, which covers `_snp.txt` and `_canSNPs.txt`. The `position` column is replaced in place by `sequence` and `sequence_position`; the other columns and comment lines are kept.
    - It raises `ConversionError` listing every row whose base on `genome_files[genome]` is neither derived nor ancestral, or whose position is beyond the end. Nothing is written in that case.
  - `convert_database(source_db, target_db, genome_files) -> int`. It copies to `target_db`, then rebuilds `snp_annotation` with the new columns through create, copy, drop and rename. It does not use `ALTER TABLE DROP COLUMN`, because the sqlite version is not fixed.
  - `genome_files(refdir, genomes, references=None) -> {genome: path}`. It uses `<refdir>/<genome>.fna` and downloads what is missing, using the accessions from the database's `snp_references` or from the `--references` file.
- CLI forms:
  - `CanSNPer2-database --convert --annotation OLD_snp.txt --references OLD_references.txt --refdir DIR -o OUTDIR` writes `OUTDIR/<basename>`.
  - `CanSNPer2-database --convert --database old.db --refdir DIR -o OUTDIR` writes `OUTDIR/old.db`.
  - Either form is refused if the output exists or is the input (Review Focus 5).

Steps:

- [ ] **Step 1: Write the failing tests:**
  - `test_converts_snp_file_on_two_sequences`: joined position 650 becomes `chr2 50`.
  - `test_converts_canSNPs_file_keeping_rank_and_score`.
  - `test_refuses_when_reference_is_in_another_order`: nothing is written.
  - `test_converts_database_and_types_the_same`: build an old database with `helpers.make_db(joined=True)` and convert it. `get_snps` then gives the expected keys, and `check_database` passes.
  - `test_refuses_to_overwrite`.
  - `test_downloads_missing_reference_once`, with a monkeypatched `retrieve`.
- [ ] **Step 2: Run** `pytest tests/test_convert_positions.py`. Expected: FAIL.
- [ ] **Step 3: Implement the interfaces.**
- [ ] **Step 4: Run the full suite.** Expected: PASS.
- [ ] **Step 5: Commit** "Add CanSNPer2-database --convert for joined positions".

### Task 6: RefSeq downloads renamed to GenBank names

**Files:**
- Modify: `CanSNPer2/DownloadGenomes.py`.
- Test: `tests/test_download_inputs_timeout.py`, or a new `tests/test_refseq_names.py`.

**Interfaces:**
- Produces:
  - `assembly_report_url(genome_url) -> str`, giving the `_assembly_report.txt` in the same NCBI folder.
  - `genbank_names(report_text) -> {refseq_accn: genbank_accn}`. It uses the `GenBank-Accn` and `RefSeq-Accn` columns of the report and skips rows whose GenBank value is `na`.
  - `rename_sequences(path, names)`, which rewrites the first word of each header in place and keeps the rest of the line. It is applied after unpacking and before `check_reference` when `source == "refseq"`.

Steps:

- [ ] **Step 1: Write the failing tests:**
  - `test_refseq_download_is_renamed_and_passes_check`.
  - `test_report_row_with_na_keeps_the_refseq_name_and_the_check_names_it`.
  - `test_genbank_download_does_not_fetch_a_report`.
  - All three use a monkeypatched `retrieve` serving the genome and the report.
- [ ] **Step 2: Run** the new tests. Expected: FAIL.
- [ ] **Step 3: Implement the interfaces.**
- [ ] **Step 4: Run** the new tests. Expected: PASS.
- [ ] **Step 5: Commit** "Rename RefSeq downloads to GenBank sequence names".

### Task 7: Docs, CI, release test, version

**Files:**
- Modify:
  - `README.md`. Rewrite the "References of several sequences" section for the new format. Add "Upgrading from 2.2", covering `--convert` and the refusal message. Update the result columns.
  - `.github/workflows/ci.yml`. Set the CanSNPer2-data checkout `ref:` and the `--branch` of the scheme download step to `per_sequence_positions` (the Part B branch). Add the comment "set back to master when 2.3.0 is released".
  - `tests/release_test.sh`: change its column count from 11 to 10, and point it at converted data.
  - `CanSNPer2/__init__.py`: `__version__ = "2.3.0"`.
  - Docstrings that still describe joined positions, in `Aligners.py:9-12`, `AlignmentReader.py:8-11`, `ParsePAF.py:6-10` and `ParseDelta.py:9-10`.
- [ ] **Step 1:** Run `grep -rn "joined\|end to end" CanSNPer2 README.md`. Only Mauve and conversion contexts may remain.
- [ ] **Step 2:** Run the full suite and the help test. Expected: PASS.
- [ ] **Step 3: Commit** "Document per-sequence positions; version 2.3.0".

---

## Part B: CanSNPer2-data

Create the worktree `CanSNPer2-data-perseq`, on branch `per_sequence_positions` from `check_reference_bases` (e38acd3). If #32 is merged first, use `origin/master` instead.

### Task 8: validate.py and check_bases.py read the new format

**Files:**
- Modify: `validate.py` (`SNP_HEADER` :29, `SNP_SET_HEADER` :438, `check_snp`, `check_snp_set`, the row index at :508), `check_bases.py` (`scheme_snps`, `check_scheme`).
- Test: `tests/conftest.py` (`SNP_HEADER`, rows), `tests/test_check_bases.py`, and the validate tests.

Steps:

- [ ] **Step 1: Write the failing tests:**
  - `test_old_header_names_the_conversion_command`.
  - `test_uniqueness_on_genome_sequence_position`.
  - In `test_check_bases.py`, put the fixture SNPs at `seq1:100`, `seq1:200`, `seq2:50` and `seq2:150`, and rename the swap test to `test_a_reference_with_its_sequences_in_another_order_passes`.
  - `test_a_missing_sequence_is_a_problem`.
- [ ] **Step 2: Run** `pytest`. Expected: FAIL.
- [ ] **Step 3: Implement.**
  - `check_scheme` looks up the base by sequence name. A missing sequence is a problem.
  - The report line keeps its form: `G: seq1 250 bp (2 SNPs), seq2 250 bp (2 SNPs)`.
- [ ] **Step 4: Run** `pytest`. Expected: PASS.
- [ ] **Step 5: Commit.**

### Task 9: Convert and rebuild the master schemes

The schemes on master are bacillus_anthracis, brucella_classical, coxiella_burnetii, francisella_tularensis and yersinia_pestis.

- [ ] **Step 1:** Install CanSNPer2 from the `CanSNPer2-perseq` worktree in editable mode. For each scheme, run `CanSNPer2-database --convert` on its `_snp.txt`, and on its `_canSNPs.txt` if present, with `--references` and a shared `--refdir` (for example `.references-fna/`, added to .gitignore). Replace the files. Expected: every scheme converts with no refused row.
- [ ] **Step 2:** In `build.py`, set `MINIMUM_CANSNPER2 = (2, 3, 0)`. Then run `python build.py --write`. Expected: databases and manifests are rewritten. `database/archive/` is left as is (old versions).
- [ ] **Step 3:** In `.github/workflows/validate.yml` and `release.yml`, install CanSNPer2 with `@dev_2.3.0` and add the comment "set back to master when 2.3.0 is released".
- [ ] **Step 4:** Update `schemes/README.md`, `schemes/TEMPLATE.md`, the `brucella_classical` README Origin section and `README.md` to describe the per-sequence positions.
- [ ] **Step 5:** Run `python validate.py && python check_bases.py && pytest`. Expected: all pass, with check_bases reporting the SNPs per sequence (brucella_classical: 22 on chrI and 2 on chrII).
- [ ] **Step 6:** Commit "Store SNP positions per reference sequence".

### Task 10: Dev-only schemes (after the release merges, Part D)

The dev-only schemes are brucella_anthropi, brucella_suis and francisellaceae.

- [ ] **Step 1:** Create a branch from `origin/dev` after master has been merged into dev. Convert those schemes as in Task 9. Expected SNP counts on chrII: brucella_suis 2 of 31, brucella_anthropi 7 of 42.
- [ ] **Step 2:** Run validate, check_bases and pytest. Commit. Open a PR into dev.

---

## Part C: Verification (before any PR)

The harness is in the scratchpad at `/private/tmp/claude-501/-Users-andreassjodin-Code-cansnp-test/994ce56b-cf43-4b50-8fdf-f189d18e38a2/scratchpad/mc/`. Copy it to `/Volumes/sekvens2/genomedata/audit/harness/` first, because the scratchpad can be cleared. The 2026-11-20 results in `mc/out/` are the old-code baseline. The audit databases are in `/Volumes/sekvens2/genomedata/audit/db`.

- [ ] **V1:** Convert the audit databases with `CanSNPer2-database --convert --database ... -o /Volumes/sekvens2/genomedata/audit/db_v23/`.
- [ ] **V2:** Make a copy of `run.sh` that takes the code folder and the database folder as variables. Mount `CanSNPer2-perseq` at `/code` for Mauve; Docker volumes go on short container paths. Run all three aligners on the 52 Brucella genomes, both the original and the reordered and reverse-complemented queries.
- [ ] **V3:** Adapt `compare.py` to join the old and new results on (sample, SNP Name). Expected: per SNP the same target base and the same call as the baseline, with 0 differences for each aligner.
- [ ] **V4:** Type with the reference sequences reordered (the `mc/refswap` references). Expected: identical to V3. The old code lost most calls here.
- [ ] **V5:** Use a reference with chromosome II removed. Expected: the run stops before aligning, with an error naming the missing accession.
- [ ] **V6:** Take 10 genomes each of F. tularensis, Y. pestis and B. anthracis from the scheme evaluation pool, and type them with old and new code on all three aligners. Expected: identical calls and bases.
- [ ] **V7:** Run the release gate `tests/release_test.sh` with progressiveMauve available, for example in the Docker image. Expected: pass. The 4 fast-tier failures known from master must not appear when the references are present.
- [ ] **V8:** Push both branches. Check that CanSNPer2 CI (data checkout `per_sequence_positions`) and CanSNPer2-data CI (CanSNPer2 `@dev_2.3.0`) are both green.

## Part D: Integration (each step needs the user's go-ahead)

1. Open PR `per_sequence_positions` into CanSNPer2 `dev_2.3.0`, and PR the CanSNPer2-data branch into its master. Do not merge yet.
2. At release, in this order and back to back:
   1. In CanSNPer2 `dev_2.3.0`, set the CI data ref back to master. Merge `dev_2.3.0` into master and tag 2.3.0.
   2. In the CanSNPer2-data PR, set the install back to `@master`, then merge it.
   3. Merge CanSNPer2-data master into dev, then do Task 10.
3. CanSNPCreator export and audit get their own plan after the release (spec "Order of work", step 5).

**Execution:** I recommend the subagent-driven method. Tasks 1-5 share interfaces, and a missed key or column index gives wrong calls without an error, so a fresh reviewer per task is worth the cost. Tasks 1 and 2 are reviewed as a pair.

After approval, copy this plan to `CanSNPer2-perseq/docs/superpowers/plans/2026-11-21-per-sequence-positions.md` and commit it there.
