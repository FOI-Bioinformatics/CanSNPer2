# Troubleshooting

This page lists the exit statuses of `cansnper2`, the common error messages with their cause and
what to do, the warnings that are expected in normal runs, and how to report a problem. It is
written for users who run CanSNPer2 by hand or in a pipeline.

## Exit status

Every verb checks its inputs before it writes or aligns anything, and a problem ends the command
with one message that names the file or database and says what to do.

| Status | Meaning |
|---|---|
| 0 | The command did what it was asked to do. Warnings may have been printed. |
| 1 | An input or a step could not be used; the message says which and why. `cansnper2 type` also exits with 1 when no sample of the run was typed (every sample was not called, failed or skipped) and no earlier result was reused. |
| 2 | A usage error found by the option parser: an unknown option, a missing verb, or options that cannot be combined (for example `--fetch-local` without `--local-references`, or `--data-dir` with `--repo` or `--branch`). Also an unexpected internal error, printed as `CanSNPer2 stopped: <type>: <message>; run with --debug for details`. |
| 130 | `cansnper2 type` was interrupted (Ctrl-C, or SIGTERM during the alignments, as sent by a cluster scheduler). The running aligner processes are stopped and the temporary folder is removed unless `--keep-temp` is given. |

A run with `--keep-going` exits with 0 when at least one sample was typed, also when other
samples failed; check `summary.html` or the log for the failed samples.

A database is never left half written: `db create`, `db import` and `scheme get` remove a
database they had started, and `db modify` and `db remove` work on a copy that replaces the
database only when every step succeeded.

## Where to look

- Warnings and errors are printed in the terminal, most of them on standard error, and written
  to the log of the command in `--log-dir` (default `logs/`), for example
  `logs/cansnper2-type-<Mon-DD-YYYY>.log`. The names of all log files are given in
  [output-formats.md](output-formats.md#log-files).
- By default the log holds warnings and errors. Run again with `-v` to log progress as well,
  or with `--debug` to log debug messages.
- The output of each aligner is in `<log-dir>/<reference>_<sample>.<tag>.log`. An alignment
  that fails is named there.
- For an unexpected error (exit status 2), run the same command with `--debug` to print the
  traceback. The traceback is also written to the log file of the command.

## Messages and what to do

The tables give the beginning of each message as it is printed; `X`, `N` and `DIR` stand for
names, numbers and folders.

### Typing (cansnper2 type)

| Message | Cause | What to do |
|---|---|---|
| `A database is required` | `-d` was not given. | Give the database with `-d`. |
| `No query files were given` | No query and no `--samplesheet`. | Give the query files, a list file or a samplesheet. |
| `The database X does not exist, fetch it with cansnper2 scheme get or create it with cansnper2 db create` | Wrong path. | Check the path, or get the database. |
| `The database X cannot be used. It stores SNP positions along the joined reference sequences` or `... Its snp_annotation table holds one canSNP per node` | A database from before 2.3.0, or from 2.3.0. | Convert it with `cansnper2 db convert`, or fetch the current scheme with `cansnper2 scheme get`. |
| `The database X has no canSNPs`, `The database X lists no reference genomes` | The database was built without its canSNP or references file. | Build it again with all scheme files. |
| `The database X has N canSNP(s) but none of them is marked typed` | No canSNP has `typed` `yes`. | Build it again from a canSNP file that marks one canSNP of each node `typed` `yes`. |
| `The database X has N SNP(s) on node ids that are not in the nodes or tree table` | The database is damaged. | Export it with `cansnper2 db export`, correct the files and build it again. |
| `The aligner X was not found` | The aligner is not in `PATH`: minimap2 by default, or the program chosen with `--aligner`. | Install it (for example minimap2 from Bioconda), give `--aligner-path`, or choose another `--aligner`. `cansnper2 selftest` lists the aligners it finds. |
| `Read typing needs minimap2 and samtools; X was not found` | The queries are reads and a program is missing. | Install minimap2 and samtools (for example from Bioconda). |
| `No reference .fna files were found in DIR` | `--ref-dir` is empty or wrong. | Run `cansnper2 ref get -d X -o DIR`, or give the right `--ref-dir`. |
| `N reference genome(s) of X have no file <genome>.fna in DIR` | Some reference genomes are missing; their canSNPs would not be compared. | Fetch them with `cansnper2 ref get -d X -o DIR`. To type without them on purpose, give `--allow-missing-references`. |
| `reference X (FILE) has no sequence named ...` | The sequence names of the file are not the names the database uses, for example RefSeq names for a scheme built on GenBank accessions. | Rename the sequences, or fetch the file with `cansnper2 ref get`. See [input-formats.md](input-formats.md#reference-genome-files). |
| `sequence S of reference X (FILE) has N bases, but the database places SNPs on it up to position P` | The file is truncated or is another genome. | Delete it and run `cansnper2 ref get` again. |
| `reference X (FILE) has neither the ancestral nor the derived base at N of M SNPs` | The file is not the genome, or not the version of its sequences, that the scheme was built on. | Delete it and run `cansnper2 ref get` again. |
| `reads and assemblies cannot be typed in one run` | FASTA and FASTQ files were given together. | Type them in separate runs. |
| `X is named as one mate of paired-end reads but its other mate was not given` | A file named `_R1`, `_2` and so on without its mate. | Give both files, or list the file in a `--samplesheet` to type it as single-end reads. |
| `X and Y give the same sample name Z` | Two read files map to one sample name. | Rename one, or use a `--samplesheet`. |
| `Samples are named after their file name, which must be unique` | Two query files with the same file name in different folders. | Rename one, or type them into separate `--outdir`. |
| `X is a query list (it does not start with a FASTA header) and a list must be the only query argument` | A list file was given together with other queries, or a query is not FASTA or FASTQ. | Put all paths in the list, or check the query file. |
| `X has no FASTA header: its lines are sequence, not file paths of a query list` | A FASTA file without its header line. | Add a header line. |
| `N query file(s) were not found or are not files` | Wrong paths; paths in a list file are read from the current folder. | Correct the paths. |
| `Error while processing S: ...` followed by `A query could not be typed and the run was stopped; use --keep-going to continue with the next query` | One query failed; the first message gives the reason, for example `X is gzip-compressed, decompress it first` or `query X has more than one sequence named ...`. | Correct that query, or give `--keep-going` to type the others. |
| `... did not finish within N minutes for genome S against reference R, it was stopped (--alignment-timeout)` | An alignment ran longer than `--alignment-timeout` (default 60 minutes). | Check the query; give a larger `--alignment-timeout`, or `0` for no limit. |
| `X failed with exitcode N for LOG` | The aligner failed. | Read the alignment log named in the message. |
| `FILE exists for a different input` | The output folder holds a result of another file with the same sample name. | Give `--rerun` to replace it, or another `--outdir`. |
| `FILE was made with a different database or settings` | The earlier result was made with another database or other settings. | Give `--rerun` to replace it, or another `--outdir`. |
| `No sample was typed - S: ...` | No sample of the run was called; exit status 1. | Read the reasons per sample; see [interpreting-results.md](interpreting-results.md). |

### Databases and schemes

| Message | Command | What to do |
|---|---|---|
| `The source files cannot be built into a database, N problem(s):` followed by `FILE:LINE: ...` | `db create`, `db import`, `scheme get` | Correct the listed rows; all problems are listed together. [input-formats.md](input-formats.md#checks-before-a-database-is-built) gives the rules. |
| `FILE:LINE: genbank_id 'X' is not an assembly accession` | as above, and `db modify --references` | Give the GenBank (`GCA_`) or RefSeq (`GCF_`) assembly accession with its version, such as `GCA_000008985.1`. |
| `FILE:LINE: genbank_id 'X' of genome 'Y' is already used by genome 'Z'` | as above | Give each genome its own assembly accession. |
| `The genomes file cannot be loaded, N problem(s):` | `db create`, `db import`, `scheme get` | Correct the genomes file; see [input-formats.md](input-formats.md#genomes-file). |
| `The source file FILE cannot be read: it is not valid UTF-8 text` | all that read scheme files | Save the file as UTF-8. |
| `The database X already exists; db create builds a new database` | `db create` | Choose a new name, or remove the file first. |
| `X already exists, use --overwrite to replace it` | `scheme get` | Give `--overwrite`, or another `-d`. |
| `The files cannot be added to X, N problem(s):` | `db modify` with `--cansnps` or `--references` | Correct the listed rows. A canSNP on a node or genome that the database does not hold, or a genome, rank, typed canSNP or position the database already holds, is refused. The database was not changed. |
| `A modification needs the branch it applies to, give the branch root with --parent!` | `db modify --mod-file` | Give `--parent`. |
| `The selected parent node (X) could not be found in the source database!`, `The node X could not be found in the database`, `X is a root of the tree` | `db modify`, `db remove` | Give a node of the database with `--parent`; the root cannot be removed. |
| `the modification would leave X with a tree that is not valid` | `db modify` | The modification gives a node a second parent or makes a cycle; correct the modification file. |
| `... could not be written to X` | `db create`, `db modify`, `db remove` | The disk is full or the database cannot be written; free space or check the permissions. |
| `Could not list the files of REF of OWNER/NAME: ... HTTP Error 403: rate limit exceeded` | `scheme list` | The GitHub API allows 60 anonymous requests per hour from one address, which shared networks reach. Set the environment variable `GITHUB_TOKEN` (or `GH_TOKEN`) to a GitHub token, or try again later. See [schemes-and-databases.md](schemes-and-databases.md#github-api-limits). |
| `There is no scheme X in ...` | `scheme get` | Check the name with `cansnper2 scheme list` (draft schemes: `--branch dev`). |
| `scheme X on OWNER/NAME@REF uses positions on the joined reference (before 2.3.0)` | `scheme get` | The branch holds the format before 2.3.0. Use another branch or commit. |
| `... was not found in ...; the files of the verified scheme X cannot be checked`, `... records no hash for ...` | `scheme get` | The files cannot be checked against the manifest of the repository. Report it to the repository, or give `--allow-unverified` when the source is trusted. |
| `FILE does not match the SHA-256 in MANIFEST` | `scheme get` | The download is damaged, or the scheme was changed without a rebuild of its manifest. Try again later or report it to the repository. |
| `The database X could not be built from the scheme Y`, `The database X could not be written` | `scheme get` | A database or file system error; free space or check the permissions of the folder. No partial file is left behind. |
| `X lists no reference genomes, there is nothing to download` | `ref get` | Build the database again with its references file. |
| `REF: the database has no genbank_id for this genome` | `ref get` | Give the assembly accession in the references file and build the database again. |
| `N of M references could not be downloaded:` | `ref get`, `scheme get` | The lines below give the reason per genome (network, NCBI, or a downloaded file that does not fit the canSNPs). Run the command again; genomes already downloaded are kept. |
| `The SNP positions of X could not be read` | `ref get`, `type` | The database is damaged; export it with `db export`, correct the files and build it again. |

## Expected warnings

These warnings appear in normal runs and do not change the exit status:

- `N warnings about the results are in the log FILE, run with --verbose to see them` (for one,
  `1 warning about the results is in the log FILE, run with --verbose to see it`): the
  warnings about canSNP positions without any alignment data
  (`N SNP positions in the tree had no alignment data at all, ...`) are written to the log only.
  Such positions are listed with `N/A` in `<sample>_not_called.txt`. A few positions without data
  are common for draft assemblies and for samples distant from the reference genomes.
- `N SNP positions within 25 bases of an alignment edge were found against R and are flagged
  'edge' in the result files`: these canSNPs are marked in the `Edge flag` column; a call whose
  path runs over one has the status `review` (see
  [method.md](method.md#positions-near-alignment-edges-and---mask)).
- `Skipping N reference file(s) in DIR that are not reference genomes with SNPs in the
  database`: other genomes in a shared reference folder are not aligned against.
- `Could not resolve REF to a commit, the files are read from REF directly and the commit is not
  recorded: ...`: the GitHub API refused the request (often the rate limit above). The database
  is built, with the commit `unknown`; build it again later, or with `GITHUB_TOKEN` set, to
  record the commit (see [reproducibility.md](reproducibility.md#the-sourcejson-file)).
- `Scheme of X: ... The scheme is a draft and has not been verified`, `... The status of the
  scheme is unknown`, `... The database was changed after it was built`: the database is not a
  verified published scheme as built. See [reproducibility.md](reproducibility.md).
- `reference X has neither the ancestral nor the derived base at N of M SNPs: ...; the scheme
  may give a wrong base for these SNPs`: fewer than a fifth of the canSNPs of a reference disagree
  with it.
- `N tree annotations name SNPs that are not in this database and were skipped`: lines of the
  `--node-labels` file name nodes that the database does not hold.
- `X is given more than once, it is typed once`.
- `S is possibly mixed: at N canSNP sites no base reached the allele fraction`: reads only; see
  [interpreting-results.md](interpreting-results.md).
- `--mask is not used for reads`.

## Reporting a problem

Report problems at the issue tracker of the repository,
<https://github.com/FOI-Bioinformatics/CanSNPer2/issues>. Please include:

- the output of `cansnper2 --version` and `cansnper2 selftest`;
- the full command line;
- the messages, and the log file of a run with `--debug`;
- the output of `cansnper2 db info -d <database>` when a database is involved;
- the operating system and how CanSNPer2 was installed (Bioconda or pip).

When the data cannot be shared, a small example that shows the same problem is helpful.

```sh
cansnper2 --version
cansnper2 selftest
cansnper2 db info -d francisella_tularensis.db
```
