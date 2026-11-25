# Interpreting results

This page is for users who read CanSNPer2 results: what the final SNP, the support and the call
status mean, what to check before using a result, and what to record in a report. The exact
rules are in [method.md](method.md); the columns of each output file are in
[output-formats.md](output-formats.md).

## The final SNP

The final SNP is the deepest node of the scheme tree whose path from the root is confirmed by
the canSNPs of the sample. It names a lineage of the scheme: the sample carries the derived
allele of that node and of the nodes above it. It does not say that the sample belongs to any
node below it; the canSNPs below the final SNP were ancestral, had no data, or did not meet the
rules (see [method.md](method.md#tree-walk-and-the-final-snp)).

A final SNP is only as specific as the scheme. A scheme resolves the lineages it has nodes for,
and a sample from a lineage the scheme does not cover is placed at the deepest node of the
scheme that it shares.

## Support

Support is written as `found/expected`: the number of derived canSNPs on the path from the root
to the final SNP, out of the number of canSNPs on that path. `12/12` means every canSNP on the
path was derived. `11/12` means one canSNP on the path had no data; the call relies on the
`--strictness` tolerance (default 0.7) for it. A canSNP on the path is never ancestral or an
other base in a called sample: either breaks the path.

Support is a count, not a probability. A short path (a node near the root) has a small
`expected` and gives less evidence than a long one, whatever the fraction.

## Call status

`summary.html` shows a status for each sample and lists the reasons for it.

| Status | Meaning | What to do |
|---|---|---|
| `clean` | A final SNP was called, every canSNP on its path has data, none is flagged at an alignment edge, no derived canSNP lies outside the path, and (for reads) the sample is not flagged as possibly mixed. | Use the call, with the checks under [When to treat a result as unreliable](#when-to-treat-a-result-as-unreliable). |
| `review` | A final SNP was called, but at least one of the conditions above does not hold. | Read the reasons, and look at the sample in the tree of `summary.html` and in its result files (sections below). |
| `not called` | No final SNP (`NA`). | See [NA results](#na-results). |
| `failed` | The alignment, the mapping or the processing of the sample failed. | Read the log of the sample in `--log-dir`. |
| `reused` | The result was read from an earlier run in the same output folder. | Type again with `--rerun` to see the sample in the tree. |

## Reasons for review

### canSNPs on the path without data

The call depends on the strictness tolerance. Check where the positions without data lie: a
canSNP without data directly above the final SNP is not allowed, so the missing ones are higher
up the path. Positions without data arise from gaps between the contigs of a fragmented
assembly, from regions that do not align to the reference genome, and, for reads, from low
depth. Several missing canSNPs on a short path weaken the call.

### Edge flags

A canSNP within 25 bases of the end of an aligned block is flagged `edge` in the result files
(assemblies only). Contig ends and block ends are where alignment errors are most likely. An
edge-flagged canSNP on the path sets the status to `review`. Check whether the call still holds
without it; typing again with `--mask 25` leaves such canSNPs out (they then count as no data)
and shows whether the call depends on them.

### Derived canSNPs outside the path

A derived canSNP on another branch than the SNP path is not expected in a single isolate typed
with a correct scheme. Possible causes are an alignment error at that position (for example a
paralogous region), a canSNP that is not canonical in the scheme, or a mixed sample. The reasons in `summary.html`
name the canSNPs; several on one other branch need a closer look than a single one.

### Possibly mixed (reads)

A read sample with at least `--mixed-sites` (default 3) canSNP sites where no base reached
`--min-allele-fraction` is flagged as possibly mixed. The flag never changes the call. Two
causes give this pattern:

- the sample holds more than one strain; the call then follows the canSNPs where one base
  dominates and may stop higher in the tree than the majority strain would;
- the sample has a region twice where the reference genome has it once, so reads from both
  copies map to one place.

The counts of each canSNP in `<sample>_site_counts.tsv` show which applies at each site. A
sample flagged as mixed should not be reported as a single lineage without this check.
CanSNPer2 does not separate the strains of a mixture.

## No-data sites

Every canSNP of the scheme appears in one of the two per-sample result files: derived canSNPs
in `<sample>_snps.txt`, the others in `<sample>_not_called.txt`, where canSNPs that no
alignment covered have the target base `N/A`. For reads, the `state` column of
`<sample>_site_counts.tsv` gives the reason for each canSNP (`no_reads`, `low_depth`, `mixed`,
`not_lifted`, `conflict`). Many no-data sites in one part of the tree, with data elsewhere,
point to a region the sample lacks or that does not align to the reference genome; no data
throughout points to a sample from outside the scheme or to an input problem.

## A final SNP shallower than expected

A call that stops above the expected lineage is a less specific call, not a wrong one, when the
canSNPs below it are without data rather than ancestral. Common causes:

- canSNPs below the call had no data (missing contigs, low read depth, reads far from the
  reference genome). For reads of lineages distant from the reference genome, see
  `--local-references` in [method.md](method.md#local-references).
- the canSNP directly above a deeper node had no data, which the rules do not allow;
- the sample belongs to a lineage the scheme does not resolve further.

Compare the tree of the sample in `summary.html`: open circles are canSNPs without data, grey
ones ancestral, blue ones an other base.

## NA results

`NA` means that no node met the rules. The reason is given in the log and in `summary.html`:

- fewer derived canSNPs in the sample than `--min-required-hits` (default 3), which is
  expected for a genome from outside the scheme's taxon and for a sample near the root of a
  scheme;
- no node confirmed: an ancestral or other base on the path, too few canSNPs with data
  (strictness), or no data directly above the deepest derived canSNP;
- two confirmed paths at the same depth: CanSNPer2 does not choose between them. The candidates
  are named in the log; in `summary.html` the reason currently reads "... was not confirmed
  with enough hits". Look at the sample in the tree.

`cansnper2 type` exits with status 1 when no sample of the run was typed. Lowering
`--min-required-hits` or `--strictness` gives calls on less evidence; when it is done, record
the values used.

## When to treat a result as unreliable

Treat a call with caution, or report it as unresolved, when:

- the scheme is a draft (`--branch dev`), its status is unknown, the database was modified
  after it was built, or the scheme lists known issues (`cansnper2 db info` shows them) and a canSNP of the path
  is among them;
  CanSNPer2 logs a notice in each case and `summary.html` shows it in its header;
- the sample is from outside the taxon the scheme was built for (see
  [limitations.md](limitations.md));
- the status is `review` and the reasons cannot be explained by the data;
- the support has a small `expected`, or several canSNPs on the path have no data;
- a read sample is flagged as possibly mixed, has a low median depth, or was typed with
  thresholds below the platform defaults;
- the result differs between aligners (see [validation.md](validation.md) for the differences
  observed so far).

## What to report

For each sample, report the final SNP, the support and the call status. For the run, record:

- the scheme name, its status, and the repository, branch and commit it was built from (the
  `Scheme:` line of `snps_summary.txt` and of the `summary.html` header);
- the sha256 of the database file (the `Database:` line, and the `Typing settings:` line of each
  `<sample>_snps.txt`);
- the CanSNPer2 version (the `CanSNPer2 version:` line);
- the aligner, or for reads the platform and read thresholds (the `Aligner:` and
  `Typing settings:` lines), and the version of the aligner or of minimap2 and samtools, which
  CanSNPer2 does not record;
- `--min-required-hits`, `--strictness` and `--mask` when they differ from the defaults, and
  whether `--local-references` was used.

`snps_summary.txt` is written with `--summary`. See [reproducibility.md](reproducibility.md) for
keeping a run reproducible.
