# Validation

These are internal evaluations made by the developers of CanSNPer2 while the features were
developed; they are not a peer-reviewed validation. They show how the aligners and read typing
behaved on the datasets named below, and they are given so that users can judge the results on
their own data. Each section states what was compared, the data, the software and the date.

None of these evaluations uses an external truth at SNP level. Read typing is compared with the
calls CanSNPer2 makes on the assemblies of the same genomes, and the aligners with each other.

## Aligners on Francisellaceae genomes

**Compared.** Typing with minimap2 and with MUMmer4 against typing with progressiveMauve
(progressiveMauve was the default aligner in earlier releases): final calls, and the base read
at each canSNP position.

**Data.** The 1,205 genomes of the family Francisellaceae in GTDB release 232: 937 genomes of
*F. tularensis*, and 268 genomes of other Francisellaceae species and genera. A subset of 77 of
the 268 (43 Francisella genomes of other species, 34 of other genera) was also typed with
progressiveMauve.

**Software and settings.** CanSNPer2 development version 2.2.2 with the minimap2 and MUMmer4
aligners as they were proposed (pull requests #58 and #61), the `francisella_tularensis`
scheme, `--min-required-hits 3`, `--strictness 0.7`. progressiveMauve was the Bioconda
package `progressivemauve` in a container with CanSNPer2 2.2.2. The tool versions were not
recorded with the results; the environments used contained minimap2 2.31 and MUMmer 4.0.1 when
this page was written. MUMmer4 `nucmer` ran with its default of two threads; `-t 1` was added
later (2026-11-21). The positions were read as positions along the joined reference sequences,
before per-sequence positions (2.3.0).

**Date.** 2026-11-17 and 2026-11-18.

**Results, 937 F. tularensis genomes.** Both aligners gave the same final call as
progressiveMauve for every genome, and no canSNP position had a different base.

| | minimap2 | MUMmer4 |
|---|---|---|
| positions with the same base as progressiveMauve | 305,513 | 305,693 |
| positions with data only from progressiveMauve | 291 | 111 |
| positions with data only from this aligner | 169 | 177 |
| time for the 937 genomes, three parallel runs on one computer | 12.9 min | 15.0 min |

MUMmer4 and minimap2 agreed on every call and base; MUMmer4 had data at 222 positions where
minimap2 had none, and minimap2 at 34 positions where MUMmer4 had none.

**Results, other Francisellaceae.** Neither minimap2 nor MUMmer4 gave a final call for any of the
268 genomes. progressiveMauve, run on the subset of 77, called 5 genomes of other genera
*F. tularensis* (two at `B.1`, one at `B.2`, two at `M.1`). minimap2 does not align genomes of
other genera; MUMmer4 does, at about 43 s per genome against under 4 s for minimap2. At 41
canSNP positions in divergent Francisella species, MUMmer4 and progressiveMauve read different
bases. These positions lie next to 1 bp gaps in alignments of 80-86 % identity, where the base
depends on where the gap is placed; no call changed.

**minimap2 secondary chains.** The minimap2 options `-N 50 -p 0.1` were added after the first
comparison. On the 208 genomes of the 937 where at least one aligner lacked data at some canSNP,
the positions with data only from progressiveMauve went from 351 to 291; 53 of the 60 recovered
positions were confirmed by a unique reference flank in the query. No position changed its base
and all 208 final calls were unchanged; the run time was about 10 % longer.

## Read typing on simulated reads

**Compared.** The final SNP from simulated reads against the final SNP of the assembly of the
same genome, at several depths and for three platforms; and the possibly-mixed flag on single
genomes and on simulated mixtures.

**Data.** 50 genomes from GTDB, 30 of the genus Brucella and 20 of *F. tularensis*, typed with
the schemes that cover them (the Brucella genomes against the reference genome *B. melitensis*
16M). Reads were simulated with wgsim (Illumina), badread (ONT) and as 15 kb reads with 0.1 %
errors (HiFi). Mixtures were made from two genomes on different branches at 50/50, 80/20 and
95/5.

**Software and settings.** The read typing branch of CanSNPer2 that became 2.6.0, minimap2 2.31,
samtools 1.24, the platform defaults of 2.6.0 (the ONT `--min-depth` default was lowered from 10
to 5 after this evaluation).

**Date.** 2026-11-24.

**Results.** No read call was on another branch than the assembly call, at any depth or
platform (430 runs). With fewer reads the call stopped higher in the tree or was `NA`:

| | 5x | 10x | 20x | 30x | 50x |
|---|---|---|---|---|---|
| Illumina, same call as the assembly | 0 of 50 | 37 of 50 | 40 of 50 | | 41 of 50 |
| ONT | | 42 of 50 | | 50 of 50 | |
| HiFi | | 46 of 50 | | 50 of 50 | |

The Illumina misses from 20x on are Brucella genomes outside the classical cluster, far from
the reference genome 16M, whose reads mapped at too few canSNPs. With a `--min-depth` of 10,
ONT reads at 10x gave no call for 42 of 50 genomes.

**Mixed sites.** Among the 50 single genomes, one (an *F. tularensis* draft with 90-96 %
identical copies of three canSNP regions on short contigs) was flagged with 3 to 5 mixed sites;
the other 49 had at most 2. Mixtures of two genomes had 8 to 17 mixed sites at 50/50, 3 to 11
at 80/20 and 0 to 3 at 95/5. No mixture was called on another branch than its majority genome;
the call stopped higher in the tree or was that of the majority genome.

## Local references

**Compared.** The final SNP from 50x Illumina reads mapped only to the scheme reference genome,
and with `--local-references`, against the final SNP of the assembly.

**Data.** The `brucella_genus` scheme (283 typed canSNPs, reference genome
*B. melitensis* 16M, 221 representative genomes). The nine Brucella genomes of the read typing
evaluation whose 50x Illumina reads missed their assembly call, all outside the classical
cluster. Five of them are themselves representatives; for these, the genome was taken out of
the local references, so that the reads were not mapped to their own genome ("held out" below).
The other four are members of a representative at 99 % average nucleotide identity or closer.

**Software.** The local references branch of CanSNPer2 that became 2.6.0, minimap2 2.31,
samtools 1.24.

**Date.** 2026-11-24 and 2026-11-25.

**Results.**

| genome | assembly call | 16M only | with local references |
|---|---|---|---|
| GCA_039674445.1 | `BR.chromiisoli` | NA | `BR.chromiisoli` |
| GCF_026256505.1 | `BR.intermedia.1.4` | NA | `BR.intermedia.1.4` |
| GCA_031999685.1 | `BR.intermedia.1.1.1` | `BR.1.1.1.1.1.1` | `BR.intermedia.1.1.1` |
| GCF_000182645.1 (held out) | `BR.intermedia.1.3.1` | NA | `BR.intermedia.1.3.1` |
| GCF_008932555.1 | `BR.intermedia` | NA | `BR.intermedia` |
| GCA_008932375.1 (held out) | `BR.haematophila` | `BR.1.1.1.1.1.2` | `BR.haematophila` |
| GCF_003664555.1 (held out) | `BR.teleogrylli.1.1` | `BR.1.1.1.1.1.2` | `BR.teleogrylli` |
| GCA_018798825.1 (held out) | `BR.sp018798825` | NA | `BR.sp018798825` |
| GCF_003049685.1 (held out) | `BR.pituitosa.1.1.2` | NA | `BR.pituitosa.1.1.2` |

Eight genomes reached the assembly call and one stopped a node above it; none was called on
another branch. The canSNPs with reads rose from 60-158 to 175-259 of 283. The counts of the
first round chose the local reference in eight samples without a screen; one needed three
screens. One sample (*B. pituitosa*) went first to a genome of another species and reached its
own lineage in a second round.

**Real reads.** Two Illumina read sets of *B. anthropi* isolates (internal data, not public),
typed with `brucella_genus` against 16M only, had reads at 169 and 171 of the 283 canSNPs and
stopped one node above `BR.anthropi`. With `--local-references` and one local round each, both
were called `BR.anthropi.1.1.1.1`, the call of the draft assembly of one of the two isolates.
