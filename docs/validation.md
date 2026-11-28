# Validation

These are internal evaluations made by the developers of CanSNPer2 while the features were
developed; they are not a peer-reviewed validation. They show how the aligners and read typing
behaved on the datasets named below, and they are given so that users can judge the results on
their own data. Each section states what was compared, the data, the software and the date.

None of these evaluations uses an external truth at SNP level. Read typing is compared with the
calls CanSNPer2 makes on the assemblies of the same genomes, and the aligners with each other.

## Aligners on Francisellaceae genomes

**Compared.** Typing with minimap2 and with MUMmer4 against typing with progressiveMauve
(progressiveMauve was the default aligner up to 2.4.0): final calls, and the base read
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
samtools 1.24, the platform defaults of 2.6.0 (`--min-depth` 5 for every platform). The ONT
default was lowered from 10 to 5 as a result of this evaluation; the table below uses 5.

**Date.** 2026-11-24.

**Results.** No read call was on another branch than the assembly call, at any depth or
platform (430 runs). With fewer reads the call stopped higher in the tree or was `NA`:

| | 5x | 10x | 20x | 30x | 50x |
|---|---|---|---|---|---|
| Illumina, same call as the assembly | 0 of 50 | 37 of 50 | 40 of 50 | | 41 of 50 |
| ONT | | 42 of 50 | | 50 of 50 | |
| HiFi | | 46 of 50 | | 50 of 50 | |

The Illumina misses from 20x on are Brucella genomes outside the classical cluster, far from
the reference genome 16M, whose reads mapped at too few canSNPs. With the earlier ONT
`--min-depth` of 10, ONT reads at 10x gave no call for 42 of 50 genomes, where `--min-depth` 5
gave the assembly call for 42 of 50.

**Mixed sites.** Among the 50 single genomes, one (an *F. tularensis* draft with 90-96 %
identical copies of three canSNP regions on short contigs) was flagged with 3 to 5 mixed sites;
the other 49 had at most 2. Mixtures of two genomes had 8 to 17 mixed sites at 50/50, 3 to 11
at 80/20 and 0 to 3 at 95/5. No mixture was called on another branch than its majority genome;
the call stopped higher in the tree or was that of the majority genome.

## Platform detection and the clr preset

**Compared.** The platform read from the input (#119) against the platform of the file, and the
call from PacBio CLR reads with the `clr` preset against the call recorded for the same reads
with the earlier settings and against the published assembly.

**Data.** Seven read files of four kinds: an Illumina run from the SRA with its read names
replaced (SRR29342398), an Illumina run with the names the instrument wrote, two PacBio CLR
subread files of PRJNA991498 and three nanopore runs of short fragments (median 217-343 bases).
For the call, SRR29384771 of PRJNA991498 filtered to reads of at least 10 kb and subsampled to
50x, against the draft brucella_anthropi scheme.

**Date.** 2026-11-27.

**Results.** All seven files were reported correctly. The nanopore runs are shorter than the
threshold for long reads and were recognized by their read names; the two CLR files by their
read lengths and the absence of base qualities.

With `--platform auto` the CLR reads were typed with the `clr` preset (42 of 42 canSNP sites
with reads, median depth 40, one mixed site) and gave `Bant.anthropi.env.1.1.1.1.1.1`, the call
recorded for all eight isolates of the project with `--platform ont --min-baseq 0` and the call
of the published assembly GCA_030518775.1. The same file with the default `illumina` platform
warns twice, that the reads look like `clr` and that `min_baseq=20` discards every base of them,
and gives `NA` as it did before.

The thresholds of the `clr` preset were measured afterwards; the next section reports it.

**Open.** Detection has not been checked on public PacBio HiFi runs, which belongs to the
validation on real long reads that is still open. `--platform` remains the setting to use when
the platform is known, and the default is unchanged (`illumina`).

## Thresholds of the clr preset

**Compared.** The call and the canSNP states from PacBio CLR reads against the call CanSNPer2
makes on the assembly of the same reads, over four minimap2 presets, six read depths and a grid
of `--min-depth` and `--min-allele-fraction`.

**Data.** The eight *B. anthropi* CLR runs of BioProject PRJNA991498 (Sequel II subreads,
2.6-3.5 Gb per run, 500-700x of the 5.1 Mb genome, mean read length 10.4 kb, every quality
value Phred 0). A nested depth series of 5x, 10x, 20x, 30x and 50x was drawn from each run with
a seeded generator, so a read kept at one depth is kept at the deeper ones; the depths obtained
are 4.9-50.1x. No length filter was applied. The full read sets were typed as well, with every
read mapped (`--max-reads 0`).

Two schemes, both drafts on the `dev` branch of CanSNPer2-data: `brucella_anthropi` (42
canSNPs, reference genome *B. anthropi* ATCC 49188, the species of the samples) and
`brucella_genus` (283 canSNPs, reference genome *B. melitensis* 16M, another species of the
genus). The comparison is the Flye assembly of the same reads: all eight assemblies are called
`Bant.anthropi.env.1.1.1.1.1.1` (8 of 8 canSNPs on the path, `clean`) and `BR.anthropi.1.1.1.1`
(16 of 16, `review`, from derived canSNPs on `BR.classical.*` branches that follow from the
reference genome). The published genome of the project, GCA_030518775.1, from another
BioSample, gives the same two calls. The assemblies are made from the same reads, so this is
not an independent truth; it is the comparison the sections above use.

**Software and settings.** CanSNPer2 2.6.0 at commit b545d3d, minimap2 2.31-r1302, samtools
1.25, `--min-required-hits` 3, `--strictness` 0.7. 146 typing runs. The base counts of a canSNP
do not depend on `--min-depth` or `--min-allele-fraction`, which are applied to the counts
afterwards, so the `<sample>_site_counts.tsv` of each run was classified again with
`ParsePileup.classify_counts` and walked again with `NewickTree.draw_tree` at six values of
`min_depth` and five of `min_allele_fraction`: 4,380 results. The harness reproduces the calls,
the support and the site counts that CanSNPer2 recorded for the eight full-depth runs on
2026-11-25, at those runs' own settings.

**Date.** 2026-11-27.

**Results, base quality.** With `--platform ont`, whose `min_baseq` is 10, the 50x file of one
isolate gave a depth of zero at all 42 canSNPs and the call `NA`: a quality line of repeated
`!` leaves no base. `min_baseq` 0 is a requirement of these files, not a preference.

**Results, minimap2 preset.** At 30x, `min_depth` 10 and `min_allele_fraction` 0.7, medians
over the eight isolates:

| preset | canSNPs with reads of 42 | depth | the assembly call | canSNPs with reads of 283 | depth | the assembly call |
|---|---|---|---|---|---|---|
| `map-pb` (`clr`) | 42 | 26 | 8 of 8 | 276 | 20 | 5 of 8, 3 a node above |
| `map-ont` | 42 | 26 | 8 of 8 | 277 | 22 | 5 of 8, 3 a node above |
| `map-hifi` | 42 | 22 | 8 of 8 | 76 | 8 | 0 of 8, 8 NA |
| `sr` | 42 | 5 | 0 of 8, 8 NA | 96 | 3 | 0 of 8, 8 NA |

Against a reference genome of the same species every long-read preset reaches every canSNP.
Against the reference genome of another species `map-hifi` keeps 76 of the 283 positions and
`sr` 96, while `map-pb` and `map-ont` keep 276 and 277. `map-pb` is confirmed, `map-ont` is
equivalent on these reads, and `map-hifi` is not a substitute for CLR reads.

**Results, depth.** At `min_depth` 10 and `min_allele_fraction` 0.7, the values the preset had
while it was being measured (`min_depth` is 5 since this measurement).
"a node above" is a node on the path from the root to the call of the assembly:

| depth of the reads | depth at the canSNPs, the two schemes | brucella_anthropi | brucella_genus |
|---|---|---|---|
| 5x | 4, 3 | 8 NA | 8 NA |
| 10x | 8, 7 | 8 NA | 2 a node above, 6 NA |
| 20x | 18, 14 | 5 the assembly call, 3 above | 3 the assembly call, 5 above |
| 30x | 26, 20 | 8 of 8 | 5 of 8, 3 above |
| 50x | 43, 35 | 8 of 8 | 6 of 8, 2 above |
| full | 523, 414 | 8 of 8 | 8 of 8 |

With `min_depth` 5 the same calls arrive from 20x, and at 10x four of the eight isolates reach
the assembly call with the dedicated scheme and five with the genus scheme; with `min_depth` 3,
six and eight of eight at 10x. No value of `min_depth` gave a call on another branch at
`min_allele_fraction` 0.7.

**Results, allele fraction.** States compared with the assembly and states differing from it,
over every depth of the `clr` arm and both schemes, and the calls at 30x of the 16 samples
(eight isolates, two schemes) at `min_depth` 10:

| `--min-allele-fraction` | states compared | states differing | the assembly call at 30x |
|---|---|---|---|
| 0.5 | 6,876 | 15 | 16 of 16 |
| 0.6 | 6,833 | 4 | 16 of 16 |
| 0.7 (`clr`) | 6,558 | 1 | 13 of 16 |
| 0.8 | 5,704 | 1 | 8 of 16 |
| 0.9 | 3,744 | 1 | 1 of 16 |

Every state that differs from the assembly comes from 5x or 10x: 18 of 2,349 comparisons at 5x
and 4 of 3,548 at 10x, and none of the 23,818 comparisons from 20x upward. One of the 4,380
results was a call on another branch: at 5x with `min_depth` 5 and `min_allele_fraction` 0.5,
one isolate was called `Bant.anthropi.T`, the branch beside `Bant.anthropi.env`, on 3 of 3
canSNPs of that path while 21 of the 42 positions were below the depth threshold.

**Results, mixed sites.** Against the reference genome of the same species no sample had a
mixed site at 20x or above. Against 16M the median rose with depth, from 4 at 20x to 6 at 30x
and 7 at 50x and at full depth, with 11 in one sample. Mixed canSNPs on the path are why
several samples stop a node above the assembly call between 20x and 50x; at full depth all 16
canSNPs of the path are derived in all eight samples. Single-pass long reads mapped to a
reference genome of another species therefore give mixed sites that are not a mixture of
strains.

**Conclusion.** `min_baseq` 0 is required, `map-pb` is the right preset, and
`min_allele_fraction` 0.7 lies in the band that gave no call on another branch (0.6 and 0.7);
0.5 gave the one off-path call and most of the differing states, 0.8 and above lose calls that
are correct. `min_depth` was 10 while the preset was being measured, which asks for about 30x.
The measurement moved it to 5, the value of the other platforms: 5 gives the same calls from 20x
and adds calls at 10x, and no value of `min_depth` gave a call on another branch at
`min_allele_fraction` 0.7. For a run of CLR reads this means about 20x, and a sample below that
should be expected to stop above its lineage or to give `NA`.

**Open.** One species, eight isolates of one experiment, two draft schemes, one instrument and
two reference distances. The depth series is drawn from one deep run per isolate, so it does not
show the uneven coverage of a run that is shallow to begin with.

## Read typing on real nanopore reads

**Compared.** The final SNP and the canSNP states from real nanopore reads against the call of
an assembly of the same isolate and against the call of its Illumina reads, on five schemes, at
several depths, with two nanopore chemistries, and grouped by the homopolymer run around each
canSNP position. The evaluation on simulated reads above used reads whose errors are close to
random; this one uses runs from instruments.

**Data.** 17 runs of seven isolates: *F. tularensis* holarctica FDC408 (PRJNA287971),
*C. burnetii* Nine Mile II RSA439 (PRJNA541584), *Y. pestis* YP_80 (PRJNA1098581), *B.
anthracis* 14RA5914 and 14RA5916, and *B. suis* 08RB3277 and 08RB3701 (PRJEB59317 and
PRJEB94163, an evaluation of nanopore sequencing by one laboratory). Eleven nanopore runs, six
Illumina runs of the same BioSamples. Four isolates were sequenced twice by that laboratory,
with R9 and with R10.4.1 chemistry, which is where the chemistry comparison comes from. Five
isolates have a complete or draft assembly in RefSeq (GCF_002102455.2, GCF_005280755.1,
GCF_039640425.1, GCF_002277915.1, GCF_017963425.1); the two *B. suis* isolates have none.

A nested depth series of 5x, 10x, 20x, 30x and 50x was drawn from each nanopore run, and the
whole run was typed as well. The Illumina pairs were typed with every read mapped. Both
comparisons are reported because neither is beyond question on its own: three of the complete
genomes are hybrid or nanopore-polished assemblies of the same reads, so they are not
independent of them, while the Illumina reads are independent but are reads as well.

**Software and settings.** CanSNPer2 2.6.0 at commit 79ee2e8, minimap2 2.31-r1302, samtools
1.25, sra-tools 3.4.1, Python 3.13.18. The platform defaults of 2.6.0, that is `--platform ont`
with `--min-depth` 5 and `--min-allele-fraction` 0.7, except in the threshold table. The five
schemes are from CanSNPer2-data master. The thresholds were varied afterwards by classifying
the kept site counts again and walking the tree again, as in the measurement of the `clr`
preset; the harness reproduced all 77 calls of the runs at the shipped thresholds.

**Date.** 2026-11-27.

**Results.** With the shipped settings no call was on another branch than the call of the
assembly, at any depth. The call either matched, stopped higher on the same path, or was `NA`,
over the seven samples that have an assembly to compare with:

| | 5x | 10x | 20x | 30x | 50x | full |
|---|---|---|---|---|---|---|
| same call as the assembly | 0 of 7 | 3 of 7 | 6 of 7 | 7 of 7 | 7 of 7 | 7 of 7 |
| stopped higher on the same path | | 2 of 7 | 1 of 7 | | | |
| `NA` | 7 of 7 | 2 of 7 | | | | |

The depth at which the assembly call first appeared differs by scheme: 10x for *C. burnetii*
(11 canSNPs) and for the two *B. anthracis* isolates on R9, 20x for *Y. pestis* (97 canSNPs) and
for the same two isolates on R10.4.1, and 30x for *F. tularensis* (328 canSNPs). At 10x
*F. tularensis* stopped at `B.3` and *Y. pestis* at `VI`, both ancestors of the assembly call
(`B.75` and `VII`). Simulated ONT reads reached the assembly call at 10x for 42 of 50 genomes;
real runs need about 20x, and a deep scheme about 30x.

The depth of the series is the depth over the genome, and the depth at the canSNP positions is
lower: the median depth there was 2 to 5 at 5x, 5 to 9 at 10x, 9 to 19 at 20x, 13 to 27 at 30x
and 23 to 45 at 50x, so about two thirds of the nominal depth after the mapping quality filter.
This is why a 10x run sits at the `--min-depth` of 5 rather than above it.

**canSNP states.** At the shipped thresholds no canSNP state from the nanopore reads differed
from the comparison: 2,752 states against the assemblies and 2,754 against the Illumina reads,
over all depths and both chemistries, none differing. One canSNP of one run was flagged as a
mixed site: *F. tularensis* FDC408 at 30x, where `A.I.5` had 4 reads of the ancestral base and 2
of a third base, so the most frequent base reached 0.67 and stayed under the allele fraction of
0.7. The call of that run was the assembly call. The measurement of the `clr` preset found 4 to
7 mixed sites per sample against a reference genome of another species; here every scheme
reference is of the species of the sample, which is the difference.

**Thresholds.** The two thresholds were varied over a grid, 3,577 results with 131,905 state
comparisons. All 39 differing states are at `--min-depth` 1, and so is the only call on another
branch. At `--min-depth` 3 and 5 no state differed and no call was off the path at any allele
fraction from 0.5 to 0.9. Calls equal to the assembly call, at allele fraction 0.7:

| `--min-depth` | 5x | 10x | 20x | 30x | 50x | full | off the path |
|---|---|---|---|---|---|---|---|
| 1 | 7 of 7 | 6 of 7 | 7 of 7 | 7 of 7 | 7 of 7 | 7 of 7 | 1, at fraction 0.5 |
| 3 | 1 of 7 | 6 of 7 | 7 of 7 | 7 of 7 | 7 of 7 | 7 of 7 | none |
| 5 (shipped) | 0 of 7 | 3 of 7 | 6 of 7 | 7 of 7 | 7 of 7 | 7 of 7 | none |
| 10 | 0 of 7 | 0 of 7 | 5 of 7 | 7 of 7 | 7 of 7 | 7 of 7 | none |
| 20 | 0 of 7 | 0 of 7 | 0 of 7 | 1 of 7 | 6 of 7 | 7 of 7 | none |

The one off-path call is *F. tularensis* at 5x with `--min-depth` 1 and `--min-allele-fraction`
0.5, which called `B.21` instead of `B.75` after reading one canSNP of the path (`B.66`) wrong.
The same reads give `NA` at the shipped thresholds.

**Homopolymer context.** The canSNP positions were grouped by the longest run of identical bases
that reaches the position or stops within one base of it, measured in the reference genome of
the scheme. No state differed from the comparison in any group, at any depth, in any of the
eleven nanopore runs. The groups are uneven, which is set by the schemes: of the 473 canSNPs of
the five schemes, 376 sit in a run shorter than three, 82 in a run of three or four and 15 in a
run of five or more, and all 15 are in the *F. tularensis* and *Y. pestis* schemes. At full
depth those 15 positions were read correctly in both samples.

**Chemistry.** The four isolates sequenced with both chemistries gave the same call from both,
and no state differed from the Illumina comparison with either. From 20x both chemistries had
usable data at the same number of canSNPs (13 of 13 for *B. anthracis*, 23 of 24 for *B. suis*).
Below that the counts differ in both directions: at 10x R9 covered more canSNPs in the two
*B. anthracis* isolates (12 and 11 against 9 and 10) and R10.4.1 more in the two *B. suis*
isolates (23 and 22 against 20 and 20). With four isolates this is not a difference worth
reading as one. The two *B. anthracis* isolates reached the assembly call at 10x on R9 and at
20x on R10.4.1, which follows the same counts.

**Platform detection.** `--platform auto` read the platform from the lengths and the base
qualities of every one of the eleven nanopore files, and reported `ont` for seven of them and
`hifi` for four. The four are the R10.4.1 runs of the one laboratory, basecalled with a `sup`
model, whose mean base quality over the first 400 reads is 25.8 to 37.3, above the 25 at which
long reads are taken for HiFi reads. The module description gives this as a known limit of the
lengths and the qualities. The calls at 30x were the same as with `--platform ont` for all
eleven files, since `map-hifi` and `map-ont` both map these reads and the deeper `hifi`
thresholds were still met at that depth.

The read names of these files still carry the nanopore read UUID, and in one case the
basecalling model, but the SRA conversion puts the run accession first and the UUID form of
`ReadPlatform.NAME_FORMS` was anchored at the start of the name, so the name signal was not used
for any file while this was measured. An unanchored form matches 400 of 400 names in each of the
ten files that still have a UUID and none of the 400 names of either Illumina file, which is why
the forms are searched from the start of the name and after each space since this measurement
(#124). Nine of the eleven nanopore files are now read from their names, the four reported as
`hifi` above among them, and they are reported as `ont`; the two whose names the SRA replaced
altogether are still read from their lengths and qualities. The base quality that separates `ont`
from `hifi` is left at 25 all the same: these four runs show that no value separates the
platforms, so the names are the signal to rely on where they survive.

**A scheme that does not cover the sample.** Neither *B. suis* isolate was called, from any
input: 22 of the 23 canSNPs with data are ancestral and only the root is derived, from the
Illumina reads at 129x as much as from the nanopore runs. The `brucella_classical` scheme has no
node for these isolates, so `not called` is the right report, and the two isolates contribute
only their canSNP states to this evaluation. See
[limitations.md](limitations.md) on taxa outside the scheme.

**Conclusion.** On real nanopore runs of four schemes the shipped `ont` settings gave the call
of the assembly and of the Illumina reads with no call on another branch, and no canSNP state
differed. About 20x over the genome is the depth to plan for, 30x for a scheme with many
canSNPs, and a run below that should be expected to stop above its lineage or to give `NA`.
Homopolymer context did not change a single state, which is the error the simulated reads could
not cover, although only 15 canSNP positions of these schemes sit in a run of five or more
bases. The errors that the thresholds hold back appear at `--min-depth` 1; 5 and 3 were both
without error here, and 3 would add the call at 10x for 6 of 7 samples instead of 3. The value
stays 5: 3 is one step from the value where every differing state and the one off-path call of
the grid are, and what 5 does at 10x is stop a node above on the same path, which is the safe
way to be wrong on thin data. The same measurement moved `clr` from 10 to 5, where the gain is
at 20x rather than at the thinnest depths and all four platforms then ask for the same depth.

**Open.** Seven isolates, five schemes, two laboratories for the nanopore runs and one
instrument type. The depth series is drawn from one deep run per isolate, so it does not show
the uneven coverage of a run that is shallow to begin with. R10.4.1 duplex reads and barcoded
runs were not separated out. Twelve of the runs were taken from the NCBI Open Data bucket and converted with
`fasterq-dump`, which drops the nanopore tags that follow the read name; the read UUID itself
survives in nine of the eleven files and is what the names are read from, so a file served with
its full nanopore tags is not covered by this data. No PacBio HiFi run is included, which remains
open together with the aligner comparison.

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
