# Limitations

This page states what CanSNPer2 is designed for and the conditions under which its results are
known to be incomplete or misleading. It is for users who decide whether CanSNPer2 fits a task
and for those who report its results.

## Intended use

CanSNPer2 assigns a bacterial isolate to a lineage of a canSNP scheme. The input is one
assembly or one read set of a single isolate of an organism that the scheme covers, and the
result is the deepest node of the scheme tree that the canSNPs of the isolate support. The
result is only as specific and as correct as the scheme: its tree, the choice of canSNPs and
their alleles.

## Out of scope

- Species identification. CanSNPer2 compares a sample only with the reference genomes of the
  scheme it is given and does not check that the sample belongs to the scheme's taxon. A genome
  from another taxon can be reported as `NA`, or as a node of the scheme (see below). Identify
  the species with another method first.
- Taxa outside the scheme. A lineage that the scheme has no node for is placed at the deepest
  node it shares with the scheme, or not called.
- Metagenomes and other samples with several organisms. Reads of other organisms that map to a
  reference genome are counted at the canSNP positions like reads of the target.
- Strain mixtures. For reads, CanSNPer2 flags a sample with several mixed canSNP sites
  (`--mixed-sites`) as possibly mixed; it does not resolve the strains or their proportions. For
  assemblies there is no mixture check: an assembly holds one base per position.

## Known failure modes

### Genomes from outside the scheme's taxon

progressiveMauve (`--aligner progressiveMauve`; the default aligner up to 2.4.0) aligns
genomes of other genera to the reference genomes and can then produce a call. In an evaluation on genomes of other Francisellaceae genera, progressiveMauve gave a
Francisella tularensis final SNP for 5 of the 77 genomes it was run on (all 5 of other genera),
while minimap2 and MUMmer4 gave no call for any of 268 genomes of other Francisellaceae species
and genera (see [validation.md](validation.md)). minimap2 with the `asm20`
preset, the default aligner, does not align genomes this distant. The final SNP of a genome whose species has not
been confirmed should not be reported.

### Reads far from the reference genomes

Reads map to a reference genome only where they are close to it, while an assembly aligns over
larger distances. In a scheme that spans a genus with one reference genome, canSNPs of lineages
far from that reference get no reads, and the call stops above them or is `NA`. In the
evaluation of read typing, the Brucella genomes outside the classical cluster were the ones
whose 50x Illumina reads missed the assembly call. `--local-references` maps the reads again
to a representative genome of the scheme closer to the sample and reduces this effect (see
[method.md](method.md#local-references) and [validation.md](validation.md)); it needs the
genomes file of the scheme and the representative genomes.

### Duplicated and paralogous regions

A canSNP in a region that occurs more than once in a genome can be read from the wrong copy.
CanSNPer2 limits this in several ways: minimap2 alignments are read only as primary alignments
with a mapping quality of at least 20, MUMmer4 alignments are filtered with `delta-filter -1`,
and a position covered by two alignments with different bases is reported as no data. For
reads, a region that the sample has twice while the reference genome has it once gives mixed
sites, which can flag an isolate as possibly mixed; in the read typing evaluation this happened
for 1 of 50 genomes. None of these measures detects a copy that the reference genome lacks and
that aligns better than the true locus.

### Low read depth and low base quality

A canSNP with fewer than `--min-depth` bases (default 5) after the mapping and base quality
filters has no data. At low depth, the call stops higher in the tree or is `NA`. In the read
typing evaluation on simulated reads, 5x Illumina coverage gave the assembly call for none of
50 genomes. The ONT defaults were checked on simulated reads whose errors are close to random;
systematic errors of real ONT reads, as at homopolymers, were not part of that evaluation.

Only the first 300000 reads (read pairs) of a sample are mapped by default (`--max-reads`). For
Illumina reads this limits the depth to about 45x on a 2 Mb genome with 2 x 150 base reads, and
less on a larger genome. A read file sorted by position is covered only in part by its first
reads; type such a file with `--max-reads 0`.

### Fragmented assemblies

canSNPs in gaps between contigs have no data, and canSNPs within 25 bases of a contig or block
end are flagged `edge`. A draft assembly can therefore give a shallower call, or a `review`
status, compared with a complete genome of the same isolate.

### Draft schemes

Draft schemes (on the `dev` branch of CanSNPer2-data, fetched with `--branch dev`) are under
evaluation and have not been verified. Their files are checked against the manifest of the
repository when it has one, as for verified schemes; this shows that the files are those
published, not that the scheme gives correct calls. CanSNPer2 warns about a draft scheme, a
scheme of unknown status and a database modified after it was built (see
[schemes-and-databases.md](schemes-and-databases.md#verified-and-draft-schemes)). Results from
such schemes are provisional.

### Differences between aligners

The three aligners (minimap2, the default; progressiveMauve; MUMmer4) do not cover the same
canSNP positions. In the comparison on 937
F. tularensis genomes, the final calls were the same with all three, and no position had a
different base, but each aligner had data at some positions where another had none (see
[validation.md](validation.md)). Near 1 bp indels in alignments of low identity, the base read
at a canSNP can depend on where the aligner places the gap. Results made with different
aligners can therefore differ in support and, near the limits of the rules, in the final SNP;
record the aligner with each result.

### Scheme limits

The tree walk does not choose between two paths confirmed at the same depth; such a sample is
`NA`. Untyped (alternative) canSNPs of a node are stored in the database but are not read for
typing.
