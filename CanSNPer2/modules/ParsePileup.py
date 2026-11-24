'''
ParsePileup reads the output of samtools mpileup at the canSNP positions of one reference,
made from sequencing reads mapped to that reference (read typing, see Aligners.ReadMapper).

An assembly gives one base at a position; reads give a count of each base. The base of a
canSNP is decided here from those counts (classify_counts):

	no reads                                   no data
	fewer than min_depth A/C/G/T bases         no data (low depth)
	derived base >= min_allele_fraction        derived
	ancestral base >= min_allele_fraction      ancestral
	another base >= min_allele_fraction        other base
	no base reaches min_allele_fraction        no data (mixed)

A mixed site gives no evidence either way to the tree walk, as a position without data does;
it is counted, because several of them in one sample suggest that the sample holds more than
one strain. mpileup has already removed reads below the mapping quality and bases below the
base quality given to it, so depth here is the number of bases that passed those filters.

Besides the SNP states and the rows of the result files, which have the columns of the
assembly typing (QuerySeq is "reads", QueryPosition the depth), the reader writes the counts
of every canSNP of the reference to <pileup>.counts (SITE_COUNT_COLUMNS), from which the
per sample site count table is made.
'''

import logging

from CanSNPer2.modules.AlignmentReader import AlignmentReader, NO_DATA, DERIVED, ANCESTRAL, OTHER_BASE

logger = logging.getLogger(__name__)

BASES = "ACGT"

'''the settings of each sequencing platform: the minimap2 preset and the thresholds of the
	pileup. Checked on reads simulated from 30 Brucella and 20 F. tularensis genomes against
	their assembly calls (2026-11-24): with these values no read call was on another branch
	than the assembly call at any depth, Illumina reached the assembly call from 10x and ONT
	from 10x (42 of 50; min_depth 10 gave no call for 42 of 50 at 10x). The ONT reads were
	simulated with badread, whose errors are close to random; real ONT errors at homopolymers
	are not covered'''
PLATFORMS = {
	"illumina": {"preset": "sr", "min_mapq": 20, "min_baseq": 20, "min_depth": 5, "min_allele_fraction": 0.8},
	"ont": {"preset": "map-ont", "min_mapq": 20, "min_baseq": 10, "min_depth": 5, "min_allele_fraction": 0.7},
	"hifi": {"preset": "map-hifi", "min_mapq": 20, "min_baseq": 20, "min_depth": 5, "min_allele_fraction": 0.8},
}
DEFAULT_PLATFORM = "illumina"

SITE_COUNT_COLUMNS = ["Name", "Reference", "Sequence", "Pos", "Ancestral base", "Derived base",
	"A", "C", "G", "T", "depth", "fraction_derived", "state", "source", "lifted_position"]

STATES = {"derived": DERIVED, "ancestral": ANCESTRAL, "other_base": OTHER_BASE,
	"mixed": NO_DATA, "low_depth": NO_DATA, "no_reads": NO_DATA, "not_lifted": NO_DATA, "conflict": NO_DATA}
COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A"}


def parse_bases(column, ref="N"):
	'''the A, C, G and T in the base column of an mpileup line, both strands together.
		"^" starts a read and is followed by its mapping quality character, "$" ends one,
		"+n"/"-n" followed by n bases is an insertion or deletion after the base, "*" and "#"
		are deleted bases, ">" and "<" skipped reference (spliced) bases; "." and "," are the
		reference base, counted only when the reference base is known'''
	counts = dict.fromkeys(BASES, 0)
	ref = ref.upper()
	i, n = 0, len(column)
	while i < n:
		ch = column[i]
		if ch == "^":
			i += 2
			continue
		if ch in "+-":
			j = i + 1
			while j < n and column[j].isdigit():
				j += 1
			i = j + int(column[i+1:j] or 0)
			continue
		if ch in ".,":
			if ref in counts:
				counts[ref] += 1
		else:
			up = ch.upper()
			if up in counts:
				counts[up] += 1
		i += 1
	return counts


def classify_counts(counts, ancestral, derived, min_depth, min_allele_fraction):
	'''(label, state, base, depth, fraction of the derived base) of a canSNP from its base
		counts; base is the base the site is called as, N when it is not decided by the
		fraction, "-" without reads'''
	depth = sum(counts.get(b, 0) for b in BASES)
	if depth == 0:
		return "no_reads", NO_DATA, "-", 0, 0.0
	fraction = counts.get(derived, 0) / depth
	if depth < min_depth:
		return "low_depth", NO_DATA, "N", depth, fraction
	top = max(BASES, key=lambda b: counts.get(b, 0))
	if counts.get(top, 0) / depth < min_allele_fraction:
		return "mixed", NO_DATA, "N", depth, fraction
	if top == derived:
		return "derived", DERIVED, top, depth, fraction
	if top == ancestral:
		return "ancestral", ANCESTRAL, top, depth, fraction
	return "other_base", OTHER_BASE, top, depth, fraction


def read_site_counts(path):
	'''the rows of a site count file as dicts by SITE_COUNT_COLUMNS'''
	with open(path) as fin:
		header = fin.readline().rstrip("\n").split("\t")
		return [dict(zip(header, line.rstrip("\n").split("\t"))) for line in fin if line.strip()]


class PileupReader(AlignmentReader):
	'''Reader of samtools mpileup output (one sample, without -f or with it); the reference
		sequence and the 1-based position of each line are the database coordinates'''
	def __init__(self, verbose=False, min_depth=5, min_allele_fraction=0.8, lift=None, source=None, **kwargs):
		'''lift: {snp_id: LiftOver.Lifted} when the reads were mapped to a local reference: the
			pileup is then on the contigs of that genome and is read at the lifted positions, the
			counts of a canSNP lifted to the "-" strand complemented. source: the genome the reads
			were mapped to, for the site count rows (default: the reference)'''
		super(PileupReader, self).__init__(verbose=verbose, **kwargs)
		self.min_depth = int(min_depth)
		self.min_allele_fraction = float(min_allele_fraction)
		self.lift = lift
		self.source = source

	def _reset(self):
		super(PileupReader, self)._reset()
		self.mixed_sites = []
		self.counts_path = None

	def read_alignment(self, path):
		'''the base counts of the canSNP positions, in observations {(sequence, position): counts}'''
		self.counts_path = getattr(self, "counts_file", None) or path + ".counts"
		lifted = {}
		if self.lift is not None:
			by_id = {snp_id: key for key, (_a, _d, snp_id) in self.snplist.items()}
			lifted = {(lf.contig, lf.position): (by_id[snp_id], lf.strand) for snp_id, lf in self.lift.items() if snp_id in by_id}
		with open(path) as fin:
			for number, line in enumerate(fin, 1):
				if not line.strip():
					continue
				fields = line.rstrip("\n").split("\t")
				if len(fields) < 5 and not (len(fields) == 4 and fields[3] == "0"):
					raise ValueError("{path} line {n} is not samtools mpileup output (sequence, position, reference base, depth, bases, qualities): {line}".format(
							path=path, n=number, line=line.strip()[:80]))
				sequence = fields[0]
				if self.lift is not None:
					site = lifted.get((sequence, int(fields[1])))
					if site is None:
						continue
					key, strand = site
					counts = parse_bases(fields[4] if len(fields) > 4 else "", fields[2])
					if strand == "-":
						counts = {COMPLEMENT[b]: n for b, n in counts.items()}
					self.observations[key] = counts
					continue
				if sequence not in self.snp_positions:
					self.unknown_sequences.add(sequence)
					continue
				key = (sequence, int(fields[1]))
				if key not in self.snplist:
					continue
				self.observations[key] = parse_bases(fields[4] if len(fields) > 4 else "", fields[2])

	def finish(self):
		'''the states, rows and called SNPs as AlignmentReader.finish gives them for an
			assembly, from the base counts; every canSNP of the reference goes into the site
			count file, also those the pileup has no line for'''
		site_rows = []
		source = self.source or self.reference
		for key in sorted(self.snplist):
			sequence, position = key
			ancestral, derived, snp_id = self.snplist[key]
			counts = self.observations.get(key, dict.fromkeys(BASES, 0))
			label, state, base, depth, fraction = classify_counts(counts, ancestral, derived, self.min_depth, self.min_allele_fraction)
			lifted = ""
			if self.lift is not None:
				lf = self.lift.get(snp_id)
				if lf is None:
					label = "not_lifted"
				else:
					lifted = "{c}:{p}{s}".format(c=lf.contig, p=lf.position, s=lf.strand)
			site_rows.append([snp_id, self.reference, sequence, str(position), ancestral, derived] +
					[str(counts[b]) for b in BASES] + [str(depth), "{0:.2f}".format(fraction), label, source, lifted])
			if label == "mixed":
				self.mixed_sites.append(snp_id)
			if depth == 0:
				continue  ## no reads: listed by the tree as a position without data, as for assemblies
			self.SNPS[snp_id] = state
			self.SNP_info.append([snp_id, self.reference, sequence, str(position), "reads", str(depth), ancestral, derived, base, ""])
			if state == DERIVED:
				self.called_snps.append(snp_id)
		if self.counts_path:
			with open(self.counts_path, "w") as out:
				print("\t".join(SITE_COUNT_COLUMNS), file=out)
				for row in site_rows:
					print("\t".join(row), file=out)
		return self.SNPS
