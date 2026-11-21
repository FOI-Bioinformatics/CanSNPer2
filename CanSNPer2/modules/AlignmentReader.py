'''
AlignmentReader holds what every alignment format shares: the SNP positions of one
reference from the database, and the classification of the query base found at each of
them (derived, ancestral, other base or no data), the edge flag, the --mask filter and
the rows of the per sample SNP tables.

A reader for an alignment format only has to walk its file and report, for each aligned
stretch of the reference, the query base at the SNP positions inside it through observe().
A SNP is identified by the reference sequence it lies on (the first word of its fasta
header, a GenBank accession for the published schemes) and its 1-based position on that
sequence, as the database stores it. minimap2 and MUMmer4 report positions per sequence;
progressiveMauve counts along the sequences joined end to end, and its reader converts them
with the sequence table written next to the alignment (write_sequence_table).
'''

import os
import logging
from bisect import bisect_left, bisect_right

from CanSNPer2.modules.DatabaseConnection import XMFAFunctions

logger = logging.getLogger(__name__)

'''file ending of the reference sequence table written next to an alignment: name and length
	of each reference sequence, in file order, for readers whose alignment format counts
	positions along the joined reference'''
SEQUENCES_SUFFIX = ".sequences"


def write_sequence_table(fasta, path):
	'''write name and length of the sequences of a fasta file, in file order'''
	rows = []
	with open(fasta) as fin:
		for line in fin:
			if line.startswith(">"):
				rows.append([line[1:].split()[0] if line[1:].split() else "", 0])
			elif rows:
				rows[-1][1] += len("".join(line.split()))
	with open(path, "w") as out:
		for name, length in rows:
			out.write("{0}\t{1}\n".format(name, length))


def read_sequence_table(path):
	'''[(name, length)] of the reference sequences, in file order, from a table written by
		write_sequence_table'''
	rows = []
	try:
		with open(path) as fin:
			for line in fin:
				if not line.strip():
					continue
				name, length = line.rstrip("\n").split("\t")
				rows.append((name, int(length)))
	except FileNotFoundError as e:
		raise FileNotFoundError("The reference sequence table {path} of the alignment is missing. It is written when the alignment is made; alignments made before CanSNPer2 2.3.0 have none, align again without --skip_alignment".format(path=path)) from e
	return rows

'''SNPs closer to either edge of an alignment block than this many reference bases are
	flagged in the output: assembly ends and read ends are the most likely place for
	alignment errors, and a SNP sitting there should be reviewed by the user instead of
	silently trusted or silently dropped (--mask, issue #48)'''
EDGE_MARGIN = 25

'''SNP states, as NewickTree reads them'''
NO_DATA = 0
DERIVED = 1
ANCESTRAL = 2
OTHER_BASE = 3


class AlignmentReader(object):
	'''Base class of the alignment readers. Subclasses implement read_alignment(path).'''
	def __init__(self, verbose=False, **kwargs):
		super(AlignmentReader, self).__init__()
		self.verbose = verbose
		'''Mask option, for distant relatives false SNPs near edges of alignments may become a problem,
			this option allows masking of SNPs placed within n bases of the edge of an alignment'''
		self.mask = int(kwargs.get("mask",0) or 0) ## the option was advertised but never read (issue #175)
		self.masked = {} ## how many positions each reference lost to the mask
		'''the database is only connected in run(), by whichever process uses this object:
			the SNP calling workers are forked from the parent and used to inherit its open
			sqlite connection and cursor, so several processes read through one cursor (issue #49)'''
		self.database_path = kwargs.get("database",False)
		self.database = False
		self.reference = None
		'''the query genome the alignments were made from, for readers that need its sequence
			names and lengths to report query positions (set per query by CanSNPer2)'''
		self.query_file = kwargs.get("query",None)
		self._reset()

	def _reset(self):
		'''state of one parsed alignment'''
		self.snplist = {}          ## {(sequence, position): (ancestral base, derived base, snp id)}
		self.snp_positions = {}    ## {sequence: sorted positions of its SNPs}
		self.observations = {}     ## {(sequence, position): observation row}, see observe()
		self.SNPS = {}
		self.SNP_info = []
		self.called_snps = []
		self.edge_snps = set() ## SNPs within EDGE_MARGIN bases of an alignment edge, flagged not dropped (issue #48)

	def get_snp_info(self):
		'''Return SNP_info'''
		return self.SNP_info

	def get_snps(self):
		'''Return snps'''
		return self.SNPS

	def get_called_snps(self):
		'''Return true snps'''
		return self.called_snps

	def snps_in(self, sequence, start, end):
		'''the registered SNP positions on reference sequence from start to end, both included,
			in ascending order; [] for a sequence without SNPs (a plasmid the scheme does not
			use, or a sequence not in the database). Looking positions up per block makes the
			result independent of the order in which the blocks come in the file: walking the
			sorted SNP list as a queue dropped every SNP that sat in a block coming after a
			block further along the reference'''
		positions = self.snp_positions.get(sequence)
		if not positions:
			return []
		return positions[bisect_left(positions,start):bisect_right(positions,end)]

	def observe(self, sequence, position, base, block_offset, block_length, query_name="-", query_pos="N/A"):
		'''Record the query base found at a SNP position of the reference.

			sequence      reference sequence of the SNP
			position      1-based position of the SNP on that sequence ((sequence, position)
			              is a key of snplist)
			base          query base on the forward strand of the reference, "-" for a gap in
			              the query; any other character than A, C, G, T, "-" or N is kept as
			              it is and reads as "other base"
			block_offset  1-based position of the SNP in the aligned stretch of the reference
			block_length  length of that stretch in reference bases; with the offset it gives
			              the distance to both edges, for --mask and the edge flag
			query_name    query sequence (contig) the base comes from
			query_pos     position of the base in that query sequence, "N/A" for a gap
		'''
		if self.mask and (block_offset <= self.mask or block_length - block_offset < self.mask):
			'''the flanking mask: a SNP this close to either edge of the alignment is not
				trustworthy, and it is dropped completely: counting it as absent would penalise
				the sample for a base that was never compared (issue #175)'''
			self.masked[self.reference] = self.masked.get(self.reference,0)+1
			return
		edge = block_offset <= EDGE_MARGIN or block_length - block_offset < EDGE_MARGIN
		base = base.upper()
		key = (sequence, position)
		previous = self.observations.get(key)
		if previous is not None:
			'''a position covered by more than one alignment block: progressiveMauve never
				produces that for the reference, other aligners can (repeats, duplicated
				contigs). The same base twice is one observation; different bases are no
				evidence either way, so the position is reported as an N'''
			if previous["base"] != base:
				logger.warning("SNP {snp} at position {pos} of {sequence} in {reference} is covered by more than one alignment with different query bases ({a}, {b}), it is treated as N".format(
						snp=self.snplist[key][2], pos=position, sequence=sequence, reference=self.reference, a=previous["base"], b=base))
				previous["base"] = "N"
			previous["edge"] = previous["edge"] or edge
			return
		self.observations[key] = {"base": base, "edge": edge, "query_name": query_name or "-", "query_pos": str(query_pos)}

	def classify(self, base, ancestral, derived):
		'''state of a SNP from the query base'''
		if base == derived:
			return DERIVED
		if base == ancestral:
			return ANCESTRAL
		if base != "-" and base != "N":
			return OTHER_BASE
		return NO_DATA

	def finish(self):
		'''turn the observations into the SNP states, the result rows and the called list, in
			(sequence, position) order. A row is [snp_id, reference, sequence, position,
			query sequence, query position, ancestral base, derived base, query base, edge flag],
			the columns of CanSNPer2.SNP_FILE_COLUMNS'''
		for key in sorted(self.observations):
			obs = self.observations[key]
			sequence, position = key
			rbase,tbase,snp_id = self.snplist[key]
			state = self.classify(obs["base"], rbase, tbase)
			self.SNPS[snp_id] = state
			edge_flag = "edge" if obs["edge"] else ""
			if edge_flag:
				self.edge_snps.add(snp_id)
			self.SNP_info.append([snp_id,self.reference,sequence,str(position),obs["query_name"],obs["query_pos"],rbase,tbase,obs["base"],edge_flag])
			if state == DERIVED:
				self.called_snps.append(snp_id)
		return self.SNPS

	def read_alignment(self, path):
		'''walk the alignment file and call observe() for the SNP positions it covers'''
		raise NotImplementedError("{cls} does not implement read_alignment".format(cls=type(self).__name__))

	def run(self, alignment, reference=False, database=False):
		'''Read one alignment against one reference and return the SNP states {snp_id: state}.
			The reference is the genome name in the database; it is taken from the file name
			(<reference>_<query>.<suffix>) only when the caller does not give it, which fails
			for reference names that contain "_".'''
		self._reset()
		'''Create the connection to the SNP database here, in the process that parses the
			alignment, so a forked worker never reads through a cursor opened by its parent
			(issue #49)'''
		if not self.database:
			database = self.database_path or database
			if not database:
				raise ValueError("No database was provided to {cls}, cannot retrieve snps from the database!".format(cls=type(self).__name__))
			self.database = XMFAFunctions(database,verbose=self.verbose)
		if not reference:
			reference = os.path.basename(alignment).split("_")[0]
		self.reference = reference
		self.snplist, self.snp_positions = self.database.get_snps(reference)
		if not self.snp_positions: ## the reference has no registered SNPs in this database, nothing to look for (issue #99)
			logger.warning("No SNPs are registered in the database for reference {reference}, no SNP calls possible".format(reference=reference))
			return self.SNPS
		try:
			self.read_alignment(alignment)
		except FileNotFoundError as e:
			'''this used to return False, which every caller read as "no SNPs in this
				alignment", so a missing alignment typed the sample NA and exited 0
				(issue #176, #179). A file the reader needs besides the alignment (the
				sequence table of progressiveMauve) is named by its own message.'''
			if os.path.isfile(alignment):
				logger.error(str(e))
			else:
				logger.error("The alignment file {f} was not found".format(f=alignment))
			raise
		self.finish()
		if self.mask and self.masked.get(self.reference):
			logger.info("{n} SNP positions against {reference} were within {mask} bases of an alignment edge and were not evaluated".format(
					n=self.masked[self.reference], reference=reference, mask=self.mask))
		'''the mask drops those positions when asked, but without it a SNP at an assembly or read
			end must still be seen and reviewed, so it is reported loudly and flagged per row (issue #48)'''
		if self.edge_snps:
			logger.warning("{n} SNP positions within {m} bases of an alignment edge were found against {reference} and are flagged 'edge' in the result files: {ids}".format(
					n=len(self.edge_snps), m=EDGE_MARGIN, reference=reference, ids=", ".join(sorted(self.edge_snps)[:20])))
		return self.SNPS
