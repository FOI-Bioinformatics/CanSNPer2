'''
ParseXMFA reads the xmfa alignments of progressiveMauve and reports the query base at
every SNP position of the reference. The classification of the bases, the --mask filter,
the edge flag and the result rows are shared with the other alignment readers, see
AlignmentReader.
'''

#__name__="ParseXMFA"
from CanSNPer2.modules.AlignmentReader import AlignmentReader, EDGE_MARGIN  # noqa: F401 (EDGE_MARGIN is imported from here by callers)
import logging
logger = logging.getLogger(__name__)

reference_genomes = ["FSC200","SCHUS4.1","SCHUS4.2","OSU18","LVS","FTNF002-00"]

class ParseXMFA(AlignmentReader):
	"""Reader of progressiveMauve xmfa files: the reference is the first sequence of every
		block and the query the second."""
	def __init__(self, verbose=False, **kwargs):
		super(ParseXMFA, self).__init__(verbose=verbose, **kwargs)
		### Define translation table for complement base
		self.rcDict = {
			 "A" : "T",
			 "T" : "A",
			 "C" : "G",
			 "G" : "C",
			 "-" : "-",
			 "N" : "N"
		}

	def find_snps(self,ref,target,head,query_head=None):
		'''Walk through one aligned sequence pair and report the query base at the SNP
			positions of the reference stretch start..end (both included, issue #39)'''
		start = int(head["start"])
		block_len = int(head["end"]) - start + 1 ## the block length in reference bases, to measure the distance to the far edge (issue #175)
		wanted = [pos - start + 1 for pos in self.snps_in(start, int(head["end"]))] ## 1-based offsets in the block
		if not wanted:
			return
		next_snp = 0
		'''will count the relative position in the reference (without gaps -)'''
		i = 0
		'''same count for the query side, so a SNP can be traced back to a contig and a
			position in the submitted assembly, not only to the reference coordinates (issue #23);
			the walk goes backwards on "-" blocks, so the total is needed to convert the
			right-to-left count into the left-to-left position the contig coordinates use'''
		q = 0
		q_total = len(target) - target.count("-")
		query_name = query_head["name"] if query_head and query_head.get("name") else "-"
		''' ii is the actual position in the sequence'''
		base_positions = range(len(ref))
		'''	If sign is "-" it means the reference sequence is the
			reverse complement, hence positions need to be reversed
		'''
		if head["sign"] == "-":
			base_positions = reversed(base_positions)
		'''Walk through the positions and check whenever the position matches a SNP'''
		for ii in base_positions:
			if ii >= len(target): ## target sequence shorter than the reference block, the rest is unaligned (issue #97)
				break
			'''q_here is the number of query bases before this column, counted in the walking
				direction, which maps the column to a position in the query contig; a gap column
				has no query base and no position (issue #23)'''
			if target[ii] != "-":
				q_here = q
				q += 1
			else:
				q_here = None
			if ref[ii] == "-":          	## if the reference contains a gap "-" do not count
				continue
			i += 1
			if i != wanted[next_snp]:
				continue
			_snp = target[ii].upper()        	## get base in target sequence make sure it is upper case!
			if head["sign"] == "-":
				'''If the sequence sign is "-" the complement base needs to be retrieved'''
				_snp = self.rcDict.get(_snp,"N") ## unknown base characters become N (issue #96)
			'''query contig and position behind the target base: the walk above counted the
				query bases up to this column, and the block header gives the contig name and
				the coordinate the block starts at; a target gap has no query position (issue #23)'''
			if q_here is None or not query_head:
				query_pos = "N/A"
			elif query_head["sign"] == "-":
				query_pos = str(int(query_head["end"]) - (q_total - q_here - 1))
			else:
				query_pos = str(int(query_head["start"]) + q_here)
			self.observe(start + i - 1, _snp, i, block_len, query_name=query_name, query_pos=query_pos)
			next_snp += 1
			if next_snp == len(wanted):
				break

	def parse_head(self,head):
		'''This help function parses the head of an xmfa sequence line and returns info.
			The mauve layout is "> 1:start-end sign name", the name being the contig the
			block comes from, needed to trace calls back to the assembly (issue #23).'''
		cols = head.split(" ") 	## Split header information into columns
		sign = cols[1]			## Save sign of read (orientation of sequence)

		'''Parse out start and end position information of sequence'''
		start,end = list(map(int,cols[0].split(":")[-1].split("-")))
		name = " ".join(cols[2:]).strip() if len(cols) > 2 else ""
		return {"sign":sign,"start":start,"end":end,"name":name}

	def read_sequence(self,seqP):
		'''Read one block of the xmfa file; blocks with only one of the two genomes are
			unaligned sequence and carry no SNP information'''
		seqLines = seqP.strip().split("> ") 	### Split pair into two sequences
		if len(seqLines) > 2:  					### Both target and reference have sequence
			refSeq = seqLines[1].split("\n")				## reference sequence
			targetSeq = seqLines[2].split("\n")				## target sequence
			try:
				refHead = self.parse_head(refSeq.pop(0))		## parse reference header info
				targetHead = self.parse_head(targetSeq.pop(0))	## parse target sequence header
			except (ValueError, IndexError): ## malformed sequence block header, skip this alignment pair (issue #98)
				logger.warning("Malformed xmfa sequence block header, skipping block: {head}".format(head=(refSeq[:1] or [""])[0]))
				return
			'''Parse aligned sequence pair, if any SNP is inside it'''
			if self.snps_in(refHead["start"], refHead["end"]):
				ref = "".join(refSeq)				 ## Make reference sequence one string
				target = "".join(targetSeq)			 ## Make target sequence one string
				'''the target head gives the query contig name for the result columns (issue #23)'''
				self.find_snps(ref,target,head=refHead,query_head=targetHead)

	def read_xmfa(self,f):
		'''read xmfa file; each aligned block is separated by a "=" line'''
		with open(f) as fin:
			for seqP in fin.read().strip().split("="):
				self.read_sequence(seqP)

	def read_alignment(self, path):
		self.read_xmfa(path)

	def run(self, xmfa, reference=False,database=False):
		'''Parse XMFA file and return SNPS matching the given database'''
		return super(ParseXMFA, self).run(xmfa, reference=reference, database=database)

if __name__=="__main__":
	import argparse

	parser = argparse.ArgumentParser(description='Parse xmfa files and extract non overlapping sequences')

	parser.add_argument('xmfa', 		metavar='', help='fasta xmfa to be parsed')
	parser.add_argument('database', 	metavar='', help='CanSNP database')
	parser.add_argument('--reference', 	metavar='', default="",		help="Specify reference", choices=reference_genomes)
	parser.add_argument('--mask', 		metavar='', default=0, type=int, 	help="Mask regions near end of alignments (nbases)")
	parser.add_argument('--verbose',	action='store_true',help="print process info, default no output")

	args = parser.parse_args()

	logger.debug(args)
	xmfa = ParseXMFA(verbose=args.verbose,mask=args.mask)
	SNPS = xmfa.run(database=args.database, xmfa=args.xmfa, reference=args.reference)
	logger.info(SNPS)
