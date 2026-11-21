'''
ParseXMFA reads the xmfa alignments of progressiveMauve and reports the query base at
every SNP position of the reference. The classification of the bases, the --mask filter,
the edge flag and the result rows are shared with the other alignment readers, see
AlignmentReader.

progressiveMauve counts reference positions along all reference sequences joined end to
end, in file order. The database stores a SNP as (sequence, position on that sequence), so
the reader converts the SNPs to joined positions with the reference sequence table that
ProgressiveMauve.prepare writes next to the alignment (<alignment>.sequences), and converts
each observed position back before it is recorded.
'''

#__name__="ParseXMFA"
from CanSNPer2.modules.AlignmentReader import AlignmentReader, EDGE_MARGIN, SEQUENCES_SUFFIX, read_sequence_table  # noqa: F401 (EDGE_MARGIN is imported from here by callers)
import os
import logging
from bisect import bisect_left, bisect_right
logger = logging.getLogger(__name__)

reference_genomes = ["FSC200","SCHUS4.1","SCHUS4.2","OSU18","LVS","FTNF002-00"]

class ParseXMFA(AlignmentReader):
	"""Reader of progressiveMauve xmfa files: the reference is the first sequence of every
		block and the query the second."""
	def __init__(self, verbose=False, **kwargs):
		super(ParseXMFA, self).__init__(verbose=verbose, **kwargs)
		'''names and start offsets of the query sequences, see load_query_map'''
		self.query_map = None
		'''sorted joined positions of the SNPs, and the (sequence, position) of each, see
			load_joined_positions'''
		self.joined = []
		self.joined_keys = {}
		### Define translation table for complement base
		self.rcDict = {
			 "A" : "T",
			 "T" : "A",
			 "C" : "G",
			 "G" : "C",
			 "-" : "-",
			 "N" : "N"
		}

	def load_joined_positions(self, table_path):
		'''Convert the SNP positions to positions along the joined reference sequences, with
			the sequence table of the alignment: a sequence starts after the lengths of all
			sequences before it in the table. A SNP on a sequence the table does not have, or
			beyond the end of its sequence, has no joined position and is left out; it would
			otherwise be read from a neighbouring sequence.'''
		offsets, lengths, total = {}, {}, 0
		for name, length in read_sequence_table(table_path):
			offsets[name] = total
			lengths[name] = length
			total += length
		self.joined, self.joined_keys = [], {}
		missing_sequences, outside = set(), []
		for (sequence, position), (_anc, _der, snp_id) in self.snplist.items():
			if sequence not in offsets:
				missing_sequences.add(sequence)
				continue
			if not 1 <= position <= lengths[sequence]:
				outside.append(snp_id)
				continue
			joined = offsets[sequence] + position
			self.joined_keys[joined] = (sequence, position)
		self.joined = sorted(self.joined_keys)
		self.unknown_sequences.update(name for name in offsets if name not in self.snp_positions)
		if missing_sequences:
			logger.warning("SNPs on {sequences} of {reference} are not evaluated: the aligned reference has no such sequence ({table})".format(
					sequences=", ".join(sorted(missing_sequences)), reference=self.reference, table=table_path))
		if outside:
			logger.warning("SNPs {ids} of {reference} lie beyond the end of their sequence in the aligned reference and are not evaluated ({table})".format(
					ids=", ".join(sorted(outside)[:20]), reference=self.reference, table=table_path))

	def joined_snps_in(self, start, end):
		'''the joined SNP positions from start to end, both included, in ascending order'''
		return self.joined[bisect_left(self.joined,start):bisect_right(self.joined,end)]

	def find_snps(self,ref,target,head,query_head=None):
		'''Walk through one aligned sequence pair and report the query base at the SNP
			positions of the reference stretch start..end (both included, issue #39)'''
		start = int(head["start"])
		block_len = int(head["end"]) - start + 1 ## the block length in reference bases, to measure the distance to the far edge (issue #175)
		wanted = [pos - start + 1 for pos in self.joined_snps_in(start, int(head["end"]))] ## 1-based offsets in the block
		if not wanted:
			return
		next_snp = 0
		'''will count the relative position in the reference (without gaps -)'''
		i = 0
		'''same count for the query side, so a SNP can be traced back to a contig and a
			position in the submitted assembly, not only to the reference coordinates (issue #23);
			the walk goes backwards on "-" reference blocks, so the total is needed to turn the
			count into the number of query bases left of the column'''
		q = 0
		q_total = len(target) - target.count("-")
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
			'''query contig and position behind the target base; a target gap has no query
				position (issue #23). q_here counts in the walking direction, which follows the
				reference strand: the query bases left of the column in the block row are q_here
				on a "+" reference block and the rest of the block on a "-" one. A "+" query row
				runs from start to end, a "-" row is reverse complemented and runs from end down
				to start. The result is a position in the query sequences joined end to end,
				which query_location turns into a sequence name and a position in it.'''
			query_name, query_pos = "-", "N/A"
			if q_here is not None and query_head:
				left = q_here if head["sign"] != "-" else q_total - q_here - 1
				if query_head["sign"] == "-":
					joined = int(query_head["end"]) - left
				else:
					joined = int(query_head["start"]) + left
				query_name, query_pos = self.query_location(joined)
			self.observe(*self.joined_keys[start + i - 1], _snp, i, block_len, query_name=query_name, query_pos=query_pos)
			next_snp += 1
			if next_snp == len(wanted):
				break

	def parse_head(self,head):
		'''This help function parses the head of an xmfa sequence line ("> 1:start-end sign
			file") and returns info'''
		cols = head.split(" ") 	## Split header information into columns
		sign = cols[1]			## Save sign of read (orientation of sequence)

		'''Parse out start and end position information of sequence'''
		start,end = list(map(int,cols[0].split(":")[-1].split("-")))
		return {"sign":sign,"start":start,"end":end}

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
			if self.joined_snps_in(refHead["start"], refHead["end"]):
				ref = "".join(refSeq)				 ## Make reference sequence one string
				target = "".join(targetSeq)			 ## Make target sequence one string
				self.find_snps(ref,target,head=refHead,query_head=targetHead)

	def read_query_map(self,path):
		'''names and start offsets of the sequences of a fasta file, in file order, the order
			progressiveMauve joins them in'''
		names, starts, total = [], [], 0
		with open(path) as fin:
			for line in fin:
				if line.startswith(">"):
					names.append(line[1:].split()[0] if line[1:].split() else "")
					starts.append(total)
				else:
					total += len(line.strip())
		return names, starts, total

	def load_query_map(self,xmfa_text,xmfa_path):
		'''The query coordinates of an xmfa file run along all query sequences joined end to
			end. The sequence table that maps them back is read from the file the xmfa header
			names (#Sequence2File, exactly what was aligned), otherwise from the query file the
			caller gave (query_file); without either the query columns are N/A (issue #23)'''
		self.query_map = None
		candidates = []
		for line in xmfa_text[:10000].splitlines():
			if line.startswith("#Sequence2File"):
				candidates.append(line.split("\t",1)[-1].strip())
				break
		if self.query_file:
			candidates.append(self.query_file)
		for path in candidates:
			if path and os.path.isfile(path):
				try:
					self.query_map = self.read_query_map(path)
					return
				except (OSError, UnicodeDecodeError) as e:
					logger.debug("Could not read the query sequences from {path}: {e}".format(path=path,e=e))
		logger.warning("The query sequences of {xmfa} could not be read ({paths}), the query sequence and position columns are N/A".format(
				xmfa=xmfa_path, paths=", ".join(candidates) or "no file given"))

	def query_location(self,joined):
		'''(sequence name, position in it) of a 1-based position in the joined query'''
		if not self.query_map:
			return "-", "N/A"
		names, starts, total = self.query_map
		if not names or not 1 <= joined <= total:
			return "-", "N/A"
		idx = bisect_right(starts, joined - 1) - 1
		return names[idx], str(joined - starts[idx])

	def read_xmfa(self,f):
		'''read xmfa file; each aligned block is separated by a "=" line. The sequence table
			is read after the alignment, so a missing alignment is reported as such.'''
		with open(f) as fin:
			text = fin.read()
		self.load_joined_positions(f + SEQUENCES_SUFFIX)
		self.load_query_map(text, f)
		for seqP in text.strip().split("="):
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
