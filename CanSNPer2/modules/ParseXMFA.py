'''
ParseXMFA is a script that parse xmfa files and extract all SNPs
	The script allows an option using a flanking score that limits
	SNPs near edges of a block to have a impact on classifying decisions.
'''

#__name__="ParseXMFA"
from CanSNPer2.modules.DatabaseConnection import XMFAFunctions
import os
import logging
logger = logging.getLogger(__name__)

reference_genomes = ["FSC200","SCHUS4.1","SCHUS4.2","OSU18","LVS","FTNF002-00"]

'''SNPs closer to either edge of an alignment block than this many reference bases are
	flagged in the output: assembly ends and read ends are the most likely place for
	alignment errors, and a SNP sitting there should be reviewed by the user instead of
	silently trusted or silently dropped (--mask, issue #48)'''
EDGE_MARGIN = 25

class ParseXMFA(object):
	"""docstring for ParseXMFA."""
	def __init__(self, verbose=False, **kwargs):
		super(ParseXMFA, self).__init__()
		### Define translation table for complement base
		self.rcDict = {
			 "A" : "T",
			 "T" : "A",
			 "C" : "G",
			 "G" : "C",
			 "-" : "-",
			 "N" : "N"
		}

		'''Define all base variables'''
		self.verbose = verbose
		self.called_snps = []

		'''SNP state is initialized here so _next_pos/_snpinfo never hit an
			unbound attribute before run() fills them (issues #48, #49)'''
		self.snplist = {}
		self.snp_positions = []
		self.current_snp = False

		## as the snps are ordered according to the reference mauve positions they can be used sequentialy without search
		'''SNPs will be stored as a sorted set containing (position, refBase, targetBase,SNPname)'''
		self.SNPS = {}
		self.SNP_info = []

		'''Mask option, for distant relatives false SNPs near edges of alignments may become a problem,
			this option allows masking of SNPs placed within n bases of the edge of an alignment'''
		self.mask = int(kwargs.get("mask",0) or 0) ## the option was advertised but never read (issue #175)
		self.masked = {} ## how many positions each reference lost to the mask
		self.edge_snps = set() ## SNPs within EDGE_MARGIN bases of an alignment edge, flagged not dropped (issue #48)
		kwargs["verbose"] = True
		'''the database is only connected in run(), by whichever process uses this object:
			the SNP calling workers are forked from the parent and used to inherit its open
			sqlite connection and cursor, so several processes read through one cursor (issue #49)'''
		self.database_path = kwargs.get("database",False)
		self.database = False

	def get_snp_info(self):
		'''Return SNP_info'''
		return self.SNP_info

	def get_snps(self):
		'''Return snps'''
		return self.SNPS

	def get_called_snps(self):
		'''Return true snps'''
		return self.called_snps

	def _next_pos(self):
		'''return the next SNP position if list is not empty'''
		if len(self.snp_positions) > 0:
			self.current_snp = self.snp_positions.pop(0)	## get next SNP
		return

	def _snpinfo(self,head):
		snppos,rbase,tbase,snp_id = self.snplist[self.current_snp]    	## Retrieve all information known about the next upcoming SNP
		snppos -= int(head["start"])-1  	## python counts from 0 not 1
		return snppos,rbase,tbase,snp_id

	def find_snps(self,ref,target,head=0):
		'''Walk through the paired sequences and save SNPs'''
		snppos,rbase,tbase,snp_id = self._snpinfo(head)    	## Retrieve all information known about the next upcoming SNP
		SNP = {snp_id:0} ## SNP not found
		masked_ids = []
		'''will count the relative position in the reference (without gaps -)'''
		i = 0
		''' ii is the actual position in the sequence'''
		base_positions = range(len(ref))
		'''the block length in reference bases, to measure the distance to the far edge (issue #175)'''
		block_len = int(head["end"]) - int(head["start"]) + 1

		'''	If sign is "-" it means the reference sequence is the
			reverse complement, hence positions need to be reversed
		'''
		if head["sign"] == "-":
			base_positions = reversed(base_positions)
		'''Walk through the positions and check whenever the position matches a SNP'''
		for ii in base_positions:
			if ii >= len(target): ## target sequence shorter than the reference block, the rest is unaligned (issue #97)
				break
			if ref[ii] != "-":          	## if the reference contains a gap "-" do not count
				i+=1
			if i == snppos:                	## if current possition contains a snp
				if self.mask and (snppos <= self.mask or block_len - snppos < self.mask):
					'''the flanking mask: a SNP this close to either edge of the alignment is
						not trustworthy, and it is dropped completely: counting it as absent
						would penalise the sample for a base that was never compared (issue #175)'''
					masked_ids.append(snp_id)
					self.masked[self.reference] = self.masked.get(self.reference,0)+1
					self._next_pos()  ## SNP handled, get next
					snppos,rbase,tbase,snp_id = self._snpinfo(head)    	## Retrieve all information known about the next upcoming SNP
					continue
				SNP[snp_id] = 0
				_snp = target[ii].upper()        	## get base in target sequence make sure it is upper case!
				if head["sign"] == "-":
					'''If the sequence sign is "-" the complement base needs to be retrieved'''
					_snp = self.rcDict.get(_snp,"N") ## unknown base characters become N (issue #96)
				'''Fetch information about snp to allow print to file; a SNP this close to either
					edge of the alignment block is still evaluated but carries the "edge" flag in
					the output, so the user can review the positions that sit on the ends of the
					assembly or of an alignment, where errors are most likely (issue #48)'''
				orig_snp_pos = str(self.snplist[self.current_snp][0])
				edge_flag = "edge" if (snppos <= EDGE_MARGIN or block_len - snppos < EDGE_MARGIN) else ""
				snpinfo = [snp_id,self.reference,orig_snp_pos,rbase,tbase,_snp,edge_flag]
				if edge_flag:
					self.edge_snps.add(snp_id)
				self.SNP_info.append(snpinfo)

				if tbase == _snp:           ## SNP is confirmed to be Derived
					SNP[snp_id] = 1 		## Derived SNP
					self.called_snps.append(snp_id)
				elif rbase == _snp:  		## SNP is confirmed to be ancestral
					SNP[snp_id] = 2			## Ancestral SNP
				elif _snp != "-" and _snp != "N":
					SNP[snp_id] = 3			## Other base than ancestral or derived was found
				self._next_pos()  ## SNP found get next
				snppos,rbase,tbase,snp_id = self._snpinfo(head)    	## Retrieve all information known about the next upcoming SNP
		'''masked positions are not part of this alignment result at all (issue #175)'''
		if masked_ids:
			return {id:state for id,state in SNP.items() if id not in masked_ids}
		return SNP

	def parse_head(self,head):
		'''This help function parses the head of an xmfa file and returns info'''
		cols = head.split(" ") 	## Split header information into columns
		sign = cols[1]			## Save sign of read (orientation of sequence)

		'''Parse out start and end position information of sequence'''
		start,end = list(map(int,cols[0].split(":")[-1].split("-")))
		return {"sign":sign,"start":start,"end":end}

	def read_sequence(self,seqP,res=None):
		'''Read information in sequence pair'''
		if res is None:
			res = {}
		seqLines = seqP.strip().split("> ") 	### Split pair into two sequences
		headinfo = seqLines[0]					### Get header info of xmfa file (All comments)
		if len(seqLines) > 2:  					### Both target and reference have sequence
			refSeq = seqLines[1].split("\n")				## reference sequence
			targetSeq = seqLines[2].split("\n")				## target sequence
			try:
				refHead = self.parse_head(refSeq.pop(0))		## parse reference header info
				targetHead = self.parse_head(targetSeq.pop(0))	## parse target sequence header
			except (ValueError, IndexError): ## malformed sequence block header, skip this alignment pair (issue #98)
				logger.warning("Malformed xmfa sequence block header, skipping block: {head}".format(head=(refSeq[:1] or [""])[0]))
				return res
			'''Parse aligned sequence pair '''
			while int(self.current_snp) < int(refHead["start"]): ## Current SNP not aligned
				# Keep iterating new canSNP-positions until the end of current alignment or until there are no more canSNPs to check
				self._next_pos() # iterate snp position
				if self.current_snp > refHead["end"] or len(self.snp_positions) == 0: # if the new snp position is beyond alignment range, then break. Also break if we run out of snp positions to check
					break
				#/

			if refHead["start"] <= self.current_snp and refHead["end"] >= self.current_snp:
				'''Check if current snp is within this mapped region else skip; the block
					covers start..end inclusive (block_len counts both ends, issue #175), so a
					SNP exactly at the first or last position is inside too and the whole block
					must not be skipped for it (issue #39)'''
				ref = "".join(refSeq)				 ## Make reference sequence one string
				target = "".join(targetSeq)			 ## Make target sequence one string
				'''For each snp within this region find it and merge with all others'''
				res = dict(**res, **self.find_snps(ref,target,head=refHead))
		return res

	def read_xmfa(self,f=False):
		'''read xmfa file'''
		if not f:
			f = self.xmfa
		try:
			with open(f) as fin:
				#### Each aligned sequence part is separated by = sign, so start with splitting the sequences into pairs of sequence
				seqPairs = fin.read().strip().split("=")
				for seqP in seqPairs:
					### Join together all SNPs found in data
					self.SNPS = dict(**self.SNPS, **self.read_sequence(seqP))
		except FileNotFoundError:
			'''this used to return False, which every caller read as "no SNPs in this
				alignment", so a missing alignment typed the sample NA and exited 0
				(issue #176, #179)'''
			logger.error("The alignment file {f} was not found".format(f=f))
			raise
		return self.SNPS

	def run(self, xmfa, reference=False,database=False):
		'''Parse XMFA file and return SNPS matching the given database'''
		self.SNPS = {}  ## For each run SNPs has to be emtpy
		self.edge_snps = set() ## one reporting set per parsed alignment (issue #48)
		'''Create the connection to the SNP database here, in the process that parses the
			alignment, so a forked worker never reads through a cursor opened by its parent
			(issue #49)'''
		if not self.database:
			database = self.database_path or database
			if not database:
				raise ValueError("No database was provided to ParseXMFA, cannot retrieve snps from the database!")
			self.database = XMFAFunctions(database,verbose=self.verbose)
		''' retrieve registered SNPs'''
		if not reference:
			reference = os.path.basename(xmfa).split("_")[0]
		self.reference=reference
		self.snplist, self.snp_positions = self.database.get_snps(reference)
		#if self.verbose: print(self.snplist)
		'''save first snp to look for'''
		if not self.snp_positions: ## the reference has no registered SNPs in this database, nothing to look for (issue #99)
			logger.warning("No SNPs are registered in the database for reference {reference}, no SNP calls possible".format(reference=reference))
			return self.SNPS
		self.current_snp = self.snp_positions.pop(0)
		snps = self.read_xmfa(xmfa)
		if self.mask and self.masked.get(self.reference):
			logger.info("{n} SNP positions against {reference} were within {mask} bases of an alignment edge and were not evaluated".format(
					n=self.masked[self.reference], reference=reference, mask=self.mask))
		'''the mask drops those positions when asked, but without it a SNP at an assembly or read
			end must still be seen and reviewed, so it is reported loudly and flagged per row (issue #48)'''
		if self.edge_snps:
			logger.warning("{n} SNP positions within {m} bases of an alignment edge were found against {reference} and are flagged 'edge' in the result files: {ids}".format(
					n=len(self.edge_snps), m=EDGE_MARGIN, reference=reference, ids=", ".join(sorted(self.edge_snps)[:20])))
		return snps

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
