#!/usr/bin/env python3 -c

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
		self.export = kwargs["export"]
		self.called_snps = []

		## as the snps are ordered according to the reference mauve positions they can be used sequentialy without search
		'''SNPs will be stored as a sorted set containing (position, refBase, targetBase,SNPname)'''
		self.SNPS = {}
		self.SNP_info = []

		'''Mask option, for distant relatives false SNPs near edges of alignments may become a problem,
			this option allows masking of SNPs placed within n bases of the edge of an alignment'''
		#self.mask = kwargs["mask"]
		kwargs["verbose"] = True
		if "database" in kwargs:
			self.database = XMFAFunctions(kwargs["database"],verbose=self.verbose)
		else:
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

	def reverse_complement(self,dna):
		'''Complement and reverse DNA string'''
		dna_rev = [ self.rcDict[x] for x in dna[::-1] ]
		return "".join(dna_rev)

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
		'''will count the relative position in the reference (without gaps -)'''
		i = 0
		''' ii is the actual position in the sequence'''
		base_positions = range(len(ref))

		'''	If sign is "-" it means the reference sequence is the
			reverse complement, hence positions need to be reversed
		'''
		if head["sign"] == "-":
			base_positions = reversed(base_positions)
		'''Walk through the positions and check whenever the position matches a SNP'''
		for ii in base_positions:
			if ref[ii] != "-":          	## if the reference contains a gap "-" do not count
				i+=1
			if i == snppos:                	## if current possition contains a snp
				SNP[snp_id] = 0
				_snp = target[ii]        	## get base in target sequence
				if head["sign"] == "-":
					'''If the sequence sign is "-" the complement base needs to be retrieved'''
					_snp = self.rcDict[_snp]
				'''Fetch information about snp to allow print to file'''
				orig_snp_pos = str(self.snplist[self.current_snp][0])
				snpinfo = [snp_id,self.reference,orig_snp_pos,rbase,tbase,_snp]
				self.SNP_info.append(snpinfo)

				if tbase == _snp:           ## SNP is confirmed to be Derived
					SNP[snp_id] = 1 		## Derived SNP
					self.called_snps.append(snp_id)
				elif rbase == _snp:  		## SNP is confirmed to be ancestral
					SNP[snp_id] = 2			## Ancestral SNP
				elif _snp != "-" or _snp != "N":
					SNP[snp_id] = 3			## Other base than ancestral or derived was found
				self._next_pos()  ## SNP found get next
				snppos,rbase,tbase,snp_id = self._snpinfo(head)    	## Retrieve all information known about the next upcoming SNP
		return SNP

	def parse_head(self,head):
		'''This help function parses the head of an xmfa file and returns info'''
		cols = head.split(" ") 	## Split header information into columns
		sign = cols[1]			## Save sign of read (orientation of sequence)

		'''Parse out start and end position information of sequence'''
		start,end = list(map(int,cols[0].split(":")[-1].split("-")))
		return {"sign":sign,"start":start,"end":end}

	def read_sequence(self,seqP,res={}):
		'''Read information in sequence pair'''
		seqLines = seqP.strip().split("> ") 	### Split pair into two sequences
		headinfo = seqLines[0]					### Get header info of xmfa file (All comments)
		if len(seqLines) > 2:  					### Both target and reference have sequence
			refSeq = seqLines[1].split("\n")				## reference sequence
			refHead = self.parse_head(refSeq.pop(0))		## parse reference header info
			targetSeq = seqLines[2].split("\n")				## target sequence
			targetHead = self.parse_head(targetSeq.pop(0))	## parse target sequence header
			'''Parse aligned sequence pair '''
			while int(self.current_snp) < int(refHead["start"]): ## Current SNP not aligned
				if not self._next_pos():
					break
			if refHead["start"] < self.current_snp and refHead["end"] > self.current_snp:
				'''Check if current snp is within this mapped region else skip'''
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
			return False
		return self.SNPS

	def run(self, xmfa, organism,reference=False,database=False):
		'''Parse XMFA file and return SNPS matching the given database'''
		self.SNPS = {}  ## For each run SNPs has to be emtpy
		'''Create connection to SNP database if it is not connected'''
		if not self.database:
			self.database = XMFAFunctions(database)
		''' retrieve registered SNPs'''
		if not reference:
			reference = os.path.basename(xmfa).split("_")[0]
		self.reference=reference
		self.snplist, self.snp_positions = self.database.get_snps(organism,reference)
		#if self.verbose: print(self.snplist)
		'''save first snp to look for'''
		self.current_snp = self.snp_positions.pop(0)
		snps = self.read_xmfa(xmfa)
		return snps

if __name__=="__main__":
	import argparse

	parser = argparse.ArgumentParser(description='Parse xmfa files and extract non overlapping sequences')

	parser.add_argument('xmfa', 		metavar='', help='fasta xmfa to be parsed')
	parser.add_argument('database', 	metavar='', help='CanSNP database')
	parser.add_argument('--organism', 	metavar='', default="Francisella", 	help="Specify organism")
	parser.add_argument('--reference', 	metavar='', default="FSC200",		help="Specify reference", choices=reference_genomes)
	parser.add_argument('--mask', 		metavar='', default=0, type=int, 	help="Mask regions near end of alignments (nbases)")
	parser.add_argument('--verbose',	action='store_true',help="print process info, default no output")

	args = parser.parse_args()

	logger.debug(args)
	xmfa = ParseXMFA(verbose=args.verbose,mask=args.mask)
	SNPS = xmfa.run(database=args.database, xmfa=args.xmfa, organism=args.organism,reference=args.reference)
	logger.info(SNPS)
