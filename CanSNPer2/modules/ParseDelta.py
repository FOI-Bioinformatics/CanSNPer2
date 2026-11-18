'''
ParseDelta reads the delta alignments of MUMmer4 nucmer (filtered with delta-filter) and
reports the query base at every SNP position of the reference. The classification of the
bases, the --mask filter, the edge flag and the result rows are shared with the other
alignment readers, see AlignmentReader.

The delta format gives, per alignment, the reference and query coordinates and the
positions of the indels, but no bases; the first line of the file names the reference and
the query fasta files, which are read for the query bases and for the reference sequence
order and lengths. The database positions run along the reference sequences joined end to
end in file order (the progressiveMauve convention), while delta positions are per sequence.

Each alignment record is "S1 E1 S2 E2 errors similarity_errors stops" followed by one
integer per indel and a terminating 0. For an integer d, |d|-1 aligned pairs come first;
then a positive d is a reference base with a gap in the query, a negative d a query base
with a gap in the reference. S2 > E2 marks an alignment to the reverse strand of the query.
'''

import logging

from CanSNPer2.modules.AlignmentReader import AlignmentReader

logger = logging.getLogger(__name__)

COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}


def read_fasta(path):
	'''[(name, sequence)] in file order, sequence upper case without whitespace'''
	records, name, parts = [], None, []
	with open(path) as fin:
		for line in fin:
			if line.startswith(">"):
				if name is not None:
					records.append((name, "".join(parts).upper()))
				name, parts = line[1:].split()[0] if line[1:].split() else "", []
			else:
				parts.append("".join(line.split()))
	if name is not None:
		records.append((name, "".join(parts).upper()))
	return records


class ParseDelta(AlignmentReader):
	"""Reader of MUMmer4 delta files"""
	def __init__(self, verbose=False, **kwargs):
		super(ParseDelta, self).__init__(verbose=verbose, **kwargs)
		self.offsets = {}
		self.query = {}

	def read_alignment(self, path):
		with open(path) as fin:
			header = fin.readline().split()
			if len(header) < 2:
				raise ValueError("{path} is not a delta file: the first line must name the reference and the query file".format(path=path))
			fin.readline()                              ## NUCMER
			try:
				reference = read_fasta(header[0])
				self.query = dict(read_fasta(header[1]))
			except FileNotFoundError as e:
				raise FileNotFoundError("The sequence files named in {path} are needed to read it: {e}".format(path=path, e=e)) from e
			self.offsets, total = {}, 0
			for name, seq in reference:
				self.offsets[name] = total
				total += len(seq)
			rname = qname = None
			for line in fin:
				if line.startswith(">"):
					rname, qname = line[1:].split()[:2]
					continue
				cols = line.split()
				if len(cols) != 7:
					continue
				s1, e1, s2, e2 = (int(c) for c in cols[:4])
				indels = []
				for value in fin:
					d = int(value)
					if d == 0:
						break
					indels.append(d)
				self.read_record(rname, qname, s1, e1, s2, e2, indels)

	def read_record(self, rname, qname, s1, e1, s2, e2, indels):
		'''report the query base at the SNP positions one alignment covers'''
		offset = self.offsets[rname]
		wanted = self.snps_in(offset + s1, offset + e1)
		if not wanted:
			return
		qseq = self.query[qname]
		step = 1 if s2 <= e2 else -1
		block_length = e1 - s1 + 1
		next_snp = 0
		r, q = s1, s2               ## next reference and query position, per sequence, 1-based

		def emit_pairs(count):
			'''count aligned pairs starting at (r, q)'''
			nonlocal next_snp
			while next_snp < len(wanted) and wanted[next_snp] - offset < r + count:
				k = wanted[next_snp] - offset - r
				qpos = q + step * k
				base = qseq[qpos - 1]
				if step < 0:
					base = COMPLEMENT.get(base, "N")
				self.observe(wanted[next_snp], base, wanted[next_snp] - offset - s1 + 1, block_length, query_name=qname, query_pos=str(qpos))
				next_snp += 1

		for d in indels:
			emit_pairs(abs(d) - 1)
			r += abs(d) - 1
			q += step * (abs(d) - 1)
			if d > 0:
				'''a reference base with a gap in the query'''
				if next_snp < len(wanted) and wanted[next_snp] - offset == r:
					self.observe(wanted[next_snp], "-", r - s1 + 1, block_length, query_name=qname, query_pos="N/A")
					next_snp += 1
				r += 1
			else:
				q += step
			if next_snp == len(wanted):
				return
		emit_pairs(e1 - r + 1)
