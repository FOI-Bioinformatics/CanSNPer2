'''
ParsePAF reads the PAF alignments of minimap2 (run with -c --cs=long) and reports the
query base at every SNP position of the reference. The classification of the bases, the
--mask filter, the edge flag and the result rows are shared with the other alignment
readers, see AlignmentReader.

The database positions run along the reference sequences joined end to end in file order
(the progressiveMauve convention), while PAF gives positions within each reference
sequence. The aligner writes the reference sequence table next to the alignment
(<alignment>.sequences: name and length per line, in file order) so the alignment can be
read again later (--skip_alignment) without the reference file.

Which alignments are used: only primary alignments (tp:A:P) with a mapping quality of at
least min_mapq. A secondary alignment or a primary one with a low mapping quality is a
repeat or paralog, whose base can belong to another copy of the region. A reference
position covered by two primary alignments with different query bases is reported as N
(AlignmentReader.observe).
'''

import re
import logging

from CanSNPer2.modules.AlignmentReader import AlignmentReader

logger = logging.getLogger(__name__)

'''the long cs tag: "=ACGT" identical bases, "*ag" reference a / query g, "+acg" bases only
	in the query, "-acg" bases only in the reference, "~gt12ag" an intron (not produced by the
	assembly presets, skipped as reference without query)'''
CS_OP = re.compile(r"=([A-Za-z]+)|\*([a-z])([a-z])|\+([a-z]+)|-([a-z]+)|~[a-z]{2}(\d+)[a-z]{2}")

SEQUENCES_SUFFIX = ".sequences"


def write_sequence_table(fasta, path):
	'''write name and length of the sequences of a fasta file, in file order'''
	rows = []
	with open(fasta) as fin:
		for line in fin:
			if line.startswith(">"):
				rows.append([line[1:].split()[0] if line[1:].split() else "", 0])
			elif rows:
				rows[-1][1] += len(line.strip())
	with open(path, "w") as out:
		for name, length in rows:
			out.write("{0}\t{1}\n".format(name, length))


class ParsePAF(AlignmentReader):
	"""Reader of minimap2 PAF files with the long cs tag; the reference is the target"""
	def __init__(self, verbose=False, min_mapq=20, **kwargs):
		super(ParsePAF, self).__init__(verbose=verbose, **kwargs)
		self.min_mapq = int(min_mapq)
		self.offsets = {}

	def read_sequence_table(self, alignment):
		'''start offset of each reference sequence in the joined reference'''
		path = alignment + SEQUENCES_SUFFIX
		offsets, total = {}, 0
		try:
			with open(path) as fin:
				for line in fin:
					name, length = line.rstrip("\n").split("\t")
					offsets[name] = total
					total += int(length)
		except FileNotFoundError as e:
			raise FileNotFoundError("The reference sequence table {path} of the alignment is missing, it is written when the alignment is made".format(path=path)) from e
		return offsets

	def read_alignment(self, path):
		self.offsets = self.read_sequence_table(path)
		used = skipped = 0
		with open(path) as fin:
			for line in fin:
				fields = line.rstrip("\n").split("\t")
				if len(fields) < 12:
					continue
				tags = fields[12:]
				if "tp:A:P" not in tags or int(fields[11]) < self.min_mapq:
					skipped += 1
					continue
				cs = next((tag[5:] for tag in tags if tag.startswith("cs:Z:")), None)
				if cs is None:
					raise ValueError("{path} has no cs tag, minimap2 must be run with -c --cs=long".format(path=path))
				if fields[5] not in self.offsets:
					raise ValueError("Reference sequence {name} of {path} is not in its sequence table".format(name=fields[5], path=path))
				used += 1
				self.read_record(fields, cs)
		logger.debug("{path}: {used} alignments used, {skipped} secondary or below mapping quality {q} skipped".format(
				path=path, used=used, skipped=skipped, q=self.min_mapq))

	def read_record(self, fields, cs):
		'''report the query base at the SNP positions one alignment covers. The cs tag is on
			the forward strand of the reference, with the query bases complemented for a "-"
			alignment, so its bases need no further complement'''
		qname, qstart, qend, strand = fields[0], int(fields[2]), int(fields[3]), fields[4]
		start = self.offsets[fields[5]] + int(fields[7]) + 1   ## 1-based, joined reference
		end = self.offsets[fields[5]] + int(fields[8])
		wanted = self.snps_in(start, end)
		if not wanted:
			return
		block_length = end - start + 1
		next_snp = 0
		at = start       ## reference position of the next reference base in the cs string
		qdone = 0        ## query bases consumed so far
		def query_pos(consumed):
			'''position in the query sequence of the query base after "consumed" bases'''
			return str(qstart + consumed + 1) if strand == "+" else str(qend - consumed)
		for match in CS_OP.finditer(cs):
			same, rbase, qbase, insertion, deletion, intron = match.groups()
			if insertion:
				qdone += len(insertion)
				continue
			if same:
				span, query_bases = len(same), same
			elif rbase:
				span, query_bases = 1, qbase
			elif deletion:
				span, query_bases = len(deletion), None
			else:
				span, query_bases = int(intron), None
			while next_snp < len(wanted) and wanted[next_snp] < at + span:
				pos = wanted[next_snp]
				k = pos - at
				if query_bases is None:
					self.observe(pos, "-", pos - start + 1, block_length, query_name=qname, query_pos="N/A")
				else:
					self.observe(pos, query_bases[k], pos - start + 1, block_length, query_name=qname, query_pos=query_pos(qdone + k))
				next_snp += 1
			at += span
			if query_bases is not None:
				qdone += span
			if next_snp == len(wanted):
				break
