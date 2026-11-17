'''
Checks that a reference genome file can hold the SNPs the database places on it.

A reference that was cut short (an interrupted download, a truncated copy) used to be
aligned as it was: every SNP beyond its end had no alignment data, which counts against a
call, so the visible result was missing or wrong genotypes and no error. The database
positions run along the reference sequences joined end to end, so the file must have at
least as many bases as the last SNP position on that reference.
'''

import logging

logger = logging.getLogger(__name__)


def sequence_length(path):
	'''number of bases in a fasta file, all sequences together (whitespace is not counted)'''
	total = 0
	with open(path) as fin:
		for line in fin:
			if not line.startswith(">"):
				total += len("".join(line.split()))
	return total


def last_snp_positions(database):
	'''{reference genome name: the largest SNP position on it} from an open database connection'''
	query = '''SELECT snp_references.genome, max(snp_annotation.position)
					FROM snp_annotation
					JOIN snp_references ON (snp_references.id = snp_annotation.genome_i)
					GROUP BY snp_references.genome'''
	return {genome: int(position) for genome, position in database.query(query).fetchall() if position is not None}


def check_reference(path, name, last_position):
	'''None when the file at path holds last_position bases or more, otherwise the reason'''
	if not last_position:
		return None
	length = sequence_length(path)
	if length >= last_position:
		return None
	return ("reference {name} ({path}) has {length} bases, but the database places SNPs up to position {last}: "
			"the file is truncated or is not the genome the database was built on; delete it and run "
			"CanSNPer2-download again").format(name=name, path=path, length=length, last=last_position)
