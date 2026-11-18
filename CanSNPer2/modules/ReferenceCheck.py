'''
Checks that a reference genome file can hold the SNPs the database places on it, and that
the sequences of a fasta file have unique names.

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
	'''None when the reference file at path has unique sequence names and holds last_position
		bases or more, otherwise the reason'''
	names = sequence_names(path)
	problem = duplicate_names(names)
	if problem:
		return "reference {name} ({path}) has {problem}".format(name=name, path=path, problem=problem)
	if not last_position:
		return None
	length = sequence_length(path)
	if length >= last_position:
		return None
	return ("reference {name} ({path}) has {length} bases, but the database places SNPs up to position {last}: "
			"the file is truncated or is not the genome the database was built on; delete it and run "
			"CanSNPer2-download again").format(name=name, path=path, length=length, last=last_position)


def sequence_names(path):
	'''the name of every sequence in a fasta file, in file order: the first word of the header
		line, which is how minimap2, MUMmer and the result files identify a sequence'''
	names = []
	with open(path) as fin:
		for line in fin:
			if line.startswith(">"):
				words = line[1:].split()
				names.append(words[0] if words else "")
	return names


def duplicate_names(names):
	'''None when every name is set and unique, otherwise the reason. Two sequences with the
		same first word are one sequence to the aligners: the reference table of the minimap2
		reader mapped both to one offset, which lost SNPs and read others from the wrong
		sequence without an error'''
	seen, duplicates, empty = set(), [], 0
	for name in names:
		if not name:
			empty += 1
		elif name in seen and name not in duplicates:
			duplicates.append(name)
		seen.add(name)
	problems = []
	if duplicates:
		problems.append("more than one sequence named {names}".format(names=", ".join(duplicates[:5]) + (" ..." if len(duplicates) > 5 else "")))
	if empty:
		problems.append("{n} sequence(s) without a name".format(n=empty))
	if not problems:
		return None
	return "{problems}; sequence names (the first word of each header line) must be set and unique".format(problems=" and ".join(problems))
