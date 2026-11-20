'''
Checks that a reference genome file can hold the SNPs the database places on it, and that
the sequences of a fasta file have unique names.

A reference that was cut short (an interrupted download, a truncated copy) used to be
aligned as it was: every SNP beyond its end had no alignment data, which counts against a
call, so the visible result was missing or wrong genotypes and no error. The database
positions run along the reference sequences joined end to end, so the file must have at
least as many bases as the last SNP position on that reference.

The length alone does not show that the sequences are in the order the positions were
counted in: a reference with its chromosomes in another order (or another genome of about
the same length) passed, and every aligner then read the bases of the wrong positions. The
reference base at each SNP position must therefore be the ancestral or the derived base of
that SNP. A correct reference meets this at every SNP; one with its sequences in another
order meets it at about half of them by chance, so a reference that fails at a fifth of
its SNPs or more is refused, and fewer mismatches (an error in the scheme) are reported as
a warning.
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


'''a reference is refused when at least this fraction of its SNPs has neither the ancestral nor
	the derived base in it'''
MISMATCH_FRACTION = 0.2


def reference_snps(database):
	'''{reference genome name: [(snp_id, position, ancestral base, derived base), ...]} from an
		open database connection'''
	query = '''SELECT snp_references.genome, snp_annotation.snp_id, snp_annotation.position,
					snp_annotation.ancestral_base, snp_annotation.derived_base
					FROM snp_annotation
					JOIN snp_references ON (snp_references.id = snp_annotation.genome_i)'''
	snps = {}
	for genome, snp_id, position, ancestral, derived in database.query(query).fetchall():
		if position is None:
			continue
		snps.setdefault(genome, []).append((snp_id, int(position), str(ancestral).upper(), str(derived).upper()))
	return snps


def joined_sequence(path):
	'''the sequences of a fasta file joined end to end in file order, upper case, without
		whitespace'''
	parts = []
	with open(path) as fin:
		for line in fin:
			if not line.startswith(">"):
				parts.append("".join(line.split()).upper())
	return "".join(parts)


def base_mismatches(path, snps):
	'''the SNPs whose position in the reference file at path holds neither the ancestral nor
		the derived base: [(snp_id, position, reference base, ancestral, derived), ...]; a
		position beyond the end is left to the length check'''
	sequence = joined_sequence(path)
	mismatches = []
	for snp_id, position, ancestral, derived in snps:
		if 0 < position <= len(sequence):
			base = sequence[position - 1]
			if base != ancestral and base != derived:
				mismatches.append((snp_id, position, base, ancestral, derived))
	return mismatches


def sequence_table(path):
	'''[(name, length), ...] of the sequences of a fasta file, in file order'''
	table = []
	with open(path) as fin:
		for line in fin:
			if line.startswith(">"):
				words = line[1:].split()
				table.append([words[0] if words else "", 0])
			elif table:
				table[-1][1] += len("".join(line.split()))
	return [(name, length) for name, length in table]


def locate(table, position):
	'''(sequence name, position on that sequence) of a position on the joined reference, or
		("-", "-") when it is not a position inside the sequences of table'''
	try:
		position = int(position)
	except (TypeError, ValueError):
		return "-", "-"
	if position < 1:
		return "-", "-"
	start = 0
	for name, length in table:
		if position <= start + length:
			return name, position - start
		start += length
	return "-", "-"


def check_bases(path, name, snps):
	'''None when the reference file at path has the ancestral or the derived base at enough of
		its SNPs, otherwise the reason; mismatches below the limit are logged as a warning'''
	if not snps:
		return None
	mismatches = base_mismatches(path, snps)
	if not mismatches:
		return None
	listed = ", ".join("{snp} at {pos} has {base} (ancestral {anc}, derived {der})".format(snp=snp, pos=pos, base=base, anc=anc, der=der)
			for snp, pos, base, anc, der in mismatches[:5]) + (" ..." if len(mismatches) > 5 else "")
	share = "{n} of {m} SNPs".format(n=len(mismatches), m=len(snps))
	if len(mismatches) >= MISMATCH_FRACTION * len(snps):
		order = ", ".join(seq for seq, _length in sequence_table(path))
		return ("reference {name} ({path}) has neither the ancestral nor the derived base at {share} ({listed}): "
				"the reference base at a SNP position must be one of them. The sequences of the file ({order}) are "
				"joined in this order to count the positions, so they must be in the order of the genome the "
				"database was built on; delete the file and run CanSNPer2-download again").format(
				name=name, path=path, share=share, listed=listed, order=order)
	logger.warning("reference {name} has neither the ancestral nor the derived base at {share}: {listed}; the scheme may give a wrong base for these SNPs".format(
			name=name, share=share, listed=listed))
	return None


def check_reference(path, name, last_position, snps=None):
	'''None when the reference file at path has unique sequence names, holds last_position
		bases or more and has the ancestral or the derived base at its SNPs (snps:
		[(snp_id, position, ancestral, derived), ...]), otherwise the reason'''
	names = sequence_names(path)
	problem = duplicate_names(names)
	if problem:
		return "reference {name} ({path}) has {problem}".format(name=name, path=path, problem=problem)
	if last_position:
		length = sequence_length(path)
		if length < last_position:
			return too_short(path, name, length, last_position)
	return check_bases(path, name, snps)


def too_short(path, name, length, last_position):
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
