'''
Checks that a reference genome file can hold the SNPs the database places on it, and that
the sequences of a fasta file have unique names.

The database stores each SNP as a reference sequence (the first word of its fasta header,
a GenBank accession for the published schemes) and a 1-based position on that sequence.
Each SNP is therefore looked up on its own sequence: the order of the sequences in the file
and sequences without SNPs (plasmids) do not matter.

A reference file is refused when
- a sequence name is empty or used twice;
- a sequence the database names is missing (a file with RefSeq names for a scheme built on
  GenBank accessions gave NA calls without an error);
- a sequence is shorter than its largest SNP position (an interrupted download, a
  truncated copy: the SNPs beyond its end had no alignment data, which counts against a
  call);
- the reference base is neither the ancestral nor the derived base at a fifth of the SNPs
  of the reference or more. A correct reference meets this at every SNP; another genome or
  another version of a sequence fails at many of them. Fewer mismatches (an error in the
  scheme) are reported as a warning.
'''

import logging

logger = logging.getLogger(__name__)

'''a reference is refused when at least this fraction of its SNPs has neither the ancestral nor
	the derived base in it'''
MISMATCH_FRACTION = 0.2

'''names listed in a message before it is cut short'''
LISTED = 5


def listing(names):
	'''the first LISTED names, comma separated, with " ..." when there are more'''
	names = list(names)
	return ", ".join(names[:LISTED]) + (" ..." if len(names) > LISTED else "")


def reference_snps(database):
	'''{reference genome name: [(snp_id, sequence, position, ancestral base, derived base), ...]}
		from an open database connection, for the canSNPs the scheme types. The alternative
		canSNPs of a node are left out, so that the bases a reference file must carry, and
		the mismatch fraction counted over them, are the same as they were in 2.3.0'''
	query = '''SELECT snp_references.genome, snp_annotation.snp_id, snp_annotation.sequence,
					snp_annotation.sequence_position, snp_annotation.ancestral_base, snp_annotation.derived_base
					FROM snp_annotation
					JOIN snp_references ON (snp_references.id = snp_annotation.genome_i)
					WHERE snp_annotation.typed = 1'''
	snps = {}
	for genome, snp_id, sequence, position, ancestral, derived in database.query(query).fetchall():
		if position is None or sequence is None:
			continue
		snps.setdefault(genome, []).append((snp_id, str(sequence), int(position), str(ancestral).upper(), str(derived).upper()))
	return snps


class FastaFileError(ValueError):
	'''a fasta file that cannot be read as text; the message names the file'''
	pass


def open_fasta(path):
	'''a fasta file opened as text. A UTF-8 byte order mark is dropped (it hid the first
		header), bytes that are not UTF-8 are replaced, as they can only be in headers and
		comments of a usable file, and a gzip-compressed file is refused by name: the
		aligners cannot read it and the decoder only reported a byte it could not decode'''
	with open(path, "rb") as fin:
		magic = fin.read(2)
	if magic == b"\x1f\x8b":
		raise FastaFileError("{path} is gzip-compressed, decompress it first".format(path=path))
	return open(path, encoding="utf-8-sig", errors="replace")


def read_sequences(path):
	'''[(name, sequence), ...] of a fasta file in file order: the name is the first word of the
		header line, the sequence is upper case and without whitespace'''
	records = []
	with open_fasta(path) as fin:
		for line in fin:
			if line.startswith(">"):
				words = line[1:].split()
				records.append((words[0] if words else "", []))
			elif records:
				records[-1][1].append("".join(line.split()).upper())
	return [(name, "".join(parts)) for name, parts in records]


def missing_sequences(path, name, sequences, snps):
	'''None when every sequence named in snps is in the file, otherwise the reason'''
	missing = sorted({sequence for _snp, sequence, _pos, _anc, _der in snps if sequence not in sequences})
	if not missing:
		return None
	return ("reference {name} ({path}) has no sequence named {missing}, on which the database places SNPs; "
			"the sequences of the file are {present}. Sequences are identified by the first word of their "
			"header line. A file with RefSeq names (NC_, NZ_) for a scheme built on GenBank accessions must "
			"be renamed to the GenBank accessions, or fetched with CanSNPer2-download").format(
			name=name, path=path, missing=listing(missing), present=listing(sequences) or "none")


def too_short(path, name, sequences, snps):
	'''None when each sequence holds its largest SNP position, otherwise the reason'''
	last = {}
	for _snp, sequence, position, _anc, _der in snps:
		last[sequence] = max(position, last.get(sequence, 0))
	for sequence in sorted(last):
		length = len(sequences[sequence])
		if length < last[sequence]:
			return ("sequence {sequence} of reference {name} ({path}) has {length} bases, but the database places "
					"SNPs on it up to position {last}: the file is truncated or is not the genome the database was "
					"built on; delete it and run CanSNPer2-download again").format(
					sequence=sequence, name=name, path=path, length=length, last=last[sequence])
	return None


def base_mismatches(sequences, snps):
	'''the SNPs whose position in the reference holds neither the ancestral nor the derived base:
		[(snp_id, sequence, position, reference base, ancestral, derived), ...]; sequences is
		{name: sequence}, a SNP on a missing sequence or beyond the end is left to the other checks'''
	mismatches = []
	for snp_id, sequence, position, ancestral, derived in snps:
		seq = sequences.get(sequence)
		if seq is not None and 0 < position <= len(seq):
			base = seq[position - 1]
			if base != ancestral and base != derived:
				mismatches.append((snp_id, sequence, position, base, ancestral, derived))
	return mismatches


def check_bases(path, name, sequences, snps):
	'''None when the reference has the ancestral or the derived base at enough of its SNPs,
		otherwise the reason; mismatches below the limit are logged as a warning. The fraction
		is counted over all SNPs of the reference.'''
	mismatches = base_mismatches(sequences, snps)
	if not mismatches:
		return None
	listed = listing("{snp} at {seq}:{pos} has {base} (ancestral {anc}, derived {der})".format(
			snp=snp, seq=seq, pos=pos, base=base, anc=anc, der=der) for snp, seq, pos, base, anc, der in mismatches)
	share = "{n} of {m} SNPs".format(n=len(mismatches), m=len(snps))
	if len(mismatches) >= MISMATCH_FRACTION * len(snps):
		return ("reference {name} ({path}) has neither the ancestral nor the derived base at {share} ({listed}): "
				"the reference base at a SNP position must be one of them. The file is probably not the genome, "
				"or not the version of its sequences, that the database was built on; delete it and run "
				"CanSNPer2-download again").format(name=name, path=path, share=share, listed=listed)
	logger.warning("reference {name} has neither the ancestral nor the derived base at {share}: {listed}; the scheme may give a wrong base for these SNPs".format(
			name=name, share=share, listed=listed))
	return None


def check_reference(path, name, snps=None):
	'''None when the reference file at path is usable for its SNPs, otherwise the reason.
		snps: [(snp_id, sequence, position, ancestral, derived), ...] of this reference.
		Checked in this order: unique and non-empty sequence names, every sequence named in
		snps present, each sequence at least as long as its largest SNP position, and the
		ancestral or derived base at enough SNPs.'''
	try:
		records = read_sequences(path)
	except FastaFileError as e:
		return "reference {name}: {problem}".format(name=name, problem=e)
	problem = duplicate_names([seq_name for seq_name, _seq in records])
	if problem:
		return "reference {name} ({path}) has {problem}".format(name=name, path=path, problem=problem)
	if not snps:
		return None
	sequences = dict(records)
	return (missing_sequences(path, name, sequences, snps)
			or too_short(path, name, sequences, snps)
			or check_bases(path, name, sequences, snps))


def sequence_names(path):
	'''the name of every sequence in a fasta file, in file order: the first word of the header
		line, which is how minimap2, MUMmer and the result files identify a sequence'''
	names = []
	with open_fasta(path) as fin:
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
		problems.append("more than one sequence named {names}".format(names=listing(duplicates)))
	if empty:
		problems.append("{n} sequence(s) without a name".format(n=empty))
	if not problems:
		return None
	return "{problems}; sequence names (the first word of each header line) must be set and unique".format(problems=" and ".join(problems))
