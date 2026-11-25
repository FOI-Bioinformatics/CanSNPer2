'''Shared test fixtures: small CanSNP databases and source file triples.

A SNP is stored as a reference sequence (its GenBank accession) and a 1-based position
on that sequence. make_db writes that table directly, without the loaders, for tests that
only need a database to read from; write_sources writes the tree, SNP and references
files that CanSNPer2-database --create reads.

pytest puts this folder on sys.path (the tests folder has no __init__.py), so the test
files import it as "from helpers import make_db".'''
import os
import sqlite3

## the snp_annotation table of the current format, as InitiateDatabase creates it: a node has
## one row per canSNP, and the canSNPs the scheme types carry typed = 1. The defaults of rank
## and typed let a test write a single typed canSNP of a node with the columns of 2.3.0
SNP_ANNOTATION_TABLE = """CREATE TABLE snp_annotation (node_id integer, snp_id varchar(6),
	sequence varchar(40), sequence_position integer, ancestral_base varchar(1), derived_base varchar(1),
	reference varchar(20), date datetime, genome_i integer, rank integer NOT NULL DEFAULT 1,
	score text, typed integer NOT NULL DEFAULT 1, UNIQUE (node_id, rank),
	UNIQUE (genome_i, sequence, sequence_position));"""

## the table of 2.3.0: one canSNP per node, without the rank, score and typed columns
SNP_TABLE_2_3 = """CREATE TABLE snp_annotation (node_id integer PRIMARY KEY, snp_id varchar(6),
	sequence varchar(40), sequence_position integer, ancestral_base varchar(1), derived_base varchar(1),
	reference varchar(20), date datetime, genome_i integer);"""

## the table of the format before 2.3.0: one position along the joined reference sequences
JOINED_SNP_ANNOTATION_TABLE = """CREATE TABLE snp_annotation (node_id integer PRIMARY KEY, snp_id varchar(6),
	position integer, ancestral_base varchar(1), derived_base varchar(1), reference varchar(20),
	date datetime, genome_i integer);"""

REFERENCES_TABLE = """CREATE TABLE snp_references (id integer PRIMARY KEY, genome varchar(6),
	strain varchar(200), genbank_id varchar(20), refseq_id varchar(20), assembly_name varchar(200));"""

## the genome_i of the single reference genome that make_db writes
GENOME_I = 5

## the header of a canSNP file, which is what --create and --import read
CANSNP_HEADER = ("node", "rank", "reference", "genome", "sequence", "sequence_position",
	"derived_base", "ancestral_base", "score", "typed")
## the SNP file of 2.3.0, still accepted by --annotation
SNP_HEADER = ("snp_id", "strain", "reference", "genome", "sequence", "sequence_position",
	"derived_base", "ancestral_base")
REFERENCES_HEADER = ("genome", "strain", "genbank_id", "refseq_id", "assembly_name")


def make_db(path, rows, genome="REF", joined=False, old_2_3=False, extra=()):
	'''Write a CanSNP database with one reference genome and one node per SNP, and return path.
		rows: (snp_id, sequence, position, ancestral_base, derived_base); node ids are 1, 2, ...
		in row order, and the tree table is left empty for the caller to fill. Each becomes the
		typed canSNP of rank 1 of its node.
		extra: (snp_id, rank, sequence, position, ancestral_base, derived_base) canSNPs that are
		not typed, of nodes that rows already created.
		joined=True writes the table of the format before 2.3.0 (a single position column), and
		old_2_3=True the table of 2.3.0 (no rank, score or typed column); the sequence item of
		each row is ignored with joined=True.'''
	path = str(path)
	if os.path.exists(path):
		os.remove(path)
	conn = sqlite3.connect(path)
	try:
		conn.executescript("""
			CREATE TABLE nodes (id integer PRIMARY KEY, name text NOT NULL);
			CREATE TABLE tree (parent integer, child integer, rank_i integer, unique (parent, child, rank_i));
			{snp_table}
			{references_table}
		""".format(snp_table=JOINED_SNP_ANNOTATION_TABLE if joined else (SNP_TABLE_2_3 if old_2_3 else SNP_ANNOTATION_TABLE),
			references_table=REFERENCES_TABLE))
		conn.execute("INSERT INTO snp_references (id, genome, strain) VALUES (?, ?, ?)", (GENOME_I, genome, genome))
		for node_id, (snp_id, sequence, position, ancestral, derived) in enumerate(rows, start=1):
			conn.execute("INSERT INTO nodes (id, name) VALUES (?, ?)", (node_id, snp_id))
			if joined:
				conn.execute("INSERT INTO snp_annotation (node_id, snp_id, position, ancestral_base, derived_base,"
					" reference, date, genome_i) VALUES (?, ?, ?, ?, ?, 'test', '06/11/26', ?)",
					(node_id, snp_id, position, ancestral, derived, GENOME_I))
			else:
				conn.execute("INSERT INTO snp_annotation (node_id, snp_id, sequence, sequence_position, ancestral_base,"
					" derived_base, reference, date, genome_i) VALUES (?, ?, ?, ?, ?, ?, 'test', '06/11/26', ?)",
					(node_id, snp_id, sequence, position, ancestral, derived, GENOME_I))
		node_ids = dict((snp_id, node_id) for node_id, (snp_id, _s, _p, _a, _d) in enumerate(rows, start=1))
		for snp_id, rank, sequence, position, ancestral, derived in extra:
			conn.execute("INSERT INTO snp_annotation (node_id, snp_id, sequence, sequence_position, ancestral_base,"
				" derived_base, reference, date, genome_i, rank, score, typed)"
				" VALUES (?, ?, ?, ?, ?, ?, 'test', '06/11/26', ?, ?, NULL, 0)",
				(node_ids[snp_id], snp_id, sequence, position, ancestral, derived, GENOME_I, rank))
		conn.commit()
	finally:
		conn.close()
	return path


def _line(row):
	return row if isinstance(row, str) else "\t".join(str(value) for value in row)


def write_sources(folder, snp_rows, tree=None, references=None, name="scheme", header=CANSNP_HEADER):
	'''Write the tree, canSNP and references files of a scheme into folder and return their paths
		as {"tree", "snp", "references"}; the "snp" key is the canSNP file.
		snp_rows: tab separated lines or tuples in the order of header, by default
		(node, rank, reference, genome, sequence, sequence_position, derived_base, ancestral_base,
		score, typed).
		tree: the tree file text; by default the first SNP is the root and every other SNP a
		child of it. references: the references file text; by default one row per genome of
		the SNP rows, with accessions of the form GCA_00000000N.1.
		header: the SNP header, a tuple of column names or a line.'''
	folder = str(folder)
	os.makedirs(folder, exist_ok=True)
	lines = [_line(row) for row in snp_rows]
	columns = header if isinstance(header, str) else "\t".join(header)
	fields = columns.split("\t")
	parsed = [dict(zip(fields, line.split("\t"))) for line in lines]
	if tree is None:
		## the node column of a canSNP file, the snp_id of a SNP file, or whatever the first
		## column is named in a header a test made up
		## a node may hold several canSNPs, so each node appears in the tree once
		names = []
		for row in parsed:
			name = row.get("node", row.get("snp_id", row[fields[0]]))
			if name not in names:
				names.append(name)
		tree = "".join([names[0] + "\n"] + ["{0}\t{1}\n".format(names[0], n) for n in names[1:]]) if names else ""
	if references is None:
		genomes = []
		for row in parsed:
			if row.get("genome") and row["genome"] not in genomes:
				genomes.append(row["genome"])
		references = "\t".join(REFERENCES_HEADER) + "\n" + "".join(
			"{g}\t{g} strain\tGCA_{n:09d}.1\tGCF_{n:09d}.1\tASM{n}v1\n".format(g=genome, n=i)
			for i, genome in enumerate(genomes, start=1))
	files = {}
	## the canSNP file is named <name>_canSNPs.txt, as CanSNPer2-data names it; a header of the
	## SNP file of 2.3.0 is written under the old name, so that both formats can be built
	snp_name = "snp" if "snp_id" in fields else "canSNPs"
	for key, text in (("tree", tree), ("snp", columns + "\n" + "".join(line + "\n" for line in lines)),
			("references", references)):
		files[key] = os.path.join(folder, "{0}_{1}.txt".format(name, snp_name if key == "snp" else key))
		with open(files[key], "w") as out:
			out.write(text)
	return files


def _revcomp(seq):
	return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def simulate_reads(genome, depth, seed=1, length=100, insert=300, paired=True, error_rate=0.0):
	'''error-free (or uniformly erroneous) reads of a genome sequence, as FASTQ records;
		returns (reads1, reads2) for paired-end reads (reads2 is None for single-end). Both
		strands are sampled; a pair spans insert bases, the second read on the other strand.'''
	import random
	rng = random.Random(seed)
	span = insert if paired else length
	n = max(1, int(depth * len(genome) / (length * (2 if paired else 1))))
	def noisy(seq):
		if not error_rate:
			return seq
		return "".join(rng.choice("ACGT".replace(b, "")) if rng.random() < error_rate else b for b in seq)
	r1, r2 = [], []
	for i in range(n):
		start = rng.randint(0, len(genome) - span)
		fragment = genome[start:start + span]
		if rng.random() < 0.5:
			fragment = _revcomp(fragment)
		first = noisy(fragment[:length])
		r1.append("@s{s}r{i}/1\n{seq}\n+\n{q}\n".format(s=seed, i=i, seq=first, q="I" * len(first)))
		if paired:
			second = noisy(_revcomp(fragment)[:length])
			r2.append("@s{s}r{i}/2\n{seq}\n+\n{q}\n".format(s=seed, i=i, seq=second, q="I" * len(second)))
	return "".join(r1), ("".join(r2) if paired else None)
