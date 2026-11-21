'''Shared test fixtures: small CanSNP databases and source file triples.

A SNP is stored as a reference sequence (its GenBank accession) and a 1-based position
on that sequence. make_db writes that table directly, without the loaders, for tests that
only need a database to read from; write_sources writes the tree, SNP and references
files that CanSNPer2-database --create reads.

pytest puts this folder on sys.path (the tests folder has no __init__.py), so the test
files import it as "from helpers import make_db".'''
import os
import sqlite3

## the snp_annotation table of the current format, as InitiateDatabase creates it
SNP_ANNOTATION_TABLE = """CREATE TABLE snp_annotation (node_id integer PRIMARY KEY, snp_id varchar(6),
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

SNP_HEADER = ("snp_id", "strain", "reference", "genome", "sequence", "sequence_position",
	"derived_base", "ancestral_base")
REFERENCES_HEADER = ("genome", "strain", "genbank_id", "refseq_id", "assembly_name")


def make_db(path, rows, genome="REF", joined=False):
	'''Write a CanSNP database with one reference genome and one node per SNP, and return path.
		rows: (snp_id, sequence, position, ancestral_base, derived_base); node ids are 1, 2, ...
		in row order, and the tree table is left empty for the caller to fill.
		joined=True writes the table of the format before 2.3.0 (a single position column),
		the sequence item of each row is then ignored.'''
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
		""".format(snp_table=JOINED_SNP_ANNOTATION_TABLE if joined else SNP_ANNOTATION_TABLE,
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
		conn.commit()
	finally:
		conn.close()
	return path


def _line(row):
	return row if isinstance(row, str) else "\t".join(str(value) for value in row)


def write_sources(folder, snp_rows, tree=None, references=None, name="scheme", header=SNP_HEADER):
	'''Write the tree, SNP and references files of a scheme into folder and return their paths
		as {"tree", "snp", "references"}.
		snp_rows: tab separated lines or tuples in the order of header, by default
		(snp_id, strain, reference, genome, sequence, sequence_position, derived_base, ancestral_base).
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
		names = [row["snp_id"] for row in parsed]
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
	for key, text in (("tree", tree), ("snp", columns + "\n" + "".join(line + "\n" for line in lines)),
			("references", references)):
		files[key] = os.path.join(folder, "{0}_{1}.txt".format(name, key))
		with open(files[key], "w") as out:
			out.write(text)
	return files
