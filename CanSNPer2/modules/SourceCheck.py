'''
Check the CanSNPer2 source files before a database is built from them.

The loaders read the tree, SNP and reference files row by row and used to accept
almost anything: a position "abc" or bases "a"/"a" were stored, a SNP on an unknown
genome or a repeated snp_id was dropped with only a warning, and a node given two
parents became two nodes. check_sources() reads the files the way the loaders do and
returns every problem it finds, so --create, --import and the scheme download can
refuse the files before a database is started.

The rules follow CanSNPer2-data/validate.py, so a scheme that passes there builds
here. Format rules that only concern the published repository (line endings,
the RefSeq accession format, a README, a SNP on every node) are left to validate.py.
The GenBank assembly accession of every reference is checked here, as CanSNPer2-download
needs it to fetch the genome.

check_modification() checks the files given to modify an existing database against the
nodes and the reference genomes that database already holds.
'''

import gzip
import re
import sqlite3
from collections import Counter, OrderedDict

from CanSNPer2.modules.DatabaseConnection import ANNOTATION_SOURCE_COLUMNS, GENOME_SOURCE_COLUMNS, JOINED_FORMAT_HINT

BASES = ("A", "C", "G", "T")
## ASCII digits only: str.isdigit() also accepts characters such as superscripts
POSITION_RE = re.compile(r"^[0-9]+$")
## trailing free-text columns that rows of the SNP file may leave out (issue #22)
COMMENT_COLUMNS = ("comment", "comments", "note", "notes", "free text", "free_text")
## a GenBank (GCA_) or RefSeq (GCF_) assembly accession with its version
ACCESSION_RE = re.compile(r"^GC[AF]_[0-9]{9}\.[0-9]+$")


def _open(path):
	if str(path).endswith(".gz"):
		return gzip.open(path, "rt")
	return open(path)


def check_tree(path, problems):
	'''Return the node names of a CanSNPer tree file (one path from the root per line,
		tab or ; separated), adding a problem for every line that does not fit one tree'''
	nodes = set()
	parents = OrderedDict()
	root = None
	with _open(path) as tree:
		for number, line in enumerate(tree, start=1):
			if line.strip() == "":
				continue
			## the same reading as ReadTaxonomyCanSNPer.parse_taxonomy
			fields = line.strip().replace("\t", ";").split(";")
			empty = [i + 1 for i, field in enumerate(fields) if field.strip() == ""]
			if empty:
				where = "at the end" if empty == [len(fields)] else "between populated columns"
				problems.append("{0}:{1}: empty column {2} {3}".format(path, number, empty[0], where))
				fields = [field for field in fields if field.strip() != ""]
				if not fields:
					continue
			if root is None:
				root = fields[0]
			elif fields[0] != root:
				problems.append("{0}:{1}: every line must start at the root '{2}' given on the "
					"first line, this line starts at '{3}'".format(path, number, root, fields[0]))
				continue
			for i, node in enumerate(fields):
				nodes.add(node)
				if i == 0:
					continue
				if node == root:
					problems.append("{0}:{1}: the root '{2}' is given a parent".format(path, number, root))
					continue
				parents.setdefault(node, OrderedDict())[fields[i - 1]] = number
	if root is None:
		problems.append("{0}: the tree file is empty".format(path))
	for node, node_parents in parents.items():
		if len(node_parents) > 1:
			problems.append("{0}: node '{1}' is given more than one parent: {2}".format(
				path, node, ", ".join(sorted(node_parents))))
	return nodes


def _header(handle):
	'''the header row as the loaders read it'''
	return handle.readline().lstrip("#").strip().split("\t")


def check_references(path, problems, known=None, database=None, known_accessions=None):
	'''Return the names a SNP may use for its genome (the genome column and both
		accessions, as the loaders accept all three), mapped to the genome name.
		known: the genome names a database already holds, which a new row must not repeat;
		known_accessions: {genbank_id: genome} of that database, an accession a new row must
		not give for another genome; database: the name of that database in the messages'''
	names = {}
	genomes = Counter()
	accessions = {}
	with open(path) as references:
		headers = _header(references)
		missing = [column for column in GENOME_SOURCE_COLUMNS if column not in headers]
		if missing:
			problems.append("{0}:1: the header has no column(s) {1}, it must name {2}".format(
				path, ", ".join(missing), ", ".join(GENOME_SOURCE_COLUMNS)))
			return names
		for number, line in enumerate(references, start=2):
			if line.strip() == "":
				continue
			data = line.rstrip("\r\n").split("\t")
			if len(data) != len(headers):
				problems.append("{0}:{1}: expected {2} columns, found {3}".format(
					path, number, len(headers), len(data)))
				continue
			row = dict(zip(headers, data))
			genome = row["genome"]
			genomes[genome] += 1
			if not ACCESSION_RE.match(row["genbank_id"]):
				problems.append("{0}:{1}: genbank_id '{2}' is not an assembly accession (GCA_ or GCF_, nine digits "
					"and a version, such as GCA_000008985.1); CanSNPer2-download needs it to fetch the genome".format(
					path, number, row["genbank_id"]))
			accession = row["genbank_id"]
			if accession in accessions:
				problems.append("{0}:{1}: genbank_id '{2}' of genome '{3}' is already used by genome '{4}' on line {5}; "
					"a source file is written per accession, so two genomes cannot share one".format(
					path, number, accession, genome, accessions[accession][0], accessions[accession][1]))
			else:
				accessions[accession] = (genome, number)
			if known_accessions and accession in known_accessions and known_accessions[accession] != genome:
				problems.append("{0}:{1}: genbank_id '{2}' of genome '{3}' is already used by genome '{4}' in {5}; "
					"a source file is written per accession, so two genomes cannot share one".format(
					path, number, accession, genome, known_accessions[accession], database))
			if known is not None and genome in known:
				problems.append("{0}:{1}: genome '{2}' is already in {3}".format(path, number, genome, database))
			for column in ("refseq_id", "genbank_id", "genome"):
				if row[column]:
					names[row[column]] = genome
	for genome, count in genomes.items():
		if count > 1:
			problems.append("{0}: genome '{1}' is listed {2} times".format(path, genome, count))
	return names


def check_snps(path, problems, tree_nodes=None, genomes=None, references_path=None, tree_name="the tree file",
		existing=None):
	'''Check every SNP row; tree_nodes and genomes are skipped when not given.
		existing: (database name, snp_ids, {(genome, sequence, position)}) of the SNPs a database
		already holds, which a new row must not repeat'''
	snp_ids = Counter()
	positions = Counter()
	with open(path) as snps:
		headers = _header(snps)
		if "position" in headers and "sequence_position" not in headers:
			problems.append("{0}:1: this SNP file cannot be used. {1}".format(path, JOINED_FORMAT_HINT))
			return
		missing = [column for column in ANNOTATION_SOURCE_COLUMNS if column not in headers]
		if missing:
			problems.append("{0}:1: the header has no column(s) {1}, it must name {2}".format(
				path, ", ".join(missing), ", ".join(ANNOTATION_SOURCE_COLUMNS)))
			return
		optional = 0
		for head in reversed(headers):
			if head.strip().lower().lstrip("#") in COMMENT_COLUMNS:
				optional += 1
			else:
				break
		required = len(headers) - optional
		for number, line in enumerate(snps, start=2):
			if line.strip() == "" or line.strip().startswith("#"):
				continue
			data = line.rstrip("\r\n").split("\t")
			if not required <= len(data) <= len(headers):
				expected = str(required) if optional == 0 else "{0} to {1}".format(required, len(headers))
				problems.append("{0}:{1}: expected {2} columns, found {3}".format(path, number, expected, len(data)))
				continue
			row = dict(zip(headers, data))
			snp_id, genome = row["snp_id"], row["genome"]
			sequence, position = row["sequence"], row["sequence_position"]
			derived, ancestral = row["derived_base"], row["ancestral_base"]
			snp_ids[snp_id] += 1
			if tree_nodes is not None and snp_id not in tree_nodes:
				problems.append("{0}:{1}: snp_id '{2}' does not appear in {3}".format(path, number, snp_id, tree_name))
			if existing is not None and snp_id in existing[1]:
				problems.append("{0}:{1}: snp_id '{2}' already has a SNP in {3}".format(path, number, snp_id, existing[0]))
			if genomes is not None:
				if genome not in genomes:
					problems.append("{0}:{1}: genome '{2}' is not listed in {3}".format(
						path, number, genome, references_path))
				else:
					genome = genomes[genome]
			if sequence.strip() == "":
				problems.append("{0}:{1}: sequence is empty, it must name a sequence of the reference genome".format(
					path, number))
			if not POSITION_RE.match(position) or int(position) < 1:
				problems.append("{0}:{1}: sequence_position '{2}' is not a positive integer".format(path, number, position))
			elif sequence.strip() != "":
				positions[(genome, sequence, int(position))] += 1
				if existing is not None and (genome, sequence, int(position)) in existing[2]:
					problems.append("{0}:{1}: position {2} on sequence '{3}' of '{4}' is already annotated in {5}".format(
						path, number, position, sequence, genome, existing[0]))
			for column, base in (("derived_base", derived), ("ancestral_base", ancestral)):
				if base not in BASES:
					problems.append("{0}:{1}: {2} '{3}' is not one of A, C, G, T".format(path, number, column, base))
			if derived == ancestral and derived in BASES:
				problems.append("{0}:{1}: derived and ancestral base are both '{2}'".format(path, number, derived))
	for snp_id, count in snp_ids.items():
		if count > 1:
			problems.append("{0}: snp_id '{1}' is listed {2} times".format(path, snp_id, count))
	for (genome, sequence, position), count in positions.items():
		if count > 1:
			problems.append("{0}: position {1} on sequence '{2}' of '{3}' is annotated {4} times".format(
				path, position, sequence, genome, count))


def check_sources(tree=None, annotation=None, references=None):
	'''Return the problems of a set of source files, an empty list when they can be built'''
	problems = []
	tree_nodes = check_tree(tree, problems) if tree else None
	genomes = check_references(references, problems) if references else None
	if annotation:
		if not references:
			problems.append("{0}: the SNP file needs the references file (--references) that lists "
				"its genomes".format(annotation))
		check_snps(annotation, problems, tree_nodes=tree_nodes, genomes=genomes, references_path=references)
	return problems


def check_modification(database, annotation=None, references=None, name=None):
	'''Return the problems of a SNP file and a references file that are to be added to an
		existing database, an empty list when they can be added. The SNPs are checked against
		the nodes of the tree in the database and against the reference genomes of the database
		and of the references file; a genome, a snp_id or a position the database already holds
		is a problem, and so is a genbank_id the database holds for another genome. name: the database as it is named in the messages (default database)'''
	name = name or database
	problems = []
	conn = sqlite3.connect(database)
	try:
		tree_nodes = set(row[0] for row in conn.execute(
			"SELECT name FROM nodes WHERE id IN (SELECT child FROM tree UNION SELECT parent FROM tree)"))
		genomes = {}
		accessions = {}
		for genome, genbank_id, refseq_id in conn.execute("SELECT genome, genbank_id, refseq_id FROM snp_references"):
			for key in (refseq_id, genbank_id, genome):
				if key:
					genomes[key] = genome
			if genbank_id:
				accessions[genbank_id] = genome
		snp_ids = set(row[0] for row in conn.execute("SELECT snp_id FROM snp_annotation"))
		positions = set((genome, str(sequence), int(position)) for genome, sequence, position in conn.execute(
			"""SELECT snp_references.genome, snp_annotation.sequence, snp_annotation.sequence_position
				FROM snp_annotation JOIN snp_references ON (snp_references.id = snp_annotation.genome_i)""")
			if sequence is not None and str(position).isdigit())
	finally:
		conn.close()
	known = set(genomes.values())
	if references:
		genomes.update(check_references(references, problems, known=known, database=name,
			known_accessions=accessions))
	if annotation:
		check_snps(annotation, problems, tree_nodes=tree_nodes, genomes=genomes,
			references_path=name if not references else "{0} or {1}".format(name, references),
			tree_name="the tree of {0}".format(name), existing=(name, snp_ids, positions))
	return problems
