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

import re
import sqlite3
from collections import Counter, OrderedDict

from CanSNPer2.modules.DatabaseConnection import (CANSNP_SOURCE_COLUMNS, SNP_SOURCE_COLUMNS, GENOME_SOURCE_COLUMNS,
	JOINED_FORMAT_HINT, TYPED_VALUES, SourceFileError, cansnp_row, open_source)

BASES = ("A", "C", "G", "T")
## ASCII digits only: str.isdigit() also accepts characters such as superscripts
POSITION_RE = re.compile(r"^[0-9]+$")
## trailing free-text columns that rows of the SNP file may leave out (issue #22)
COMMENT_COLUMNS = ("comment", "comments", "note", "notes", "free text", "free_text")
## a GenBank (GCA_) or RefSeq (GCF_) assembly accession with its version
ACCESSION_RE = re.compile(r"^GC[AF]_[0-9]{9}\.[0-9]+$")
## the largest position that fits the INTEGER column of the database (a larger value is stored as a REAL)
MAX_POSITION = 2**31 - 1
## a genome name is used in file names: letters, digits, dot, underscore and hyphen, no leading dot
GENOME_NAME_RE = re.compile(r"^[A-Za-z0-9_-][A-Za-z0-9._-]*$")
## a rank counts from 1 within its node, without leading zeros
RANK_RE = re.compile(r"^[1-9][0-9]*$")
## a score as CanSNPCreator writes it; an empty score is allowed as well, a canSNP taken from
## a publication was not scored
SCORE_RE = re.compile(r"^-?[0-9]+(\.[0-9]+)?$")
## CanSNPer2 calls one canSNP per node; a later release that calls several raises this
TYPED_PER_NODE = 1


def _open(path):
	return open_source(path)


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
			## every name is stripped, as the loader does, so "A " and "A" are one node
			fields = [field.strip() for field in line.strip().replace("\t", ";").split(";")]
			empty = [i + 1 for i, field in enumerate(fields) if field == ""]
			if empty:
				where = "at the end" if empty == [len(fields)] else "between populated columns"
				problems.append("{0}:{1}: empty column {2} {3}".format(path, number, empty[0], where))
				fields = [field for field in fields if field != ""]
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
	'''the header row as the loaders read it, None when the file is empty'''
	line = handle.readline()
	if line.strip() == "":
		return None
	return line.lstrip("#").strip().split("\t")


def split_row(line, columns):
	'''The fields of a data row. A row is read as the header is: trailing tabs after the last
		column are dropped, so a file that ends every line with a tab is read alike. Only empty fields
		beyond the number of header columns are dropped, an empty last value (such as a reference
		without an assembly name) is kept.'''
	data = [value.strip() for value in line.rstrip("\r\n").split("\t")]
	while len(data) > columns and data[-1] == "":
		data.pop()
	return data


def duplicate_columns(headers):
	'''the column names that the header gives more than once'''
	return [name for name, count in Counter(headers).items() if count > 1]


def register_genome_names(owners, row, source):
	'''Record the names a SNP may use for a genome (the genome column and both accessions) in
		owners, {name: genome}, and raise SourceFileError when a name already belongs to another
		genome: a SNP that names that genome would be attached to the wrong one.
		source: the file or database named in the message'''
	genome = row["genome"]
	for column in ("genome", "genbank_id", "refseq_id"):
		name = row[column]
		if not name:
			continue
		if owners.get(name, genome) != genome:
			raise SourceFileError("{0}: the name '{1}' ({2} of genome '{3}') is already used for genome '{4}'; "
				"a genome, a genbank_id and a refseq_id must not repeat the name of another genome".format(
				source, name, column, genome, owners[name]))
		owners[name] = genome


def check_references(path, problems, known=None, database=None, known_accessions=None, known_names=None):
	'''Return the names a SNP may use for its genome (the genome column and both
		accessions, as the loaders accept all three), mapped to the genome name.
		known: the genome names a database already holds, which a new row must not repeat;
		known_accessions: {genbank_id: genome} of that database, an accession a new row must
		not give for another genome; known_names: {name: genome} of the names that database accepts
		(genome and both accessions), a name a new row must not use for another genome;
		database: the name of that database in the messages'''
	names = {}
	genomes = Counter()
	accessions = {}
	with open_source(path) as references:
		headers = _header(references)
		if headers is None:
			problems.append("{0}: the references file is empty, it needs a header row and one row per genome".format(path))
			return names
		repeated = duplicate_columns(headers)
		if repeated:
			problems.append("{0}:1: the header names the column(s) {1} more than once".format(
				path, ", ".join("'{0}'".format(name) for name in repeated)))
			return names
		missing = [column for column in GENOME_SOURCE_COLUMNS if column not in headers]
		if missing:
			problems.append("{0}:1: the header has no column(s) {1}, it must name {2}".format(
				path, ", ".join(missing), ", ".join(GENOME_SOURCE_COLUMNS)))
			return names
		for number, line in enumerate(references, start=2):
			if line.strip() == "":
				continue
			data = split_row(line, len(headers))
			if len(data) != len(headers):
				problems.append("{0}:{1}: expected {2} columns, found {3}".format(
					path, number, len(headers), len(data)))
				continue
			row = dict(zip(headers, data))
			genome = row["genome"]
			genomes[genome] += 1
			if not GENOME_NAME_RE.match(genome):
				problems.append("{0}:{1}: genome name '{2}' is not allowed, use letters, digits, '.', '_' and '-' "
					"only and do not start with '.' (the name is used in file names)".format(path, number, genome))
			if not ACCESSION_RE.match(row["genbank_id"]):
				problems.append("{0}:{1}: genbank_id '{2}' is not an assembly accession (GCA_ or GCF_, nine digits "
					"and a version, such as GCA_000008985.1); cansnper2 ref get needs it to fetch the genome".format(
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
			for column in ("genome", "genbank_id", "refseq_id"):
				name = row[column]
				if not name:
					continue
				owner = names.get(name)
				if owner is None and known_names:
					owner = known_names.get(name)
				if owner is not None and owner != genome:
					problems.append("{0}:{1}: the name '{2}' ({3} of genome '{4}') is already used for genome '{5}'; "
						"a genome, a genbank_id and a refseq_id must not repeat the name of another genome".format(
						path, number, name, column, genome, owner))
				else:
					names[name] = genome
	for genome, count in genomes.items():
		if count > 1:
			problems.append("{0}: genome '{1}' is listed {2} times".format(path, genome, count))
	return names


def check_snps(path, problems, tree_nodes=None, genomes=None, references_path=None, tree_name="the tree file",
		existing=None):
	'''Check every canSNP row of a canSNP file, or of a SNP file of 2.3.0 and earlier, which is
		read as one typed canSNP of rank 1 per node; tree_nodes and genomes are skipped when not
		given. existing: (database name, {(node, rank)}, {node with a typed canSNP},
		{(genome, sequence, position)}) of the canSNPs a database already holds, which a new row
		must not repeat'''
	ranks = OrderedDict()
	typed_rows = Counter()
	positions = Counter()
	with open_source(path) as snps:
		headers = _header(snps)
		if headers is None:
			problems.append("{0}: the SNP file is empty, it needs a header row and one row per SNP".format(path))
			return
		repeated = duplicate_columns(headers)
		if repeated:
			problems.append("{0}:1: the header names the column(s) {1} more than once".format(
				path, ", ".join("'{0}'".format(name) for name in repeated)))
			return
		if "position" in headers and "sequence_position" not in headers:
			problems.append("{0}:1: this SNP file cannot be used. {1}".format(path, JOINED_FORMAT_HINT))
			return
		## a header with an snp_id column is the SNP file of 2.3.0; anything else is read as a
		## canSNP file, so that a header naming neither is reported against the canSNP columns
		required = SNP_SOURCE_COLUMNS if ("snp_id" in headers and "node" not in headers) else CANSNP_SOURCE_COLUMNS
		missing = [column for column in required if column not in headers]
		if missing:
			problems.append("{0}:1: the header has no column(s) {1}, it must name {2}".format(
				path, ", ".join(missing), ", ".join(CANSNP_SOURCE_COLUMNS)))
			return
		optional = 0
		for head in reversed(headers):
			if head.strip().lower().lstrip("#") in COMMENT_COLUMNS:
				optional += 1
			else:
				break
		required = len(headers) - optional
		rows = 0
		for number, line in enumerate(snps, start=2):
			if line.strip() == "" or line.strip().startswith("#"):
				continue
			data = split_row(line, len(headers))
			if not required <= len(data) <= len(headers):
				expected = str(required) if optional == 0 else "{0} to {1}".format(required, len(headers))
				problems.append("{0}:{1}: expected {2} columns, found {3}".format(path, number, expected, len(data)))
				continue
			rows += 1
			row = cansnp_row(dict(zip(headers, data)))
			node, genome = row["node"], row["genome"]
			sequence, position = row["sequence"], row["sequence_position"]
			derived, ancestral = row["derived_base"], row["ancestral_base"]
			rank, score, typed = row["rank"].strip(), row["score"].strip(), row["typed"].strip().lower()
			if not RANK_RE.match(rank):
				problems.append("{0}:{1}: rank '{2}' is not a positive integer".format(path, number, row["rank"]))
			else:
				ranks.setdefault(node, []).append(int(rank))
			if score != "" and not SCORE_RE.match(score):
				problems.append("{0}:{1}: score '{2}' is not a number".format(path, number, row["score"]))
			if typed not in TYPED_VALUES:
				problems.append("{0}:{1}: typed '{2}' is not yes or no".format(path, number, row["typed"]))
			elif TYPED_VALUES[typed]:
				typed_rows[node] += 1
			if row["reference"].strip() == "":
				problems.append("{0}:{1}: reference is empty, it must name where the canSNP is "
					"published or NA".format(path, number))
			if tree_nodes is not None and node not in tree_nodes:
				problems.append("{0}:{1}: node '{2}' does not appear in {3}".format(path, number, node, tree_name))
			if existing is not None and RANK_RE.match(rank) and (node, int(rank)) in existing[1]:
				problems.append("{0}:{1}: node '{2}' already has a canSNP of rank {3} in {4}".format(
					path, number, node, rank, existing[0]))
			if existing is not None and typed in TYPED_VALUES and TYPED_VALUES[typed] and node in existing[2]:
				problems.append("{0}:{1}: node '{2}' already has a typed canSNP in {3}".format(
					path, number, node, existing[0]))
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
			elif int(position) > MAX_POSITION:
				problems.append("{0}:{1}: sequence_position {2} is larger than {3}, the largest position the "
					"database can hold".format(path, number, position, MAX_POSITION))
			elif sequence.strip() != "":
				positions[(genome, sequence, int(position))] += 1
				if existing is not None and (genome, sequence, int(position)) in existing[3]:
					problems.append("{0}:{1}: position {2} on sequence '{3}' of '{4}' is already annotated in {5}".format(
						path, number, position, sequence, genome, existing[0]))
			for column, base in (("derived_base", derived), ("ancestral_base", ancestral)):
				if base not in BASES:
					problems.append("{0}:{1}: {2} '{3}' is not one of A, C, G, T".format(path, number, column, base))
			if derived == ancestral and derived in BASES:
				problems.append("{0}:{1}: derived and ancestral base are both '{2}'".format(path, number, derived))
	if rows == 0:
		problems.append("{0}: the canSNP file has a header but no canSNP rows, a scheme without "
			"canSNPs cannot be used".format(path))
	for node, listed in ranks.items():
		if listed != list(range(1, len(listed) + 1)):
			problems.append("{0}: the ranks of node '{1}' are {2}, not 1, 2, 3 ... in order".format(
				path, node, ", ".join(str(rank) for rank in listed)))
	for node in ranks:
		count = typed_rows[node]
		if count == 0:
			problems.append("{0}: node '{1}' has no canSNP marked typed yes, CanSNPer2 would not "
				"call it".format(path, node))
		elif count > TYPED_PER_NODE:
			problems.append("{0}: node '{1}' has {2} canSNPs marked typed yes, CanSNPer2 calls one "
				"canSNP per node".format(path, node, count))
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
		and of the references file; a genome, a canSNP of a rank the node already has, a second
		typed canSNP of a node, or a position the database already holds is a problem, and so is a genbank_id the database holds for another genome. name: the database as it is named in the messages (default database)'''
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
		held_ranks = set((row[0], int(row[1])) for row in conn.execute(
			"SELECT snp_id, snp_annotation.rank FROM snp_annotation"))
		typed_nodes = set(row[0] for row in conn.execute(
			"SELECT DISTINCT snp_id FROM snp_annotation WHERE typed = 1"))
		positions = set((genome, str(sequence), int(position)) for genome, sequence, position in conn.execute(
			"""SELECT snp_references.genome, snp_annotation.sequence, snp_annotation.sequence_position
				FROM snp_annotation JOIN snp_references ON (snp_references.id = snp_annotation.genome_i)""")
			if sequence is not None and str(position).isdigit())
	finally:
		conn.close()
	known = set(genomes.values())
	if references:
		genomes.update(check_references(references, problems, known=known, database=name,
			known_accessions=accessions, known_names=dict(genomes)))
	if annotation:
		check_snps(annotation, problems, tree_nodes=tree_nodes, genomes=genomes,
			references_path=name if not references else "{0} or {1}".format(name, references),
			tree_name="the tree of {0}".format(name), existing=(name, held_ranks, typed_nodes, positions))
	return problems
