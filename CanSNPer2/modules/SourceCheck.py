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
accession formats, a README, a SNP on every node) are left to validate.py.
'''

import gzip
import re
from collections import Counter, OrderedDict

from CanSNPer2.modules.DatabaseConnection import ANNOTATION_SOURCE_COLUMNS, GENOME_SOURCE_COLUMNS

BASES = ("A", "C", "G", "T")
## ASCII digits only: str.isdigit() also accepts characters such as superscripts
POSITION_RE = re.compile(r"^[0-9]+$")
## trailing free-text columns that rows of the SNP file may leave out (issue #22)
COMMENT_COLUMNS = ("comment", "comments", "note", "notes", "free text", "free_text")


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


def check_references(path, problems):
	'''Return the names a SNP may use for its genome (the genome column and both
		accessions, as the loaders accept all three), mapped to the genome name'''
	names = {}
	genomes = Counter()
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
			for column in ("refseq_id", "genbank_id", "genome"):
				if row[column]:
					names[row[column]] = genome
	for genome, count in genomes.items():
		if count > 1:
			problems.append("{0}: genome '{1}' is listed {2} times".format(path, genome, count))
	return names


def check_snps(path, problems, tree_nodes=None, genomes=None, references_path=None):
	'''Check every SNP row; tree_nodes and genomes are skipped when not given'''
	snp_ids = Counter()
	positions = Counter()
	with open(path) as snps:
		headers = _header(snps)
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
			snp_id, genome, position = row["snp_id"], row["genome"], row["position"]
			derived, ancestral = row["derived_base"], row["ancestral_base"]
			snp_ids[snp_id] += 1
			if tree_nodes is not None and snp_id not in tree_nodes:
				problems.append("{0}:{1}: snp_id '{2}' does not appear in the tree file".format(path, number, snp_id))
			if genomes is not None:
				if genome not in genomes:
					problems.append("{0}:{1}: genome '{2}' is not listed in {3}".format(
						path, number, genome, references_path))
				else:
					genome = genomes[genome]
			if not POSITION_RE.match(position) or int(position) < 1:
				problems.append("{0}:{1}: position '{2}' is not a positive integer".format(path, number, position))
			else:
				positions[(genome, int(position))] += 1
			for column, base in (("derived_base", derived), ("ancestral_base", ancestral)):
				if base not in BASES:
					problems.append("{0}:{1}: {2} '{3}' is not one of A, C, G, T".format(path, number, column, base))
			if derived == ancestral and derived in BASES:
				problems.append("{0}:{1}: derived and ancestral base are both '{2}'".format(path, number, derived))
	for snp_id, count in snp_ids.items():
		if count > 1:
			problems.append("{0}: snp_id '{1}' is listed {2} times".format(path, snp_id, count))
	for (genome, position), count in positions.items():
		if count > 1:
			problems.append("{0}: position {1} on '{2}' is annotated {3} times".format(path, position, genome, count))


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
