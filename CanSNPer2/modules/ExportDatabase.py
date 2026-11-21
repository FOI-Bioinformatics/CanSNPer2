'''
Export a CanSNPer2 database to the source files that CanSNPer2-database --create reads:

	<name>_tree.txt        one line per node, the tab separated path from the root to the node
	<name>_snp.txt         the SNP annotation, one row per SNP
	<name>_references.txt  the reference genomes

These are the files kept in CanSNPer2-data/source. Creating a database from an export gives
the same tree, SNPs and references. The children of a node are written in the order of their
node ids, so a database built from source files also gets the same node ids back.
'''

import logging
import os
import sqlite3

from CanSNPer2.modules.DatabaseConnection import GENOME_SOURCE_COLUMNS

logger = logging.getLogger(__name__)

## the SNP header of CanSNPer2-data (validate.py SNP_HEADER), so an export can be added there
## as it is. The database keeps no strain per SNP: the strain of the SNP's reference genome is
## written, as CanSNPCreator does. The date column is not written, it is not part of that
## header; a database built from the export gets the build date.
SNP_COLUMNS = ("snp_id", "strain", "reference", "genome", "sequence", "sequence_position", "derived_base",
	"ancestral_base")


class ExportError(Exception):
	'''the database cannot be written as source files'''


def text(value):
	return "" if value is None else str(value)


def tree_paths(conn):
	'''Return the path from the root to every node, in the order the tree file lists them.
		The source format writes each node with its full path, so every node needs exactly one
		parent and must be reachable from the root, the node that is linked to itself.'''
	names = dict(conn.execute("SELECT id,name FROM nodes"))
	parent_of, children, roots = {}, {}, []
	for child, parent in conn.execute("SELECT child,parent FROM tree ORDER BY child"):
		if child == parent:
			roots.append(child)
			continue
		if child in parent_of:
			raise ExportError("node {node} has more than one parent ({a}, {b}), the tree file can only "
				"describe a tree".format(node=names.get(child, child), a=names.get(parent_of[child], parent_of[child]),
				b=names.get(parent, parent)))
		parent_of[child] = parent
		children.setdefault(parent, []).append(child)
	if len(roots) != 1:
		raise ExportError("the database has {n} root nodes (a node linked to itself), the tree file needs "
			"exactly one".format(n=len(roots)))
	missing = sorted((set(parent_of) | set(children) | set(roots)) - set(names))
	if missing:
		raise ExportError("the tree links {n} node id(s) that are not in the nodes table: {ids}".format(
			n=len(missing), ids=", ".join(map(str, missing[:10]))))
	paths = []
	stack = [(roots[0], [names[roots[0]]])]
	seen = set()
	while stack:
		node, path = stack.pop()
		if node in seen:
			raise ExportError("the tree has a cycle through node {node}".format(node=names[node]))
		seen.add(node)
		paths.append(path)
		for child in reversed(children.get(node, [])):
			stack.append((child, path + [names[child]]))
	unreachable = sorted(names[n] for n in set(parent_of) - seen)
	if unreachable:
		raise ExportError("{n} node(s) are not connected to the root {root}: {names}".format(
			n=len(unreachable), root=names[roots[0]], names=", ".join(unreachable[:10])))
	return paths


def export_source(database, outdir, name=None):
	'''Write the tree, SNP and reference source files of database into outdir and return
		their paths. name is the file prefix, by default the database file name.'''
	if name is None:
		name = os.path.basename(database).rsplit(".", 1)[0]
	conn = sqlite3.connect(database)
	try:
		paths = tree_paths(conn)
		references = conn.execute("SELECT {cols} FROM snp_references ORDER BY id".format(
			cols=",".join(GENOME_SOURCE_COLUMNS))).fetchall()
		snps = conn.execute("""SELECT a.snp_id, g.strain, a.reference, g.genome, a.sequence, a.sequence_position,
				a.derived_base, a.ancestral_base
			FROM snp_annotation a LEFT JOIN snp_references g ON g.id = a.genome_i
			ORDER BY a.node_id""").fetchall()
	finally:
		conn.close()
	without_genome = [row[0] for row in snps if row[3] is None]
	if without_genome:
		raise ExportError("{n} SNP(s) point to a reference genome that is not in the snp_references table: "
			"{snps}".format(n=len(without_genome), snps=", ".join(map(str, without_genome[:10]))))

	os.makedirs(outdir, exist_ok=True)
	files = {key: os.path.join(outdir, "{name}_{key}.txt".format(name=name, key=key))
		for key in ("tree", "snp", "references")}
	with open(files["tree"], "w") as out:
		for path in paths:
			print(*path, sep="\t", file=out)
	with open(files["snp"], "w") as out:
		print(*SNP_COLUMNS, sep="\t", file=out)
		for row in snps:
			print(*[text(value) for value in row], sep="\t", file=out)
	with open(files["references"], "w") as out:
		print(*GENOME_SOURCE_COLUMNS, sep="\t", file=out)
		for row in references:
			print(*[text(value) for value in row], sep="\t", file=out)
	logger.info("Exported {nodes} nodes, {snps} SNPs and {refs} reference genomes to {outdir}".format(
		nodes=len(paths), snps=len(snps), refs=len(references), outdir=outdir))
	return files
