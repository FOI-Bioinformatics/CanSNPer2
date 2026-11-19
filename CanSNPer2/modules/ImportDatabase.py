'''
Import a CanSNPer2 database into the current database format.

``--import`` reads an existing database (a legacy CanSNPer2 database, one shipped before
the root-node fixes, or one produced by an older release) and writes a new database that
uses the schema and tree structure the current ``--create`` produces: exactly one root, no
orphan nodes, every SNP pointing at a reference genome that exists.

The migration path a user would otherwise have to do by hand:

    old.db  --import-->  new.db

is the same data the current export/rebuild round trip produces, but directly database to
database, so the old ``.dmp`` tables and the legacy databases stay usable. The source is
never modified.
'''

import logging
import os
import sqlite3

from CanSNPer2.modules.DatabaseConnection import check_database, ConnectionError
from CanSNPer2.modules.InitiateDatabase import CreateSNPDatabase
from CanSNPer2.modules.ExportDatabase import tree_paths, ExportError

logger = logging.getLogger(__name__)


class DatabaseImportError(Exception):
	'''the source database could not be imported into the current format'''


# copied verbatim where the column exists in both databases; ids are preserved so node_id,
# genome_i and rank_i keep pointing at the same rows in the new database.
IMPORT_TABLES = ("rank", "nodes", "snp_references", "tree", "snp_annotation", "genomes")


def _columns(conn, table):
	return [row[1] for row in conn.execute("PRAGMA table_info([{0}])".format(table))]


def _find_orphan_roots(conn):
	'''Legacy databases carry an extra root: the old tree reader added a node named "root"
		above the tree's own root and linked it to itself, so the database has two self-linked
		nodes (issues #32, #34, #149). An orphan root is a self-linked node that connects to
		nothing and carries no annotation; it is dropped so the result has a single root.
	'''
	roots = [row[0] for row in conn.execute("SELECT child FROM tree WHERE child=parent")]
	if len(roots) <= 1:
		return set()
	connected, annotated = set(), set()
	for parent, child in conn.execute("SELECT parent, child FROM tree WHERE parent <> child"):
		connected.add(parent)
		connected.add(child)
	for node_id in conn.execute("SELECT DISTINCT node_id FROM snp_annotation"):
		annotated.add(node_id[0])
	orphans = set()
	for root in roots:
		if root not in connected and root not in annotated:
			orphans.add(root)
	return orphans


def _copy_table(source, dest, table, drop_rows=None):
	'''Copy every row of table from source to dest, keeping only the columns the two
		databases share so a column added or renamed between releases never breaks an import.'''
	source_columns = _columns(source, table)
	dest_columns = set(_columns(dest, table))
	shared = [column for column in source_columns if column in dest_columns]
	if not shared:
		logger.debug("table {table} has no columns in common, nothing copied".format(table=table))
		return 0
	rows = source.execute("SELECT {0} FROM {1}".format(",".join(shared), table)).fetchall()
	if drop_rows:
		rows = [row for row in rows if not drop_rows(row, source_columns)]
	if not rows:
		return 0
	indexes = [source_columns.index(column) for column in shared]
	dest.executemany("INSERT INTO {table} ({cols}) VALUES ({marks})".format(
		table=table, cols=",".join(shared), marks=",".join("?" * len(shared))),
		[tuple(row[i] for i in indexes) for row in rows])
	return len(rows)


def _touches(orphans, row, columns):
	'''a tree, annotation or node row that points at (or is) a dropped orphan node is not copied'''
	for column in ("id", "parent", "child", "node_id"):
		if column in columns and row[columns.index(column)] in orphans:
			return True
	return False


def import_database(source, dest, verbose=False):
	'''Write a copy of the database at source to dest in the current CanSNPer2 format and
		return a dict with the row counts. dest must not exist yet; source is never changed.'''
	if os.path.abspath(source) == os.path.abspath(dest):
		raise DatabaseImportError("--import needs a different destination database than the source")
	if os.path.exists(dest):
		raise DatabaseImportError("the destination database {0} already exists, choose a new name "
			"or remove it first".format(dest))
	try:
		check_database(source)
	except ConnectionError as e:
		raise DatabaseImportError(str(e))

	# build the destination with the schema the current create produces, so the import can
	# also upgrade a database made by an older release whose schema had drifted
	maker = CreateSNPDatabase(dest, create=True, verbose=verbose)
	maker.conn.close()

	source_conn = sqlite3.connect(source)
	dest_conn = sqlite3.connect(dest)
	counts = {}
	try:
		orphans = _find_orphan_roots(source_conn)
		if orphans:
			names = {row[0]: row[1] for row in source_conn.execute("SELECT id,name FROM nodes")}
			logger.info("Dropping {0} orphan root node(s): {1}".format(
				len(orphans), ", ".join("{0} (id {1})".format(names.get(o, "?"), o) for o in sorted(orphans))))
		drop_tree = lambda row, cols: _touches(orphans, row, cols)
		for table in IMPORT_TABLES:
			counts[table] = _copy_table(source_conn, dest_conn, table,
				drop_rows=drop_tree if table in ("tree", "snp_annotation", "nodes") else None)
		dest_conn.commit()

		# the new database must itself be exportable: one root, one parent per node, every SNP
		# tied to a real reference genome. if it is not, the source was not a clean scheme.
		try:
			tree_paths(dest_conn)
		except ExportError as e:
			raise DatabaseImportError("the imported database is not a single rooted tree: {0}; the "
				"source has more than one connected root or an unlinked node, edit its tree so a "
				"single root remains and import again".format(e))
	except Exception:
		# never leave a half written destination behind
		dest_conn.close()
		source_conn.close()
		try:
			os.remove(dest)
		except OSError:
			pass
		raise
	else:
		source_conn.close()
		dest_conn.close()

	counts["dropped_orphans"] = len(orphans)
	logger.info("Imported {db}: {nodes} nodes, {tree} links, {snp_annotation} SNPs, "
		"{snp_references} reference genomes, {dropped_orphans} orphan root node(s) dropped".format(db=dest, **counts))
	return counts
