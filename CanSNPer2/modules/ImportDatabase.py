'''
Build a CanSNPer2 database from the base CanSNPer2 source files.

``--import`` reads the source triple - the tree, SNP and reference files that
``--export`` writes and ``--create`` reads - and writes a database in the current
format to ``--database``:

    tree.txt + snp.txt + references.txt  --import-->  database

It is the exact counterpart of ``--export``: export a database, edit the text,
import it back. The same loaders ``--create`` uses do the loading, so the schema
can never drift from what the current release produces, and the result is then
checked with the export validator (one root, one parent per node, every SNP on
a real node tied to a real reference genome) - a triple that does not form a
clean scheme is refused and leaves no half-written database behind.
'''

import logging
import os
import sqlite3

from CanSNPer2.modules.InitiateDatabase import CanSNPDatabase
from CanSNPer2.modules.ExportDatabase import tree_paths, ExportError

logger = logging.getLogger(__name__)


class DatabaseImportError(Exception):
	'''the source files could not be imported into a current-format database'''


def import_from_sources(tree, annotation, references, database, tmpdir="/tmp", verbose=False, scheme_info=None, genomes=None):
	'''Build the database from the three source files and return a dict with the
		row counts. The destination must not exist yet; the source files are never
		changed. scheme_info: the known keys of the scheme_info table (scheme, status,
		source); format, built_with and built are always written. genomes: the optional genomes
		file (genome, node, role) of the scheme.'''
	if os.path.exists(database):
		raise DatabaseImportError("the destination database {0} already exists, choose a new name "
			"or remove it first".format(database))
	try:
		CanSNPDatabase(database=database, create=True, source_type="CanSNPer",
			tree=tree, annotation=annotation, references=references,
			tmpdir=tmpdir, verbose=verbose, scheme_info=scheme_info, genomes=genomes)
	except Exception as e:
		'''the loaders may refuse a source file (SourceFileError, InputError) or die on
			a broken one; either way a started database must not stay behind'''
		if os.path.exists(database):
			os.remove(database)
		raise DatabaseImportError(str(e))

	conn = sqlite3.connect(database)
	try:
		'''the imported database must itself be exportable: this is the same validator
			--export runs, so an export always round-trips through --import'''
		try:
			tree_paths(conn)
		except ExportError as e:
			raise DatabaseImportError("the imported files do not form a single rooted tree: {0}; "
				"edit the tree file so a single root remains and import again".format(e))
		counts = {table: conn.execute("SELECT COUNT(*) FROM {0}".format(table)).fetchone()[0]
			for table in ("nodes", "tree", "snp_annotation", "snp_references")}
		counts["typed"] = conn.execute("SELECT COUNT(*) FROM snp_annotation WHERE typed = 1").fetchone()[0]
	except Exception:
		conn.close()
		os.remove(database)
		raise
	conn.close()

	logger.info("Imported {db}: {nodes} nodes, {tree} links, {snp_annotation} canSNPs ({typed} typed), "
		"{snp_references} reference genomes".format(db=database, **counts))
	return counts
