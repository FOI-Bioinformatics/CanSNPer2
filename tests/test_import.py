"""CanSNPer2-database --import converts a legacy database into the current format.

Legacy databases carry an extra orphan root: the old tree reader added a node linked to
itself above the tree's own root (issues #32, #34, #149), so they have two self-linked
nodes and --export refuses them. --import drops the orphan and writes a database that the
current code can export and rebuild."""
import os
import sqlite3
import subprocess
import sys

import pytest

from CanSNPer2.modules.ExportDatabase import ExportError, export_source
from CanSNPer2.modules.ImportDatabase import DatabaseImportError, import_database
from CanSNPer2.modules.InitiateDatabase import CanSNPDatabase

TREE = "R\nR\tA.1\nR\tA.1\tA.2\nR\tB.1\nR\tB.1\tB.2\n"
SNPS = ("snp_id\treference\tgenome\tposition\tderived_base\tancestral_base\n"
	"R\tPaper2020\tREF1\t10\tA\tG\n"
	"A.1\tPaper2020\tREF1\t100\tT\tC\n"
	"A.2\tPaper2021\tREF1\t200\tG\tA\n"
	"B.1\tPaper2020\tREF2\t400\tA\tC\n")
REFERENCES = ("genome\tstrain\tgenbank_id\trefseq_id\tassembly_name\n"
	"REF1\ts1\tGCA_1.1\tGCF_1.1\tASM1v1\n"
	"REF2\ts2\tGCA_2.1\tGCF_2.1\t\n")


def clean_db(tmp_path):
	files = {}
	for key, content in (("tree", TREE), ("snp", SNPS), ("references", REFERENCES)):
		files[key] = str(tmp_path / ("src_" + key + ".txt"))
		open(files[key], "w").write(content)
	db = str(tmp_path / "scheme.db")
	CanSNPDatabase(database=db, create=True, tree=files["tree"], annotation=files["snp"],
		references=files["references"], source_type="CanSNPer", verbose=False)
	return db


def add_orphan_root(db):
	'''reproduce the legacy structure: a self-linked node connected to nothing (the old "root")'''
	conn = sqlite3.connect(db)
	conn.execute("INSERT INTO nodes (name) VALUES ('root')")
	orphan = conn.execute("SELECT last_insert_rowid()").fetchone()[0]
	conn.execute("INSERT INTO tree (parent, child, rank_i) VALUES (?, ?, 1)", (orphan, orphan))
	conn.commit()
	roots = conn.execute("SELECT count(*) FROM tree WHERE child=parent").fetchone()[0]
	conn.close()
	return roots


def content(db):
	conn = sqlite3.connect(db)
	result = {
		"names": sorted(n for (n,) in conn.execute("SELECT name FROM nodes")),
		"links": sorted((ch, p) for ch, p in conn.execute(
			"SELECT c.name, p.name FROM tree t JOIN nodes c ON c.id=t.child JOIN nodes p ON p.id=t.parent "
			"WHERE t.child<>t.parent")),
		"snps": sorted((s, pos) for s, pos in conn.execute("SELECT snp_id, position FROM snp_annotation")),
		"refs": sorted(g for (g,) in conn.execute("SELECT genome FROM snp_references")),
	}
	conn.close()
	return result


def test_a_clean_database_imports_to_the_same_database(tmp_path):
	db = clean_db(tmp_path)
	target = str(tmp_path / "imported.db")
	counts = import_database(db, target)
	assert counts["dropped_orphans"] == 0
	assert content(target) == content(db)
	assert sqlite3.connect(target).execute("SELECT count(*) FROM tree WHERE child=parent").fetchone()[0] == 1


def test_a_legacy_orphan_root_is_dropped(tmp_path):
	db = clean_db(tmp_path)
	original = content(db)
	assert add_orphan_root(db) == 2, "the legacy database should have two roots before import"
	with pytest.raises(ExportError, match="2 root nodes"):
		export_source(db, str(tmp_path / "export"))   # legacy databases cannot be exported

	target = str(tmp_path / "imported.db")
	counts = import_database(db, target)
	assert counts["dropped_orphans"] == 1
	assert sqlite3.connect(target).execute("SELECT count(*) FROM tree WHERE child=parent").fetchone()[0] == 1
	assert content(target) == original   # equals the database before the orphan was added
	# the imported database round trips: it exports and the orphan node is gone
	files = export_source(target, str(tmp_path / "export"))
	assert "root" not in content(target)["names"]
	assert len(open(files["tree"]).read().splitlines()) == len(content(target)["names"])


def test_the_source_is_never_modified(tmp_path):
	db = clean_db(tmp_path)
	add_orphan_root(db)
	before = content(db)
	import_database(db, str(tmp_path / "imported.db"))
	assert content(db) == before
	assert sqlite3.connect(db).execute("SELECT count(*) FROM tree WHERE child=parent").fetchone()[0] == 2


def test_the_destination_must_not_already_exist(tmp_path):
	db = clean_db(tmp_path)
	target = str(tmp_path / "imported.db")
	open(target, "w").close()
	with pytest.raises(DatabaseImportError, match="already exists"):
		import_database(db, target)


def test_the_source_must_be_a_database(tmp_path):
	with pytest.raises(DatabaseImportError, match="does not exist"):
		import_database(str(tmp_path / "nope.db"), str(tmp_path / "imported.db"))


def test_the_command_line_import(tmp_path):
	db = clean_db(tmp_path)
	add_orphan_root(db)
	target = tmp_path / "cli.db"
	result = subprocess.run([sys.executable, "-m", "CanSNPer2.SNPDatabase", "--import",
		"--import_db", db, "--database", str(target), "--logs", str(tmp_path / "logs")],
		capture_output=True, text=True, cwd=str(tmp_path))
	assert result.returncode == 0, result.stderr
	assert "Imported" in result.stdout and "orphan root node(s) dropped" in result.stdout
	assert sqlite3.connect(str(target)).execute("SELECT count(*) FROM tree WHERE child=parent").fetchone()[0] == 1
