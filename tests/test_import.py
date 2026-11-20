"""CanSNPer2-database --import builds a database from the source triple.

--import is the counterpart of --export: it reads the tree, SNP and reference
files (--tree, --annotation, --references) and writes a current-format database.
export -> edit -> import has to round-trip, and a broken triple is refused
without leaving a half-written database behind."""
import os
import sqlite3
import subprocess
import sys

import pytest

from CanSNPer2.modules.ExportDatabase import export_source
from CanSNPer2.modules.ImportDatabase import DatabaseImportError, import_from_sources
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


def sources(tmp_path):
	files = {}
	for key, content in (("tree", TREE), ("snp", SNPS), ("references", REFERENCES)):
		files[key] = str(tmp_path / ("src_" + key + ".txt"))
		open(files[key], "w").write(content)
	return files


def clean_db(tmp_path):
	files = sources(tmp_path)
	db = str(tmp_path / "scheme.db")
	CanSNPDatabase(database=db, create=True, tree=files["tree"], annotation=files["snp"],
		references=files["references"], source_type="CanSNPer", verbose=False)
	return db


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


def roots(db):
	return sqlite3.connect(db).execute("SELECT count(*) FROM tree WHERE child=parent").fetchone()[0]


def test_export_then_import_round_trips(tmp_path):
	db = clean_db(tmp_path)
	files = export_source(db, str(tmp_path / "export"))
	target = str(tmp_path / "imported.db")
	counts = import_from_sources(files["tree"], files["snp"], files["references"], target)
	assert content(target) == content(db)
	assert roots(target) == 1
	assert counts["snp_annotation"] == 4 and counts["snp_references"] == 2


def test_import_of_the_source_files_equals_create(tmp_path):
	files = sources(tmp_path)
	created = str(tmp_path / "created.db")
	imported = str(tmp_path / "imported.db")
	CanSNPDatabase(database=created, create=True, tree=files["tree"], annotation=files["snp"],
		references=files["references"], source_type="CanSNPer", verbose=False)
	import_from_sources(files["tree"], files["snp"], files["references"], imported)
	assert content(imported) == content(created)


def test_the_destination_must_not_already_exist(tmp_path):
	files = sources(tmp_path)
	target = str(tmp_path / "imported.db")
	open(target, "w").close()
	with pytest.raises(DatabaseImportError, match="already exists"):
		import_from_sources(files["tree"], files["snp"], files["references"], target)


def test_a_broken_triple_leaves_no_half_written_database(tmp_path):
	files = sources(tmp_path)
	'''the tree loads fine and the nodes are already in the new database when the reference
		file refuses its header: the started database must not stay behind'''
	open(files["references"], "w").write("wrong\theaders\there\nx\ty\tz\n")
	target = str(tmp_path / "imported.db")
	with pytest.raises(DatabaseImportError):
		import_from_sources(files["tree"], files["snp"], files["references"], target)
	assert not os.path.exists(target)


def test_a_missing_source_file_is_refused(tmp_path):
	files = sources(tmp_path)
	os.remove(files["snp"])
	with pytest.raises(DatabaseImportError):
		import_from_sources(files["tree"], files["snp"], files["references"], str(tmp_path / "x.db"))


def test_the_command_line_import(tmp_path):
	files = sources(tmp_path)
	target = tmp_path / "cli.db"
	result = subprocess.run([sys.executable, "-m", "CanSNPer2.SNPDatabase", "--import",
		"--tree", files["tree"], "--annotation", files["snp"], "--references", files["references"],
		"--database", str(target), "--logs", str(tmp_path / "logs")],
		capture_output=True, text=True, cwd=str(tmp_path))
	assert result.returncode == 0, result.stderr
	assert "Imported" in result.stdout and "4 SNPs" in result.stdout
	assert roots(str(target)) == 1


def test_the_command_line_import_demands_all_three_files(tmp_path):
	files = sources(tmp_path)
	result = subprocess.run([sys.executable, "-m", "CanSNPer2.SNPDatabase", "--import",
		"--tree", files["tree"], "--database", str(tmp_path / "cli.db"),
		"--logs", str(tmp_path / "logs")],
		capture_output=True, text=True, cwd=str(tmp_path))
	assert result.returncode == 1
	assert "needs all three source files" in result.stderr and "--annotation" in result.stderr
