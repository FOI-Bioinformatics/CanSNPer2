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

import helpers
from CanSNPer2.modules.ExportDatabase import CANSNP_COLUMNS, export_source
from CanSNPer2.modules.ImportDatabase import DatabaseImportError, import_from_sources
from CanSNPer2.modules.InitiateDatabase import CanSNPDatabase

TREE = "R\nR\tA.1\nR\tA.1\tA.2\nR\tB.1\nR\tB.1\tB.2\n"
SNP_HEADER = helpers.CANSNP_HEADER
## the scores cover what the round trip must keep verbatim: none at all, an integer, a trailing
## zero and a negative value, none of which survives being stored as a number
SNP_ROWS = [("R", 1, "Paper2020", "REF1", "chrI", 10, "A", "G", "", "yes"),
	("A.1", 1, "Paper2020", "REF1", "chrII", 10, "T", "C", "520", "yes"),
	("A.1", 2, "Paper2020", "REF1", "chrII", 60, "A", "G", "0.980", "no"),
	("A.1", 3, "Paper2020", "REF1", "chrII", 70, "C", "T", "-0.5", "no"),
	("A.2", 1, "Paper2021", "REF1", "chrI", 200, "G", "A", "1", "yes"),
	("B.1", 1, "Paper2020", "REF2", "CP000002.1", 400, "A", "C", "", "yes")]
REFERENCES = ("genome\tstrain\tgenbank_id\trefseq_id\tassembly_name\n"
	"REF1\ts1\tGCA_000000001.1\tGCF_000000001.1\tASM1v1\n"
	"REF2\ts2\tGCA_000000002.1\tGCF_000000002.1\t\n")


def sources(tmp_path):
	return helpers.write_sources(tmp_path, SNP_ROWS, tree=TREE, references=REFERENCES, name="src",
		header=SNP_HEADER)


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
		"snps": sorted(conn.execute("SELECT snp_id, sequence, sequence_position FROM snp_annotation")),
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
	assert counts["snp_annotation"] == 6 and counts["typed"] == 4 and counts["snp_references"] == 2


def test_export_import_round_trip_keeps_sequences(tmp_path):
	'''two canSNPs at position 10, one on chrI and one on chrII, and the ranks, scores and
		typed flags of a node with three canSNPs, survive export, import and a second export
		unchanged'''
	db = clean_db(tmp_path)
	first = export_source(db, str(tmp_path / "first"), name="scheme")
	target = str(tmp_path / "imported.db")
	import_from_sources(first["tree"], first["snp"], first["references"], target)
	second = export_source(target, str(tmp_path / "second"), name="scheme")
	for key in ("tree", "snp", "references"):
		assert open(first[key], "rb").read() == open(second[key], "rb").read(), key
	lines = open(first["snp"]).read().splitlines()
	assert tuple(lines[0].split("\t")) == CANSNP_COLUMNS
	assert {tuple(line.split("\t")[4:6]) for line in lines[1:]} >= {("chrI", "10"), ("chrII", "10")}
	assert ("A.1", "chrII", 10) in content(target)["snps"]
	## the scores are written back as they were given, not as a number would print them
	scores = [line.split("\t")[8] for line in lines[1:]]
	assert scores == ["", "520", "0.980", "-0.5", "1", ""], scores
	assert [line.split("\t")[9] for line in lines[1:]] == ["yes", "yes", "no", "no", "yes", "yes"]
	assert [line.split("\t")[1] for line in lines[1:]] == ["1", "1", "2", "3", "1", "1"]


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
	assert "Imported" in result.stdout and "6 canSNPs (4 typed)" in result.stdout
	assert roots(str(target)) == 1


def test_the_command_line_import_demands_all_three_files(tmp_path):
	files = sources(tmp_path)
	result = subprocess.run([sys.executable, "-m", "CanSNPer2.SNPDatabase", "--import",
		"--tree", files["tree"], "--database", str(tmp_path / "cli.db"),
		"--logs", str(tmp_path / "logs")],
		capture_output=True, text=True, cwd=str(tmp_path))
	assert result.returncode == 1
	assert "needs all three source files" in result.stderr and "--annotation" in result.stderr
