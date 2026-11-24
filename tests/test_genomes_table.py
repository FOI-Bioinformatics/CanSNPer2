"""The genomes of a scheme (<scheme>_genomes.txt: genome, node, role) in the database.
Every published scheme ships this file, but no database held it; read typing with local
references takes its candidate genomes from it."""
import os
import sqlite3
import subprocess
import sys

import pytest

import helpers
from CanSNPer2.modules.ExportDatabase import export_source
from CanSNPer2.modules.ImportDatabase import DatabaseImportError, import_from_sources
from CanSNPer2.modules.SchemeGenomes import read_genomes

TREE = "R\nR\tA.1\nR\tA.1\tA.2\nR\tB.1\n"
SNP_ROWS = [("R", 1, "P", "REF1", "chrI", 10, "A", "G", "", "yes"),
	("A.1", 1, "P", "REF1", "chrI", 20, "T", "C", "", "yes"),
	("A.2", 1, "P", "REF1", "chrI", 30, "G", "A", "", "yes"),
	("B.1", 1, "P", "REF1", "chrI", 40, "A", "C", "", "yes")]
GENOMES = ("genome\tnode\trole\tbasis\tsource\n"
	"GCA_000000011.1\tA.2\trepresentative\tpublished\tx\n"
	"GCA_000000012.1\tA.2\tmember\t\t\n"
	"GCA_000000013.1\tB.1\trepresentative\t\t\n"
	"GCA_000000014.1\t\toutgroup\t\t\n")


def build(tmp_path, genomes=GENOMES, name="db.sqlite"):
	files = helpers.write_sources(tmp_path, SNP_ROWS, tree=TREE, name="src")
	path = tmp_path / "genomes.txt"
	path.write_text(genomes)
	db = str(tmp_path / name)
	import_from_sources(files["tree"], files["snp"], files["references"], db, genomes=str(path))
	return db, files


def test_genomes_are_stored_with_node_and_role(tmp_path):
	db, _ = build(tmp_path)
	assert read_genomes(db) == [("GCA_000000011.1", "A.2", "representative"), ("GCA_000000012.1", "A.2", "member"),
		("GCA_000000013.1", "B.1", "representative"), ("GCA_000000014.1", "", "outgroup")]


def test_a_database_without_genomes_reads_as_none(tmp_path):
	files = helpers.write_sources(tmp_path, SNP_ROWS, tree=TREE, name="src")
	db = str(tmp_path / "plain.sqlite")
	import_from_sources(files["tree"], files["snp"], files["references"], db)
	assert read_genomes(db) == []
	legacy = helpers.make_db(tmp_path / "legacy.db", [("A.1", "chr", 5, "A", "G")])
	assert read_genomes(legacy) == [], "a database of an earlier version, genomes table without role or none"


@pytest.mark.parametrize("text, message", [
	("genome\tnode\trole\nG1\tNOPE\trepresentative\n", "NOPE"),
	("genome\tnode\trole\nG1\tA.1\trepresentative\nG1\tB.1\tmember\n", "G1.*more than once"),
	("genome\tnode\trole\nG1\tA.1\treference\n", "role reference"),
	("genome\tnode\nG1\tA.1\n", "header"),
	("genome\tnode\trole\nG1\tA.1\toutgroup\n", "outgroup.*no node"),
	("genome\tnode\trole\nG1\t\trepresentative\n", "G1.*no node"),
])
def test_bad_genomes_files_are_refused_and_leave_no_database(tmp_path, text, message):
	with pytest.raises(DatabaseImportError, match=message):
		build(tmp_path, genomes=text)
	assert not os.path.exists(tmp_path / "db.sqlite")


def test_export_writes_the_genomes_file_and_it_round_trips(tmp_path):
	db, _ = build(tmp_path)
	out = tmp_path / "export"
	out.mkdir()
	files = export_source(db, str(out), name="x")
	assert open(files["genomes"]).read() == ("genome\tnode\trole\n"
		"GCA_000000011.1\tA.2\trepresentative\nGCA_000000012.1\tA.2\tmember\n"
		"GCA_000000013.1\tB.1\trepresentative\nGCA_000000014.1\t\toutgroup\n")
	again = str(tmp_path / "again.sqlite")
	import_from_sources(files["tree"], files["snp"], files["references"], again, genomes=files["genomes"])
	assert read_genomes(again) == read_genomes(db)


def test_export_without_genomes_writes_no_genomes_file(tmp_path):
	files = helpers.write_sources(tmp_path, SNP_ROWS, tree=TREE, name="src")
	db = str(tmp_path / "plain.sqlite")
	import_from_sources(files["tree"], files["snp"], files["references"], db)
	out = tmp_path / "export"
	out.mkdir()
	assert "genomes" not in export_source(db, str(out), name="x")
	assert not (out / "x_genomes.txt").exists()


def test_command_line_create_with_genomes(tmp_path):
	files = helpers.write_sources(tmp_path, SNP_ROWS, tree=TREE, name="src")
	(tmp_path / "g.txt").write_text(GENOMES)
	db = str(tmp_path / "cli.sqlite")
	run = subprocess.run([sys.executable, "-m", "CanSNPer2.SNPDatabase", "--create", "--database", db, "--tree", files["tree"],
		"--annotation", files["snp"], "--references", files["references"], "--genomes", str(tmp_path / "g.txt"),
		"--logs", str(tmp_path / "logs")], capture_output=True, text=True)
	assert run.returncode == 0, run.stderr
	assert len(read_genomes(db)) == 4
