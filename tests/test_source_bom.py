"""A UTF-8 byte order mark at the start of a source file is ignored.

A file saved as "UTF-8 with BOM" starts with U+FEFF, which was read as part of the first
header name or the root name and gave misleading errors."""
import os
import sqlite3

import pytest

import helpers
from test_create_input_checks import REFERENCES, SNP_HEADER, SNP_ROWS, TREE, create, write_sources
from test_db_build_and_inputs import annotation_file, build_db, modify_arguments, run_database_cli, tree_file
from CanSNPer2.modules.DatabaseConnection import SourceFileError
from CanSNPer2.modules.ImportDatabase import DatabaseImportError, import_from_sources
from CanSNPer2.modules.SourceCheck import check_modification, check_sources

BOM = b"\xef\xbb\xbf"


def add_bom(path):
	with open(path, "rb") as handle:
		data = handle.read()
	with open(path, "wb") as handle:
		handle.write(BOM + data)


def set_bytes(path, old, new):
	with open(path, "rb") as handle:
		data = handle.read()
	assert old in data
	with open(path, "wb") as handle:
		handle.write(data.replace(old, new, 1))


@pytest.mark.parametrize("key", ["tree", "snp", "references"])
def test_bom_in_a_source_file_is_ignored(tmp_path, key):
	files = write_sources(tmp_path)
	add_bom(files[key])
	assert check_sources(tree=files["tree"], annotation=files["snp"], references=files["references"]) == []
	db = create(tmp_path, files)
	conn = sqlite3.connect(db)
	try:
		assert conn.execute("SELECT COUNT(*) FROM snp_annotation").fetchone()[0] == 4
		assert conn.execute("SELECT COUNT(*) FROM nodes WHERE name LIKE '%R'").fetchone()[0] == 1
		assert conn.execute("SELECT COUNT(*) FROM nodes WHERE name = 'R'").fetchone()[0] == 1
	finally:
		conn.close()


def test_bom_in_all_source_files_imports(tmp_path):
	files = write_sources(tmp_path)
	for path in files.values():
		add_bom(path)
	counts = import_from_sources(files["tree"], files["snp"], files["references"], str(tmp_path / "i.db"))
	assert counts["snp_annotation"] == 4


def test_bom_in_modification_files_is_ignored(tmp_path):
	db = build_db(tmp_path, "m.db", tree_file("\t"), annotation_file())
	mod = tmp_path / "mod.txt"
	mod.write_bytes(BOM + b"parent\tchild\nB.1\tC.1\n")
	result = run_database_cli(tmp_path, modify_arguments(db, mod, "B.1"))
	assert result.returncode == 0, result.stderr
	snp = tmp_path / "snp.txt"
	snp.write_bytes(BOM + b"\t".join(s.encode() for s in helpers.SNP_HEADER) + b"\n"
		b"C.1\tfrancisella\tBirdsell2014\tREF1\tCP000001.1\t300\tT\tC\n")
	assert check_modification(db, annotation=str(snp)) == []
