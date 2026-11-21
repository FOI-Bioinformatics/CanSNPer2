"""Encoding of the source files: a UTF-8 byte order mark is accepted, other bytes are named.

A file saved as "UTF-8 with BOM" starts with U+FEFF, which was read as part of the first
header name or the root name and gave misleading errors. A file in another encoding such as
latin-1 ended the run with a bare UnicodeDecodeError that did not name the file."""
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
def test_non_utf8_source_file_is_named(tmp_path, key):
	files = write_sources(tmp_path)
	## "s1" is a strain in the SNP and references files, "A.1" a name in all of them
	## "A.1" is a name in the tree and SNP files, "s1" a strain in the references file
	set_bytes(files[key], b"s1" if key == "references" else b"A.1", b"\xe9\xe8")
	with pytest.raises(SourceFileError) as error:
		check_sources(tree=files["tree"], annotation=files["snp"], references=files["references"])
	assert files[key] in str(error.value) and "UTF-8" in str(error.value)
	with pytest.raises(SourceFileError) as error:
		create(tmp_path, files)
	assert files[key] in str(error.value) and "UTF-8" in str(error.value)
	assert not os.path.exists(str(tmp_path / "scheme.db"))


def test_non_utf8_source_file_is_named_by_import(tmp_path):
	files = write_sources(tmp_path)
	set_bytes(files["snp"], b"Paper2020", b"P\xe4per")
	with pytest.raises(DatabaseImportError) as error:
		import_from_sources(files["tree"], files["snp"], files["references"], str(tmp_path / "i.db"))
	assert files["snp"] in str(error.value) and "UTF-8" in str(error.value)
	assert not os.path.exists(str(tmp_path / "i.db"))


def test_non_utf8_modification_file_leaves_the_database_as_it_was(tmp_path):
	db = build_db(tmp_path, "m.db", tree_file("\t"), annotation_file())
	mod = tmp_path / "mod.txt"
	mod.write_bytes(b"parent\tchild\nB.1\tC\xe9.1\n")
	before = open(db, "rb").read()
	result = run_database_cli(tmp_path, modify_arguments(db, mod, "B.1"))
	assert result.returncode == 1
	assert str(mod) in result.stderr and "UTF-8" in result.stderr
	assert "Traceback" not in result.stderr
	assert open(db, "rb").read() == before


def test_non_utf8_modification_snp_file_is_named(tmp_path):
	db = build_db(tmp_path, "m.db", tree_file("\t"), annotation_file())
	snp = tmp_path / "snp.txt"
	snp.write_bytes(b"\t".join(s.encode() for s in helpers.SNP_HEADER) + b"\n"
		b"B.1\tfranc\xe9sella\tBirdsell2014\tREF1\tCP000001.1\t300\tT\tC\n")
	with pytest.raises(SourceFileError) as error:
		check_modification(db, annotation=str(snp))
	assert str(snp) in str(error.value) and "UTF-8" in str(error.value)
