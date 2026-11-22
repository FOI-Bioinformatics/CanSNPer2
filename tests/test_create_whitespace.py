"""Whitespace in the source files and SNPs that are not on a node of the tree.

The source check compared the names as written while the loaders stripped the tree
names, so a node such as "A " and its SNP were dropped with only a warning. Every field
is now stripped the same way in the checks and in the loaders."""
import sqlite3

import pytest

from test_create_input_checks import (SNP_HEADER, SNP_ROWS, REFERENCES, create, write_sources)
from CanSNPer2.modules.DatabaseConnection import SourceFileError
from CanSNPer2.modules.InitiateDatabase import CanSNPDatabase
from CanSNPer2.modules.SourceCheck import check_sources


def test_whitespace_in_tree_names_is_ignored_consistently(tmp_path):
	# "A " and "A" were two names to the source check and one to the loader, which then
	# dropped a node and its SNP with a warning
	tree = "R\nR;A ;B\n"
	files = write_sources(tmp_path, tree=tree, snp=SNP_HEADER + "R\t1\tP\tREF1\tCP000001.1\t10\tA\tG\t\tyes\n"
		"A\t1\tP\tREF1\tCP000001.1\t100\tT\tC\t\tyes\nB\t1\tP\tREF1\tCP000001.1\t200\tG\tA\t\tyes\n")
	assert check_sources(tree=files["tree"], annotation=files["snp"], references=files["references"]) == []
	db = create(tmp_path, files)
	conn = sqlite3.connect(db)
	try:
		assert conn.execute("SELECT COUNT(*) FROM snp_annotation").fetchone()[0] == 3
		assert conn.execute("SELECT COUNT(*) FROM nodes WHERE name LIKE '% '").fetchone()[0] == 0
	finally:
		conn.close()


def test_whitespace_makes_the_same_node_with_two_parents_visible(tmp_path):
	# "A " under B and "A" under R is one node with two parents once names are stripped
	files = write_sources(tmp_path, tree="R\nR;B\nR;X;A\nR;B;A ;C\n")
	problems = check_sources(tree=files["tree"])
	assert any("more than one parent" in p and "'A'" in p for p in problems), problems


def test_whitespace_in_snp_and_references_fields_is_ignored(tmp_path):
	snp = SNP_HEADER + SNP_ROWS.replace("A.1\t1", " A.1 \t1").replace("\tREF2\t", "\t REF2 \t")
	references = REFERENCES.replace("REF2\ts2", "REF2 \ts2")
	files = write_sources(tmp_path, snp=snp, references=references)
	assert check_sources(tree=files["tree"], annotation=files["snp"], references=files["references"]) == []
	db = create(tmp_path, files)
	conn = sqlite3.connect(db)
	try:
		assert conn.execute("SELECT COUNT(*) FROM snp_annotation").fetchone()[0] == 4
		assert conn.execute("SELECT COUNT(*) FROM snp_references WHERE genome = 'REF2'").fetchone()[0] == 1
	finally:
		conn.close()


def test_loader_raises_for_a_snp_whose_node_is_not_in_the_tree(tmp_path):
	# a SNP dropped with a warning left a database without it and exit status 0
	from CanSNPer2.modules.InitiateDatabase import CreateSNPDatabase
	files = write_sources(tmp_path)
	db = create(tmp_path, files)
	other = tmp_path / "other.txt"
	other.write_text(SNP_HEADER + "Z\t1\tP\tREF1\tCP000001.1\t10\tA\tG\t\tyes\n")
	loader = CanSNPDatabase(database=db, create=False, tree=None, annotation=None, references=None,
		source_type="CanSNPer", verbose=False)
	with pytest.raises(SourceFileError) as error:
		loader.load_annotation_file(str(other))
	assert "'Z'" in str(error.value)
