"""Regression tests for modules/NewickTree.py, issues #27 and #19.

#27: a SNP node that no alignment carried a base for was skipped with a
continue in the drawing loop, which also skipped set_style, so the node came
out of the pdf as an almost invisible default dot and the position was absent
from both per sample result files.
#19: an ete3 without a Qt binding (an ete3 from the etetoolkit conda channel
without pyqt) has no drawing classes, and drawing died with "module 'ete3' has
no attribute 'NodeStyle'" in the middle of the node loop."""
import re
import sqlite3

import pytest

from CanSNPer2.modules import NewickTree as nwk
from CanSNPer2.modules.NewickTree import NewickTree

needs_qt = pytest.mark.skipif(not nwk.ETENEWICK_OK, reason="ete3 has no drawing classes, no Qt binding installed")


def make_db(path):
	'''root -> A.1 -> (B.1, B.2) and root -> C.1, annotated the way the shipped databases are:
		the root node has no snp_annotation row, every other node has its snp_id'''
	conn = sqlite3.connect(path)
	conn.executescript("""
		CREATE TABLE nodes (id integer PRIMARY KEY, name text NOT NULL);
		CREATE TABLE tree (parent integer, child integer, rank_i integer, unique (parent, child, rank_i));
		CREATE TABLE snp_annotation (node_id integer PRIMARY KEY, snp_id varchar(6), position integer,
			ancestral_base varchar(1), derived_base varchar(1), reference varchar(20),
			date datetime, genome_i integer);
		CREATE TABLE snp_references (id integer PRIMARY KEY, genome varchar(6), strain varchar(200));
		INSERT INTO snp_references (id, genome, strain) VALUES (5, 'REF', 'ref strain');
		INSERT INTO nodes (id, name) VALUES (1, 'root'), (2, 'A.1'), (3, 'B.1'), (4, 'B.2'), (5, 'C.1');
		INSERT INTO tree (parent, child, rank_i) VALUES (1, 1, 1), (1, 2, 1), (2, 3, 1), (2, 4, 1), (1, 5, 1);
	""")
	for node_id, snp_id, pos in ((2, "A.1", 10), (3, "B.1", 20), (4, "B.2", 30), (5, "C.1", 40)):
		conn.execute("INSERT INTO snp_annotation (node_id, snp_id, position, ancestral_base, derived_base,"
			" reference, date, genome_i) VALUES (?, ?, ?, 'C', 'T', 'REF', '06/11/26', 5)", (node_id, snp_id, pos))
	conn.commit()
	conn.close()


def tree(tmp_path):
	db = str(tmp_path / "t.db")
	make_db(db)
	return NewickTree(db, name="t", outdir=str(tmp_path))


@needs_qt
def test_snp_without_any_data_is_drawn_and_reported(tmp_path):
	# B.2 was never covered by an alignment, so it is missing from the snp list: it used to be
	# skipped by the drawing loop and left at the ete3 default size, invisible in the pdf (#27)
	newick_tree = tree(tmp_path)
	newick_tree.tree_file = str(tmp_path / "t_tree.svg")
	snplist = {"A.1": 1, "B.1": 1, "C.1": 2}
	final, msg, called, no_data = newick_tree.draw_ete3_tree(snplist, ["A.1", "B.1"], save_tree=True)
	assert [row[0] for row in no_data] == ["B.2"], "the position without data has to be reported"

	radii = re.findall(r'<circle[^>]*\br="([\d.]+)"', open(newick_tree.tree_file).read())
	# 1.5 is the ete3 default size the old code left the missing node at, 5 is the grey
	# "no data" dot, 7.5 the green derived ones; the six nodes are ROOT and the unannotated
	# "-root-" label (both without a style of their own, #52), A.1, B.1, B.2 (no data) and C.1
	assert radii.count("1.5") == 2, "ROOT and the unannotated root label keep no dot, got %s" % sorted(radii)
	assert radii.count("5") == 2, "B.2 (no data) and C.1 (ancestral) are grey: %s" % sorted(radii)
	assert radii.count("7.5") == 2, "A.1 and B.1 are derived: %s" % sorted(radii)


@needs_qt
def test_no_data_rows_carry_the_database_annotation(tmp_path):
	# the same columns as the SNP_info rows written to <sample>_not_called.txt, the target
	# base being N/A because there is no base at that position at all (#27)
	newick_tree = tree(tmp_path)
	rows = newick_tree.get_no_data_info({"B.2", "C.1"})
	assert rows == [["B.2", "REF", "30", "C", "T", "N/A", ""], ["C.1", "REF", "40", "C", "T", "N/A", ""]], \
		"the empty seventh column keeps the shape of the SNP_FILE_COLUMNS tables (edge flag, issue #48)"
	assert newick_tree.get_no_data_info([]) == []


@needs_qt
def test_every_snp_node_of_the_tree_ends_up_in_one_of_the_two_lists(tmp_path):
	# the two files together have to cover the tree: everything not evaluated is reported,
	# nothing is silently dropped (#27, and the coverage the reporter asked for in #35)
	newick_tree = tree(tmp_path)
	snplist = {"A.1": 1, "B.1": 2, "C.1": 0}
	_, _, called, no_data = newick_tree.draw_ete3_tree(snplist, ["A.1"], save_tree=False)
	evaluated = set(snplist)
	assert set(row[0] for row in no_data) == {"B.2"}, "B.2 is the only node without data"
	assert not (evaluated & set(row[0] for row in no_data)), "a position cannot be listed twice"


@needs_qt
def test_nodes_without_snp_annotation_are_not_listed_as_snps(tmp_path):
	# the "-root-" label that get_nodes invents for a root node without snp_annotation is not a
	# SNP position: it belongs neither in the tree dots nor in the list of positions without
	# data, which the #27 fix would otherwise have started reporting (#52)
	newick_tree = tree(tmp_path)
	_, _, called, no_data = newick_tree.draw_ete3_tree({"A.1": 1, "B.1": 1}, ["A.1"], save_tree=False)
	names = set(row[0] for row in no_data)
	assert names == {"B.2", "C.1"}, names
	assert "-root-" not in names and "ROOT" not in names


def test_ete3_without_qt_raises_a_readable_error(monkeypatch, tmp_path):
	# #19: ete3 imported but its drawing classes are missing (etetoolkit channel without pyqt):
	# the failure has to name the cause instead of dying with AttributeError on ete3.NodeStyle
	newick_tree = tree(tmp_path)
	monkeypatch.setattr(nwk, "ETENEWICK_OK", False)
	monkeypatch.setattr(nwk, "ETE_TREEVIEW_ERROR", "No module named 'PyQt5'")
	with pytest.raises(ImportError) as error:
		newick_tree.draw_ete3_tree({"A.1": 1}, ["A.1"], save_tree=False)
	text = str(error.value)
	assert "no drawing classes" in text and "PyQt5" in text


def test_drawing_classes_are_imported_from_where_they_live():
	# #19: the names are imported explicitly at module level, so no attribute of ete3 is
	# probed in the middle of the drawing loop any more
	if nwk.ETENEWICK_OK:
		for name in (nwk.NodeStyle, nwk.TreeStyle, nwk.AttrFace, nwk.faces):
			assert name is not None
		assert nwk.NodeStyle.__module__ == "ete3.treeview.main"
		assert nwk.AttrFace.__module__ == "ete3.treeview.faces"
