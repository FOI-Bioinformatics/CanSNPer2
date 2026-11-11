"""Regression tests for modules/NewickTree.py.

#27: a SNP node that no alignment carried a base for was skipped with a
continue in the old ete3 drawing loop, which also skipped set_style, so the node came
out of the pdf as an almost invisible default dot and the position was absent
from both per sample result files; the html viewer keeps that rule: the node is
drawn grey and reported.
#52: ROOT and the unannotated "-root-" label are not SNPs and must not get SNP dots.
issue #184: ete3 could not import on python >=3.13 (it imports cgi), and its ete4
replacement still needs PyQt6 or selenium+Chrome to write image files, so the tree
viewer is written by CanSNPer2 itself: a self-contained html file with an embedded
json payload and a plain javascript svg renderer."""
import json
import re
import sqlite3

from CanSNPer2.modules.NewickTree import NewickTree


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


def payload(newick_tree):
	'''the json the viewer embedded in the written html, the same node data the old
		ete3 NodeStyles carried'''
	html = open(newick_tree.tree_file).read()
	match = re.search(r"var DATA = (\{.*\});", html)
	assert match, "the html viewer must embed its json payload"
	return json.loads(match.group(1))


def test_snp_without_any_data_is_drawn_and_reported(tmp_path):
	# B.2 was never covered by an alignment, so it is missing from the snp list: it used to be
	# skipped by the drawing loop and left at the ete3 default size, invisible in the pdf (#27)
	newick_tree = tree(tmp_path)
	snplist = {"A.1": 1, "B.1": 1, "C.1": 2}
	final, msg, called, no_data = newick_tree.draw_tree(snplist, ["A.1", "B.1"], save_tree=True)
	assert [row[0] for row in no_data] == ["B.2"], "the position without data has to be reported"

	nodes = {n["name"]: n for n in payload(newick_tree)["nodes"]}
	# r keeps the ete3 scale: size was a diameter, 10 drew as r=5 and 15 as r=7.5
	assert nodes["B.2"]["fill"] == newick_tree.snp_colors["non_aligned"] and nodes["B.2"]["r"] == 5, \
		"B.2 without alignment data is drawn as a grey dot, never dropped"
	assert nodes["C.1"]["fill"] == newick_tree.snp_colors["ancestral"] and nodes["C.1"]["r"] == 5
	assert nodes["A.1"]["fill"] == newick_tree.snp_colors["derived"] and nodes["A.1"]["r"] == 7.5
	assert nodes["B.1"]["fill"] == newick_tree.snp_colors["derived"] and nodes["B.1"]["r"] == 7.5


def test_root_and_label_nodes_carry_no_snp_dot(tmp_path):
	# the ROOT wrapper and the "-root-" label that get_nodes invents for a root node without
	# snp_annotation are not SNPs and must not be drawn as grey SNP dots (#52, #32, #34)
	newick_tree = tree(tmp_path)
	newick_tree.draw_tree({"A.1": 1}, save_tree=True)
	nodes = {n["name"]: n for n in payload(newick_tree)["nodes"]}
	assert nodes["ROOT"]["dot"] is False
	assert nodes["-root-"]["dot"] is False
	'''the old ete default drew ROOT and the label node without a set_style: the viewer
		keeps them as plain junctions, and B.2 stays a SNP dot (#27)'''
	assert nodes["B.2"]["dot"] is True


def test_payload_mirrors_the_newick_structure(tmp_path):
	# the flattened tree is what the viewer draws: the ROOT wrapper above the database root,
	# preorder ids, and depths that reproduce the ete distances (unit branches) the SNP
	# calling sorted dlist by
	newick_tree = tree(tmp_path)
	newick_tree.draw_tree({}, save_tree=True)
	nodes = payload(newick_tree)["nodes"]
	entries = {n["name"]: n for n in nodes}
	by_id = {n["id"]: n for n in nodes}
	assert nodes[0]["name"] == "ROOT" and nodes[0]["parent"] is None and nodes[0]["depth"] == 0
	assert entries["-root-"]["parent"] == 0 and entries["-root-"]["depth"] == 1
	assert entries["A.1"]["name"] == "A.1" and by_id[entries["A.1"]["parent"]]["name"] == "-root-"
	assert entries["B.1"]["depth"] == 3 and entries["B.2"]["depth"] == 3 and entries["C.1"]["depth"] == 2
	parents = [entries[name]["parent"] for name in ("B.1", "B.2", "C.1", "A.1", "-root-")]
	assert all(by_id[n]["id"] < entries[name]["id"]
		for name, n in zip(("B.1", "B.2", "C.1", "A.1", "-root-"), parents)), "parents come before children (preorder)"


def test_no_data_rows_carry_the_database_annotation(tmp_path):
	# the same columns as the SNP_info rows written to <sample>_not_called.txt, the target
	# base being N/A because there is no base at that position at all (#27)
	newick_tree = tree(tmp_path)
	rows = newick_tree.get_no_data_info({"B.2", "C.1"})
	assert rows == [["B.2", "REF", "30", "C", "T", "N/A", ""], ["C.1", "REF", "40", "C", "T", "N/A", ""]], \
		"the empty seventh column keeps the shape of the SNP_FILE_COLUMNS tables (edge flag, issue #48)"
	assert newick_tree.get_no_data_info([]) == []


def test_every_snp_node_of_the_tree_ends_up_in_one_of_the_two_lists(tmp_path):
	# the two files together have to cover the tree: everything not evaluated is reported,
	# nothing is silently dropped (#27, and the coverage the reporter asked for in #35)
	newick_tree = tree(tmp_path)
	snplist = {"A.1": 1, "B.1": 2, "C.1": 0}
	_, _, called, no_data = newick_tree.draw_tree(snplist, ["A.1"], save_tree=False)
	evaluated = set(snplist)
	assert set(row[0] for row in no_data) == {"B.2"}, "B.2 is the only node without data"
	assert not (evaluated & set(row[0] for row in no_data)), "a position cannot be listed twice"


def test_nodes_without_snp_annotation_are_not_listed_as_snps(tmp_path):
	# the "-root-" label that get_nodes invents for a root node without snp_annotation is not a
	# SNP position: it belongs neither in the tree dots nor in the list of positions without
	# data, which the #27 fix would otherwise have started reporting (#52)
	newick_tree = tree(tmp_path)
	_, _, called, no_data = newick_tree.draw_tree({"A.1": 1, "B.1": 1}, ["A.1"], save_tree=False)
	names = set(row[0] for row in no_data)
	assert names == {"B.2", "C.1"}, names
	assert "-root-" not in names and "ROOT" not in names


def test_viewer_is_self_contained(tmp_path):
	# the html must carry everything it needs to render: embedded json and javascript, no
	# library or network references (issue #184: the ete3/Qt chain is gone for good)
	newick_tree = tree(tmp_path)
	newick_tree.draw_tree({"A.1": 1}, ["A.1"], save_tree=True)
	html = open(newick_tree.tree_file).read()
	assert "var DATA =" in html and "render()" in html
	assert "ete3" not in html and "ete4" not in html
	assert not re.search(r"src\s*=\s*[\"']https?://", html), "no external javascript may be referenced"


def test_node_names_cannot_close_the_script_block(tmp_path):
	# the payload is embedded in a <script> block, and json.dumps does not escape "</": a node
	# name from the database such as 'B.1</script><script>...' ended the block and the rest ran
	# as script when the tree was opened; the name must reach the viewer unchanged as data
	db = str(tmp_path / "t.db")
	make_db(db)
	name = 'B.1</script><script>document.title="injected"</script><!--'
	conn = sqlite3.connect(db)
	conn.execute("UPDATE snp_annotation SET snp_id = ? WHERE node_id = 3", (name,))
	conn.commit()
	conn.close()
	newick_tree = NewickTree(db, name="t", outdir=str(tmp_path))
	newick_tree.draw_tree({name: 1}, [name], save_tree=True)
	html = open(newick_tree.tree_file).read()
	assert html.lower().count("</script") == 1, "only the viewer's own script block may be closed"
	assert "<!--" not in html.split("var DATA =", 1)[1].split(";", 1)[0]
	assert name in [n["name"] for n in payload(newick_tree)["nodes"]]


def test_sample_name_is_escaped_in_the_title(tmp_path):
	# the title is the tree file name, which comes from the query file name; markup in it must
	# be shown as text, and a placeholder in it must not pull the payload into the title
	db = str(tmp_path / "t.db")
	make_db(db)
	newick_tree = NewickTree(db, name='<img src=x onerror=alert(1)>__PAYLOAD__', outdir=str(tmp_path))
	newick_tree.draw_tree({}, save_tree=True)
	html = open(newick_tree.tree_file).read()
	title = re.search(r"<title>(.*)</title>", html).group(1)
	assert "<img" not in title and "&lt;img" in title
	assert "__PAYLOAD__" in title and '"nodes"' not in title
	assert html.count("var DATA =") == 1


def call(tmp_path, snplist, called, strictness=0.7, min_required_hits=1):
	'''final call for the test tree: (snp, confirmed, called count, message)'''
	db = str(tmp_path / "t.db")
	if not (tmp_path / "t.db").exists():
		make_db(db)
	newick_tree = NewickTree(db, name="t", outdir=str(tmp_path), min_required_hits=min_required_hits, strictness=strictness)
	result, msg, _dlist, _no_data = newick_tree.draw_tree(snplist, called, save_tree=False)
	if not result:
		return None, False, 0, msg
	return result[1], result[2][0], result[2][1], msg


def make_deep_db(path):
	'''a single lineage R.1 -> R.2 -> R.3 -> R.4, the root R.1 annotated'''
	conn = sqlite3.connect(path)
	conn.executescript("""
		CREATE TABLE nodes (id integer PRIMARY KEY, name text NOT NULL);
		CREATE TABLE tree (parent integer, child integer, rank_i integer, unique (parent, child, rank_i));
		CREATE TABLE snp_annotation (node_id integer PRIMARY KEY, snp_id varchar(6), position integer,
			ancestral_base varchar(1), derived_base varchar(1), reference varchar(20),
			date datetime, genome_i integer);
		CREATE TABLE snp_references (id integer PRIMARY KEY, genome varchar(6), strain varchar(200));
		INSERT INTO snp_references (id, genome, strain) VALUES (5, 'REF', 'ref strain');
		INSERT INTO nodes (id, name) VALUES (1, 'R.1'), (2, 'R.2'), (3, 'R.3'), (4, 'R.4');
		INSERT INTO tree (parent, child, rank_i) VALUES (1, 1, 1), (1, 2, 1), (2, 3, 1), (3, 4, 1);
	""")
	for node_id in (1, 2, 3, 4):
		conn.execute("INSERT INTO snp_annotation (node_id, snp_id, position, ancestral_base, derived_base,"
			" reference, date, genome_i) VALUES (?, ?, ?, 'C', 'T', 'REF', '06/11/26', 5)", (node_id, "R.%d" % node_id, node_id * 10))
	conn.commit()
	conn.close()


def deep_call(tmp_path, snplist, strictness):
	db = str(tmp_path / "deep.db")
	if not (tmp_path / "deep.db").exists():
		make_deep_db(db)
	newick_tree = NewickTree(db, name="deep", outdir=str(tmp_path), min_required_hits=1, strictness=strictness)
	called = [name for name, state in snplist.items() if state == 1]
	result, msg, _dlist, _no_data = newick_tree.draw_tree(snplist, called, save_tree=False)
	return (result[1], result[2][0], result[2][1], msg) if result else (None, False, 0, msg)


def test_snp_without_data_on_the_path_is_tolerated_by_strictness(tmp_path):
	# R.4 derived, the root SNP R.1 has no alignment data: three of the four SNPs on the path
	# are called, so the call passes a strictness of 0.7 and fails one of 0.8 (one unaligned
	# SNP near the root made every sample NA, whatever --strictness said)
	snps = {"R.2": 1, "R.3": 1, "R.4": 1}
	snp, confirmed, count, msg = deep_call(tmp_path, snps, strictness=0.7)
	assert (snp, confirmed, count) == ("R.4", True, 3)
	assert "R.1" in msg, "the SNPs without data on the path are reported"
	snp, confirmed, count, msg = deep_call(tmp_path, snps, strictness=0.8)
	assert snp != "R.4" or confirmed is False


def test_unaligned_base_counts_as_missing_data(tmp_path):
	# state 0 is a gap or N in the query at the SNP position: no evidence either way
	snp, confirmed, count, msg = deep_call(tmp_path, {"R.1": 0, "R.2": 1, "R.3": 1, "R.4": 1}, strictness=0.7)
	assert (snp, confirmed) == ("R.4", True)


def test_ancestral_or_other_base_on_the_path_still_breaks_it(tmp_path):
	# an ancestral base (2) or a third base (3) is evidence against the path, not a gap
	for state in (2, 3):
		snp, confirmed, count, msg = call(tmp_path, {"A.1": state, "B.1": 1}, ["B.1"], strictness=0.1)
		assert confirmed is False, "state {0} on the path must not be tolerated".format(state)


def test_unannotated_root_label_does_not_break_the_path(tmp_path):
	# the "-root-" label node of a database whose root has no snp_annotation is not a SNP: a
	# complete called path below it was never confirmed, so such a database typed nothing
	snp, confirmed, count, msg = call(tmp_path, {"A.1": 1, "B.1": 1}, ["A.1", "B.1"], strictness=1.0, min_required_hits=2)
	assert (snp, confirmed, count) == ("B.1", True, 2)


def test_min_required_hits_counts_called_snps_only(tmp_path):
	snp, confirmed, count, msg = call(tmp_path, {"B.1": 1}, ["B.1"], strictness=0.5, min_required_hits=2)
	assert confirmed is False


def test_snp_directly_above_the_final_snp_must_be_called(tmp_path):
	# a gap is tolerated higher up the path, but not right above the final SNP: one stray
	# derived SNP on a neighbouring branch whose parent had no data would otherwise outrank the
	# sample's real lineage. Here B.1 is derived but its parent A.1 has no data, so B.1 is not
	# final-called even though 1/2 passes a strictness of 0.5
	snp, confirmed, count, msg = call(tmp_path, {"B.1": 1}, ["B.1"], strictness=0.5)
	assert confirmed is False
