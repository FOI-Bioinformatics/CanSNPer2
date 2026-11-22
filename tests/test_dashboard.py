"""Results overview dashboard (summary.html).

A run of 100+ genomes had no summarized starting point: stdout only listed final SNPs and
the per-sample html trees were one file each. The dashboard lists every sample, sorts the
ones that need a manual look (failed, not called, calls leaning on the strictness tolerance
over positions without data, paths over SNPs flagged at an alignment edge) to the top, and
embeds the database tree once so clicking a row repaints the tree with that sample's states.
The topology is the same for every sample, so only the SNP states are stored per sample."""
import json
import re
import sqlite3

import helpers
from CanSNPer2.modules.NewickTree import NewickTree


def make_db(path):
	"""root -> A.1 -> (B.1, B.2) and root -> C.1, same shape as test_newick_tree.py"""
	conn = sqlite3.connect(path)
	conn.executescript(helpers.SNP_ANNOTATION_TABLE + helpers.REFERENCES_TABLE + """
		CREATE TABLE nodes (id integer PRIMARY KEY, name text NOT NULL);
		CREATE TABLE tree (parent integer, child integer, rank_i integer, unique (parent, child, rank_i));
		INSERT INTO snp_references (id, genome, strain) VALUES (5, 'REF', 'ref strain');
		INSERT INTO nodes (id, name) VALUES (1, 'root'), (2, 'A.1'), (3, 'B.1'), (4, 'B.2'), (5, 'C.1');
		INSERT INTO tree (parent, child, rank_i) VALUES (1, 1, 1), (1, 2, 1), (2, 3, 1), (2, 4, 1), (1, 5, 1);
	""")
	for node_id, snp_id, pos in ((2, "A.1", 10), (3, "B.1", 20), (4, "B.2", 30), (5, "C.1", 40)):
		conn.execute("INSERT INTO snp_annotation (node_id, snp_id, sequence, sequence_position, ancestral_base, derived_base,"
			" reference, date, genome_i) VALUES (?, ?, 'chrI', ?, 'C', 'T', 'REF', '06/11/26', 5)", (node_id, snp_id, pos))
	conn.commit()
	conn.close()


REPORTS = [
	{"name": "clean_sample", "severity": "clean", "final": "B.1", "path": "A.1;B.1",
	 "reasons": [], "states": {"A.1": 1, "B.1": 1}, "called": ["A.1", "B.1"], "edge": [], "no_data": []},
	{"name": "edge_sample", "severity": "review", "final": "B.2", "path": "A.1;B.2",
	 "reasons": ["SNP flagged at an alignment edge on the SNP path: B.2"],
	 "states": {"A.1": 1, "B.2": 1}, "called": ["A.1", "B.2"], "edge": ["B.2"], "no_data": []},
	{"name": "failed_sample", "severity": "not called", "final": "NA", "path": "NA",
	 "reasons": ["Two or more confirmed SNP paths compete at the same depth"],
	 "states": {"A.1": 1}, "called": ["A.1"], "edge": [], "no_data": []},
]


def dashboard_payload(tmp_path):
	db = str(tmp_path / "t.db")
	make_db(db)
	nt = NewickTree(db, name="t", outdir=str(tmp_path))
	out = nt.write_dashboard(REPORTS, str(tmp_path / "summary.html"))
	page = open(out).read()
	match = re.search(r"var DATA = (\{.*\});", page)
	assert match, "the dashboard must embed its json payload"
	return page, json.loads(match.group(1))


def test_dashboard_lists_every_sample(tmp_path):
	page, data = dashboard_payload(tmp_path)
	assert [s["name"] for s in data["samples"]] == ["clean_sample", "edge_sample", "failed_sample"]
	for sample in REPORTS:
		assert sample["name"] in page


def test_dashboard_embeds_topology_once(tmp_path):
	"""the tree sits once in the file; per-sample data is only the SNP states, so a run of
		100+ samples does not multiply the tree"""
	_page, data = dashboard_payload(tmp_path)
	names = sorted(n["name"] for n in data["nodes"])
	assert names == ["-root-", "A.1", "B.1", "B.2", "C.1", "ROOT"]
	assert all("states" in s for s in data["samples"])
	for sample in data["samples"]:
		assert "nodes" not in sample


def test_dashboard_severity_survives(tmp_path):
	_page, data = dashboard_payload(tmp_path)
	sev = {s["name"]: s["severity"] for s in data["samples"]}
	assert sev == {"clean_sample": "clean", "edge_sample": "review", "failed_sample": "not called"}
	assert data["samples"][1]["edge"] == ["B.2"]


def test_dashboard_escaping(tmp_path):
	"""a database/sample name containing "<" must not close the embedded script block"""
	db = str(tmp_path / "t.db")
	make_db(db)
	nt = NewickTree(db, name="t", outdir=str(tmp_path))
	out = nt.write_dashboard([{"name": "x</script><b>", "severity": "failed", "final": "NA",
		"path": "NA", "reasons": [], "states": {}, "called": [], "edge": [], "no_data": []}],
		str(tmp_path / "summary.html"))
	page = open(out).read()
	assert "</script><b>" not in page
	match = re.search(r"var DATA = (\{.*\});", page)
	assert json.loads(match.group(1))["samples"][0]["name"] == "x</script><b>"


def test_dashboard_wires_only_existing_elements(tmp_path):
	"""the first dashboard died in the browser without rendering anything: its script wired
		zoomIn/zoomOut buttons that the header did not contain, and the TypeError stopped the
		script before the table and tree were drawn"""
	page, _data = dashboard_payload(tmp_path)
	head, script = page.split("<script>")
	for element_id in set(re.findall(r'getElementById\("([^"]+)"\)', script)):
		assert 'id="' + element_id + '"' in head, "the script wires #" + element_id + ", the page has no such element"
	assert "window.onerror" in script, "the page must report its own script failures"
	assert "stayed empty" in script, "the page must notice that it stayed empty"


def test_dashboard_table_keeps_every_column_in_view(tmp_path):
	# long sample and SNP names pushed the final SNP and call columns out of the 400 px side panel;
	# the columns have fixed widths and long names wrap inside their column (the start of long
	# names is often shared, so cutting them would hide what tells the samples apart)
	page, data = dashboard_payload(tmp_path)
	assert "table-layout: fixed" in page
	assert "<colgroup>" in page
	assert "overflow-wrap: anywhere" in page
	assert 'title="' in page.split("tr.innerHTML", 1)[1].split(";", 1)[0]
