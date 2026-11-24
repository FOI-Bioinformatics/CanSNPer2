"""What the result pages show: summary.html (the results dashboard), <sample>_tree.html and
summary_tree.html. The page scripts are run under node with a small stand-in for the browser
document (page_harness.js), so the tests see what a browser would draw: the header with the
scheme and the run, the sample table, the status key, the legend, the SNP path list, the svg
of the tree, and the files the export buttons give."""
import html
import json
import os
import re
import shutil
import sqlite3
import subprocess
import xml.etree.ElementTree as ET

import pytest

import helpers
from CanSNPer2.modules import NewickTree as NT
from CanSNPer2.modules.NewickTree import NewickTree
from CanSNPer2.modules.SchemeInfo import write_scheme_info
from test_dashboard import REPORTS, make_db

NODE = shutil.which("node")
HARNESS = os.path.join(os.path.dirname(__file__), "page_harness.js")
needs_node = pytest.mark.skipif(NODE is None, reason="node is not installed")

RUN_INFO = {"version": "9.9.9", "date": "2026-11-22 10:00", "database_sha256": "ab" * 32, "aligner": "minimap2",
	"min_required_hits": 3, "strictness": 0.7, "mask": 0, "save_tree": True, "snp_files": True}


def run_page(path, after=""):
	result = subprocess.run([NODE, HARNESS, str(path), after], capture_output=True, text=True)
	assert result.returncode == 0, result.stderr
	out = json.loads(result.stdout)
	assert out["errors"] == [], out["errors"]
	return out


def scheme_db(tmp_path, status="verified"):
	db = str(tmp_path / "t.db")
	make_db(db)
	conn = sqlite3.connect(db)
	write_scheme_info(conn, {"scheme": "toy_scheme", "status": status, "source": "owner/data", "commit": "1234567890abcdef",
		"format": "2.3", "built_with": "2.3.0", "built": "2026-11-21"})
	conn.close()
	return db


def reports():
	extended = [dict(r) for r in REPORTS]
	extended[0].update({"support": "2/2", "path_no_data": [], "path_snps": [
		{"name": "A.1", "pos": "chrI:10", "ancestral": "C", "derived": "T", "target": "T"},
		{"name": "B.1", "pos": "chrI:20", "ancestral": "C", "derived": "T", "target": "T"}]})
	extended[1].update({"support": "2/2", "path_no_data": []})
	return extended


def dashboard(tmp_path, status="verified", annotations=None, run_info=RUN_INFO, samples=None):
	db = scheme_db(tmp_path, status)
	nt = NewickTree(db, name="t", outdir=str(tmp_path))
	return nt.write_dashboard(samples or reports(), str(tmp_path / "summary.html"), annotations, run_info=run_info)


def text_of(fragment):
	return re.sub(r"\s+", " ", html.unescape(re.sub(r"<[^>]+>", "", fragment)))


# 1. provenance

@needs_node
def test_the_dashboard_header_records_the_scheme_and_the_run(tmp_path):
	out = run_page(dashboard(tmp_path))
	prov = text_of(out["html"]["prov"])
	for part in ("toy_scheme (verified)", "owner/data@1234567", "REF", "t.db", "abababababab", "9.9.9", "minimap2",
			"min_required_hits 3", "strictness 0.7", "mask 0", "2026-11-22 10:00"):
		assert part in prov, part
	assert "CanSNPer2 results: toy_scheme" in out["html"]["title"]
	page = open(str(tmp_path / "summary.html")).read()
	assert "<title>CanSNPer2 results (toy_scheme)</title>" in page


@needs_node
@pytest.mark.parametrize("status", ["draft", "unknown"])
def test_a_scheme_that_is_not_verified_gets_a_banner(tmp_path, status):
	out = run_page(dashboard(tmp_path, status), "window.noticeShown = document.getElementById('notice').style.display;"
		"document.getElementById('prov').innerHTML += '<i>' + window.noticeShown + '</i>';")
	assert "<i>block</i>" in out["html"]["prov"]
	level, text = NT.SchemeInfo.scheme_notice(str(tmp_path / "t.db"))
	assert html.unescape(out["html"]["notice"]) == text, "the banner has the wording of the scheme notice"
	assert ("draft" if status == "draft" else "status of the scheme is unknown") in text


@needs_node
def test_a_verified_scheme_has_no_banner(tmp_path):
	out = run_page(dashboard(tmp_path), "document.getElementById('prov').innerHTML += '<i>' + document.getElementById('notice').style.display + '</i>';")
	assert "<i>undefined</i>" in out["html"]["prov"]


def test_the_dashboard_works_without_run_info(tmp_path):
	db = scheme_db(tmp_path)
	nt = NewickTree(db, name="t", outdir=str(tmp_path))
	page = open(nt.write_dashboard(reports(), str(tmp_path / "summary.html"))).read()
	data = json.loads(re.search(r"var DATA = (\{.*\});", page).group(1))
	assert data["info"]["scheme"] == "toy_scheme" and data["info"]["references"] == ["REF"]


@needs_node
def test_the_tree_page_has_a_title_and_a_header(tmp_path):
	db = scheme_db(tmp_path)
	nt = NewickTree(db, name="sample1", outdir=str(tmp_path), min_required_hits=1)
	info = [["A.1", "REF", "chrI", "10", "q", "1", "C", "T", "T", ""], ["B.1", "REF", "chrI", "20", "q", "2", "C", "T", "T", ""]]
	nt.draw_tree({"A.1": 1, "B.1": 1, "B.2": 2, "C.1": 2}, ["A.1", "B.1"], True, snp_info=info, run_info=RUN_INFO)
	page = open(nt.tree_file).read()
	assert "<title>CanSNPer2 tree: sample1 (toy_scheme)</title>" in page
	out = run_page(nt.tree_file)
	meta = text_of(out["html"]["meta"])
	assert "Sample sample1" in meta and "Final SNP B.1" in meta and "Support 2/2" in meta and "toy_scheme" in meta
	assert "minimap2" in text_of(out["html"]["prov"])
	items = re.findall(r'<li class="([^"]*)" data-snp="([^"]+)"[^>]*>(.*?)</li>', out["html"]["pathPanel"])
	assert [name for _cls, name, _ in items] == ["A.1", "B.1"], "the SNP path, root to final SNP"
	assert "final" in items[-1][0]
	assert "C&gt;T, found T" in items[0][2]


@needs_node
def test_the_summary_tree_page_names_the_final_snps(tmp_path):
	db = scheme_db(tmp_path)
	nt = NewickTree(db, name="summary", outdir=str(tmp_path))
	nt.draw_tree([], ["B.1", "C.1"], True, summary=True, run_info=RUN_INFO)
	assert "<title>CanSNPer2 tree: summary (toy_scheme)</title>" in open(nt.tree_file).read()
	out = run_page(nt.tree_file)
	assert "Final SNPs of the run: B.1, C.1" in text_of(out["html"]["meta"])


# 2. strength of each call

@needs_node
def test_the_table_shows_support_and_positions_without_data(tmp_path):
	out = run_page(dashboard(tmp_path))
	head = open(str(tmp_path / "summary.html")).read().split("<tbody")[0]
	assert re.findall(r"<th[^>]*>([^<]+)</th>", head) == ["Sample", "Final SNP", "Support", "No data", "Status"]
	cells = {}
	for row in out["rows"]:
		values = [text_of(c).strip() for c in re.findall(r"<td[^>]*>(.*?)</td>", row["html"])]
		cells[values[0]] = values
	assert cells["clean_sample"][1:4] == ["B.1", "2/2", "0"]
	assert cells["failed_sample"][1:4] == ["NA", "NA", "NA"]
	assert cells["failed_sample"][4].startswith("not called")


def test_the_status_key_explains_every_status(tmp_path):
	page = open(dashboard(tmp_path)).read()
	key = page.split('id="statusKey"', 1)[1].split("</details>", 1)[0]
	for status in ("clean", "review", "not called", "failed", "reused"):
		assert '">' + status + "</span>" in key, status
	assert "alignment edge" in key and "outside" in key and "no alignment data" in key


def test_the_run_puts_support_and_no_data_into_the_report(tmp_path):
	from test_typing_inputs import project, write, fake_alignment
	q = write(tmp_path / "x.fasta")
	obj, db = project(tmp_path, [q], min_required_hits=2)
	fake_alignment(obj, ["A.1", "B.1"], {"A.1": 1, "B.1": 1, "B.2": 2, "C.1": 2}, q)
	obj.run(database=db)
	report = obj.sample_reports["x"]
	assert report["final"] == "B.1"
	assert report["support"] == "2/2"
	assert report["path_no_data"] == []
	assert [p["name"] for p in report["path_snps"]] == ["A.1", "B.1"]
	assert report["path_snps"][0] == {"name": "A.1", "pos": "chr1:1", "ancestral": "C", "derived": "T", "target": "T"}


# 3. the tree figure

def svg_parts(svg):
	texts = []
	for attrs, content in re.findall(r"<text ([^>]*)>(.*?)</text>", svg):
		x = float(re.search(r'x="([-0-9.]+)"', attrs).group(1))
		y = float(re.search(r'y="([-0-9.]+)"', attrs).group(1))
		label = html.unescape(re.sub(r"<[^>]+>", "", content))
		w = len(label) * 6.0
		if 'text-anchor="middle"' in attrs:
			x -= w / 2
		elif 'text-anchor="end"' in attrs:
			x -= w
		texts.append((label, x, y - 9, x + w, y + 3))
	triangles = []
	for d in re.findall(r'<path d="(M [-0-9. ]+ L [-0-9. ]+ L [-0-9. ]+ Z)"', svg):
		v = [float(t) for t in re.findall(r"[-0-9.]+", d)]
		triangles.append(("triangle", v[0], min(v[3], v[5]), v[2], max(v[3], v[5])))
	lines = []
	for px, py, ny, nx in re.findall(r'<path d="M ([-0-9.]+) ([-0-9.]+) V ([-0-9.]+) H ([-0-9.]+)"', svg):
		px, py, ny, nx = float(px), float(py), float(ny), float(nx)
		lines.append(("line", px, min(py, ny), px, max(py, ny)))
		lines.append(("line", px, ny, nx, ny))
	circles = re.findall(r'<circle class="node" data-id="(\d+)" cx="([-0-9.]+)" cy="([-0-9.]+)" r="([-0-9.]+)" fill="([^"]+)" stroke="([^"]+)"', svg)
	return texts, triangles, circles, lines


def overlaps(a, b):
	return a[1] < b[3] and b[1] < a[3] and a[2] < b[4] and b[2] < a[4]


def big_db(path, leaves_per_lineage=14, lineages=8):
	'''a tree wider than the 60 leaves of the overview, with long SNP names: Root, lineages
		below it, SNPs below each lineage and two SNPs below the first SNP of each lineage'''
	rows = [("Root", "chrI", 1, "C", "T")]
	links = []
	pos = 2
	for k in range(lineages):
		lineage = "Lineage_with_a_long_name_{0}".format(k)
		rows.append((lineage, "chrI", pos, "C", "T")); pos += 1
		links.append((1, len(rows)))
		parent = len(rows)
		for j in range(leaves_per_lineage):
			rows.append(("SNP_{0}_{1}_named_at_length".format(k, j), "chrI", pos, "C", "T")); pos += 1
			links.append((parent, len(rows)))
			if j == 0:
				child = len(rows)
				for m in range(2):
					rows.append(("Deep_{0}_{1}".format(k, m), "chrI", pos, "C", "T")); pos += 1
					links.append((child, len(rows)))
	helpers.make_db(path, rows)
	conn = sqlite3.connect(path)
	conn.execute("INSERT INTO tree (parent, child, rank_i) VALUES (1, 1, 1)")
	conn.executemany("INSERT INTO tree (parent, child, rank_i) VALUES (?, ?, 1)", links)
	conn.commit()
	conn.close()
	return [r[0] for r in rows]


def big_tree_page(tmp_path):
	db = str(tmp_path / "big.db")
	names = big_db(db)
	nt = NewickTree(db, name="big", outdir=str(tmp_path), min_required_hits=1)
	states = {name: 2 for name in names}
	path = ["Root", "Lineage_with_a_long_name_3", "SNP_3_0_named_at_length", "Deep_3_1"]
	states.update({name: 1 for name in path})
	del states["SNP_5_4_named_at_length"]
	nt.draw_tree(states, path, True)
	return nt.tree_file, path


@needs_node
@pytest.mark.parametrize("action", ["", "document.getElementById('expandAll').onclick();"], ids=["overview", "expanded"])
def test_labels_and_triangles_do_not_overlap(tmp_path, action):
	page, _path = big_tree_page(tmp_path)
	out = run_page(page, action)
	texts, triangles, _circles, lines = svg_parts(out["html"]["tree"])
	if not action:
		assert triangles, "the wide tree opens collapsed"
	boxes = texts + triangles
	for i, a in enumerate(boxes):
		for b in boxes[i + 1:]:
			assert not overlaps(a, b), (a, b)
	## no branch line runs through a label
	for text in texts:
		inner = (text[0], text[1] + 1, text[2] + 1, text[3] - 1, text[4] - 1)
		for line in lines:
			assert not overlaps(inner, (line[0], line[1], line[2], line[3] + 0.01, line[4] + 0.01)), (text, line)


@needs_node
def test_the_states_have_distinct_marks(tmp_path):
	page, _path = big_tree_page(tmp_path)
	out = run_page(page, "document.getElementById('expandAll').onclick();")
	data = json.loads(re.search(r"var DATA = (\{.*\});", open(page).read()).group(1))
	colors = data["colors"]
	by_id = {n["id"]: n for n in data["nodes"]}
	marks = {}
	for node_id, _x, _y, r, fill, stroke in svg_parts(out["html"]["tree"])[2]:
		marks[by_id[int(node_id)]["name"]] = (fill, stroke, float(r))
	assert marks["Deep_3_1"][0] == colors["derived"]
	assert marks["Deep_3_0"][0] == colors["ancestral"]
	nodata = marks["SNP_5_4_named_at_length"]
	assert nodata[0] == "#ffffff" and nodata[1] == colors["outline"], "no alignment data is an open circle"
	assert len({colors["derived"], colors["ancestral"], colors["other_base"], "#ffffff"}) == 4


@needs_node
def test_untested_snps_of_the_summary_tree_do_not_look_ancestral(tmp_path):
	db = scheme_db(tmp_path)
	nt = NewickTree(db, name="summary", outdir=str(tmp_path))
	nt.draw_tree([], ["B.1"], True, summary=True)
	out = run_page(nt.tree_file)
	data = json.loads(re.search(r"var DATA = (\{.*\});", open(nt.tree_file).read()).group(1))
	by_id = {n["id"]: n["name"] for n in data["nodes"]}
	fills = {by_id[int(c[0])]: c[4] for c in svg_parts(out["html"]["tree"])[2]}
	assert fills["B.1"] == data["colors"]["derived"]
	assert fills["C.1"] == "#ffffff" and fills["C.1"] != data["colors"]["ancestral"]


@needs_node
def test_the_tree_is_fitted_with_the_snp_path_in_view(tmp_path):
	page, path = big_tree_page(tmp_path)
	out = run_page(page)
	x, y, w, h = [float(v) for v in out["viewBox"].split()]
	data = json.loads(re.search(r"var DATA = (\{.*\});", open(page).read()).group(1))
	by_name = {n["name"]: n["id"] for n in data["nodes"]}
	circles = {int(c[0]): (float(c[1]), float(c[2])) for c in svg_parts(out["html"]["tree"])[2]}
	for name in path:
		cx, cy = circles[by_name[name]]
		assert x <= cx <= x + w and y <= cy <= y + h, name
	assert "stroke-width=\"3\"" in out["html"]["tree"], "the SNP path is drawn as a thick line"


@needs_node
def test_the_legend_shows_the_states_on_both_pages(tmp_path):
	out = run_page(dashboard(tmp_path))
	legend = text_of(out["html"]["legend"])
	for word in ("Derived", "Ancestral", "Other base", "No alignment data", "SNP path"):
		assert word in legend
	page, _path = big_tree_page(tmp_path)
	legend = text_of(run_page(page)["html"]["legend"])
	for word in ("Derived", "Ancestral", "Other base", "No alignment data", "SNP path"):
		assert word in legend


# 4. the selected sample

@needs_node
def test_the_selected_sample_lists_its_snp_path_and_short_links(tmp_path):
	out = run_page(dashboard(tmp_path), "paint(0);")
	detail = out["html"]["detail"]
	assert [n for n in re.findall(r'data-snp="([^"]+)"', detail)] == ["A.1", "B.1"]
	assert re.search(r'<a href="samples/clean_sample_snps.txt" title="samples/clean_sample_snps.txt">called SNPs</a>', detail)
	assert re.search(r'<a href="samples/clean_sample_not_called.txt" title="samples/clean_sample_not_called.txt">not called</a>', detail)
	assert re.search(r'<a href="samples/clean_sample_tree.html" title="samples/clean_sample_tree.html">tree</a>', detail)
	assert "A.1;B.1" not in detail, "the path is a list, not a semicolon string"


@needs_node
def test_an_annotation_label_is_shown_next_to_the_final_snp(tmp_path):
	annotations = [{"node": "B.1", "kind": "node", "rank": "type", "label": "Type B"},
		{"node": "A.1", "kind": "clade", "rank": "group", "label": "Group A"}]
	out = run_page(dashboard(tmp_path, annotations=annotations), "paint(0);")
	line = out["html"]["detail"].split('<div class="line">', 1)[1].split("</div>", 1)[0]
	assert "Type B" in line and "Group A" in line


# 5. export

@needs_node
def test_the_table_downloads_as_tsv_with_the_run_as_comment_lines(tmp_path):
	out = run_page(dashboard(tmp_path), "document.getElementById('downloadTsv').onclick();")
	assert len(out["downloads"]) == 1
	download = out["downloads"][0]
	assert download["name"] == "cansnper2_results_toy_scheme.tsv"
	lines = download["text"].splitlines()
	comments = [l for l in lines if l.startswith("#")]
	assert any(l.startswith("# CanSNPer2 version: 9.9.9") for l in comments)
	assert any(l.startswith("# Scheme: toy_scheme (verified)") for l in comments)
	table = [l.split("\t") for l in lines if not l.startswith("#")]
	assert table[0] == ["sample", "final_snp", "support", "no_data", "status", "snp_path", "notes"]
	rows = {r[0]: r for r in table[1:]}
	assert rows["clean_sample"][1:6] == ["B.1", "2/2", "0", "clean", "A.1;B.1"]
	assert len(rows) == 3


@needs_node
@pytest.mark.parametrize("which", ["dashboard", "tree"])
def test_the_tree_saves_as_svg(tmp_path, which):
	if which == "dashboard":
		page = dashboard(tmp_path)
	else:
		page, _path = big_tree_page(tmp_path)
	out = run_page(page, "document.getElementById('saveSvg').onclick();")
	download = out["downloads"][0]
	assert download["name"].endswith("_tree.svg") and download["type"] == "image/svg+xml"
	root = ET.fromstring(download["text"].split("?>", 1)[1])
	assert root.tag == "{http://www.w3.org/2000/svg}svg"
	assert root.get("width") and root.get("viewBox")
	assert root.findall(".//{http://www.w3.org/2000/svg}circle")


def test_the_pages_have_a_print_stylesheet(tmp_path):
	page = open(dashboard(tmp_path)).read()
	assert "@media print" in page
	print_css = page.split("@media print", 1)[1]
	assert "display: none" in print_css and "#000" in print_css


# 6. wording

def visible_text(page):
	head = page.split("<script>")[0]
	return text_of(re.sub(r"<style>.*?</style>", "", head, flags=re.S))


def test_buttons_and_labels_are_in_sentence_case(tmp_path):
	page = open(dashboard(tmp_path)).read()
	for label in re.findall(r"<button[^>]*>([^<]+)</button>", page) + re.findall(r"<span>([^<]+)</span></label>", page):
		label = html.unescape(label)
		if label[0].isalpha():
			assert label[0].isupper() and label[1:] == label[1:].lower() or "TSV" in label or "SVG" in label or "SNP" in label, label
	for word in ("Overview", "Expand all", "Labels", "SNP paths only"):
		assert word in page


@pytest.mark.parametrize("which", ["dashboard", "tree"])
def test_the_pages_say_snp_path(tmp_path, which):
	if which == "dashboard":
		page = open(dashboard(tmp_path)).read()
	else:
		page = open(big_tree_page(tmp_path)[0]).read()
	text = page.split("var DATA =")[0] + page.split("var DATA =")[1].split("\n", 1)[1]
	assert "called path" not in text.lower()
	assert "chain" not in text.lower()


def test_the_templates_hold_no_backslash():
	# a backslash in the page script was unescaped by Python once and blanked both pages
	for template in (NT.BASE_CSS, NT.TREE_JS, NT.TREE_HTML_TEMPLATE, NT.DASHBOARD_HTML_TEMPLATE):
		assert "\\" not in template


def test_the_pages_are_self_contained(tmp_path):
	page = open(dashboard(tmp_path)).read()
	assert not re.search(r"(src|href)\s*=\s*[\"']https?://", page)
	assert "@import" not in page and "<link" not in page
