'''
Module to read and write newick trees

'''
from CanSNPer2.modules.DatabaseConnection import CanSNPdbFunctions
from CanSNPer2.modules import SchemeInfo
from urllib.request import pathname2url
import html
import sqlite3
import json
import os
import re
import logging
'''the tree drawing used ete3 (issue #19: it silently lost its drawing classes when no Qt
	binding was installed, issues #27/#52: its node loop hid positions without alignment data)
	and ete3 3.1.3 cannot even import on python >=3.13 (it imports cgi, removed in 3.13, issue
	#184). ete4's replacement renderer needs PyQt6 or selenium+Chrome for static files, so the
	package no longer depends on either: draw_tree walks the NewickNode structure the database
	already built and writes a self-contained html file (embedded JSON plus a plain javascript
	SVG viewer). The file opens in any browser on any system, shows the same SNP colors, and
	can toggle the branches without called SNPs - the dynamic view the ete pdfs never allowed.
'''
logger = logging.getLogger(__name__)


_WARNED_NODES = set() ## nodes without a SNP that were already warned about in this process


def load_annotations(path):
	'''Read a tree annotation file: one annotation per line, tab separated
		<snp name>\t<node|clade>\t<rank>\t<label>[\t<color>]
		"node" draws a colored badge under the SNP dot, "clade" draws a dashed box around the
		whole subtree below it (a feature that one branch of the tree separates, like an
		antibiotic resistance marker at the split of two lineages).
		The rank is a free-text level carried by the annotation itself and is NEVER derived
		from the depth in the tree: the database tree has uneven resolution, so B.1 and B.16
		can both sit at subspecies level while standing at different depths.
		Empty lines and lines starting with # are skipped; rank, label and color of the
		header line ("snp", "snp_id", "node", "kind"...) are recognized and skipped.'''
	annotations = []
	first = True
	with open(path) as handle:
		for number, row in enumerate(handle, start=1):
			row = row.strip()
			if not row or row.startswith("#"):
				continue
			fields = row.split("\t")
			if first and fields[0].lower() in ("snp", "snp_id", "node", "name") and len(fields) > 1:
				first = False
				continue
			first = False
			if len(fields) < 4:
				raise ValueError("{path} line {n}: need <snp name>\\t<node|clade>\\t<rank>\\t<label>, found: {row}".format(
						path=path, n=number, row=row))
			kind = fields[1].strip().lower()
			if kind not in ("node", "clade"):
				raise ValueError("{path} line {n}: kind must be node or clade, found {kind}".format(
						path=path, n=number, kind=fields[1]))
			annotation = {"node": fields[0].strip(), "kind": kind, "rank": fields[2].strip(), "label": fields[3].strip()}
			if len(fields) > 4 and fields[4].strip():
				annotation["color"] = fields[4].strip()
			annotations.append(annotation)
	return annotations


def database_provenance(path):
	'''What a result page records about the database at path: the scheme name, status and
		description from its scheme_info table (SchemeInfo.describe), the notice that
		SchemeInfo.scheme_notice logs when the scheme is not verified (empty otherwise), the
		database file name and its reference genomes'''
	info = {"database": os.path.basename(str(path)), "scheme": "", "status": "", "scheme_line": "",
		"notice": "", "references": []}
	try:
		scheme = SchemeInfo.read_scheme_info(path)
		level, notice = SchemeInfo.scheme_notice(path)
	except sqlite3.Error as e:
		logger.debug("No scheme information read from {db}: {e}".format(db=path, e=e))
		return info
	info["scheme"] = scheme.get("scheme") or ""
	info["status"] = scheme.get("status") or ""
	info["scheme_line"] = SchemeInfo.describe(scheme) if scheme else SchemeInfo.MISSING_TEXT
	if level >= logging.WARNING:
		info["notice"] = notice
	try:
		conn = sqlite3.connect("file:{path}?mode=ro".format(path=pathname2url(os.path.abspath(str(path)))), uri=True)
		try:
			rows = conn.execute("SELECT genome, genbank_id, refseq_id FROM snp_references ORDER BY id").fetchall()
		finally:
			conn.close()
	except sqlite3.Error:
		rows = []
	for genome, genbank, refseq in rows:
		accession = genbank or refseq
		info["references"].append("{genome} ({acc})".format(genome=genome, acc=accession) if accession else str(genome))
	return info


## the drawing class of a SNP state: 1 derived, 2 ancestral, 3 another base, 0 no data
SNP_KINDS = {1: "derived", 2: "ancestral", 3: "other_base", 0: "no_data"}


def path_details(path, snp_info, snplist=None):
	'''The SNPs of a SNP path with what the result file gives for each: the reference sequence
		and position, the ancestral and derived base and the base found in the query (target).
		snp_info: rows in the columns of SNP_FILE_COLUMNS (CanSNPer2.py)'''
	rows = {}
	for row in snp_info or []:
		if len(row) > 8:
			rows.setdefault(row[0], row)
	details = []
	for name in path:
		row = rows.get(name)
		item = {"name": name}
		if row:
			item.update({"pos": "{seq}:{pos}".format(seq=row[2], pos=row[3]), "ancestral": row[6],
				"derived": row[7], "target": row[8]})
		details.append(item)
	return details


def provenance_lines(info):
	'''The run as "key: value" lines, in the order the result pages and snps_summary.txt
		show them; keys that info does not hold are left out'''
	info = info or {}
	lines = []
	if info.get("version"):
		lines.append("CanSNPer2 version: {v}".format(v=info["version"]))
	if info.get("date"):
		lines.append("Date: {d}".format(d=info["date"]))
	if info.get("database"):
		lines.append("Database: {db}{sha}".format(db=info["database"],
			sha=" sha256:" + info["database_sha256"] if info.get("database_sha256") else ""))
	if info.get("scheme_line"):
		lines.append("Scheme: {line}".format(line=info["scheme_line"]))
	if info.get("references"):
		lines.append("Reference genomes: {refs}".format(refs=", ".join(info["references"])))
	if info.get("aligner"):
		lines.append("Aligner: {a}".format(a=info["aligner"]))
	if "min_required_hits" in info:
		lines.append("Typing settings: min_required_hits={hits} strictness={strict} mask={mask}".format(
			hits=info["min_required_hits"], strict=info.get("strictness"), mask=info.get("mask")))
	if info.get("notice"):
		lines.append("Note: {n}".format(n=info["notice"]))
	return lines


__version__ = "2.0.1"
__author__ = "David Sundell"
__credits__ = ["David Sundell"]
__license__ = "GPLv3"
__maintainer__ = "FOI bioinformatics group"
__email__ = ["bioinformatics@foi.se", "david.sundell@foi.se"]
__date__ = "2019-06-07"
__status__ = "Production"
__partof__ = "CanSNPer2"

'''The result pages (<sample>_tree.html, summary_tree.html and summary.html) share their style
	sheet (BASE_CSS) and the tree drawing (TREE_JS); render_template puts both into a page
	template. The colored tree arrives as JSON in __PAYLOAD__, everything else is plain
	javascript drawing an svg, no libraries, no network, so the html files work the same on any
	machine and from file:// (issue #184 removed the ete3/Qt dependency that produced the old
	pdfs). The view is interactive: clicking a dot collapses or expands its subtree,
	collapse/expand/filters relayout the tree and refit it into the window, the viewBox pans by
	dragging and zooms with the wheel, and trees with more than 60 leaves open pre-collapsed to
	an overview. Branches that contain derived SNPs are never collapsed (auto or manual), and
	"Overview" returns to the starting auto-collapsed depth.
	Drawing: derived SNPs are filled green, ancestral grey, another base blue; a position
	without alignment data is an open circle and a SNP that was not tested (the summary tree,
	results reused from an earlier run) a small open circle, so it is not read as ancestral.
	The SNP path is drawn as a thick green line and listed next to the tree. The column width
	follows the longest SNP name and a collapsed branch keeps its triangle inside the one to
	three rows it is given, so labels and triangles do not overlap.
	The templates are raw strings and the javascript holds no backslash: in a normal string
	Python unescaped them (\\" became "), which once blanked both pages.'''
BASE_CSS = r''':root { --ink: #1f2937; --muted: #4b5563; --faint: #6b7280; --line: #e3e6e9; --panel: #fbfcfd;
        --accent: #2563eb; --amber-bg: #fff7e6; --amber-line: #f0c36d; --amber-ink: #7a4a00; }
* { box-sizing: border-box; }
html, body { height: 100%; margin: 0; background: #fff; }
body { font-family: system-ui, -apple-system, "Segoe UI", Roboto, Helvetica, Arial, sans-serif; color: var(--ink);
       font-size: 13px; display: flex; flex-direction: column; }
#head { padding: 10px 16px 8px 16px; border-bottom: 1px solid var(--line); background: #fff; }
#head h1 { margin: 0; font-size: 17px; font-weight: 600; }
#head .meta { margin-top: 3px; color: var(--ink); }
#head .meta b { font-weight: 600; }
#prov { display: flex; flex-wrap: wrap; gap: 2px 16px; margin-top: 5px; font-size: 12px; color: var(--muted); }
#prov .k { color: var(--faint); }
#notice { display: none; margin-top: 7px; padding: 6px 10px; border: 1px solid var(--amber-line); border-radius: 6px;
          background: var(--amber-bg); color: var(--amber-ink); font-size: 12px; }
#bar { display: flex; align-items: center; gap: 6px; flex-wrap: wrap; padding: 6px 16px; background: var(--panel);
       border-bottom: 1px solid var(--line); }
.group { display: flex; align-items: center; gap: 5px; flex-wrap: wrap; }
.sep { width: 1px; height: 22px; background: var(--line); margin: 0 4px; }
button { appearance: none; border: 1px solid #d4d9de; background: #fff; color: var(--ink); border-radius: 6px;
         padding: 4px 10px; font: inherit; font-size: 12px; cursor: pointer; }
button:hover { background: #f1f4f7; border-color: #c2c8cf; }
button:active { background: #e8ecf0; }
#zoomIn, #zoomOut { width: 30px; padding: 4px 0; }
label.check { display: inline-flex; align-items: center; gap: 5px; font-size: 12px; color: #374151; cursor: pointer;
              padding: 3px 7px; border-radius: 6px; }
label.check:hover { background: #f1f4f7; }
input[type="checkbox"] { width: 13px; height: 13px; accent-color: #2e8b57; }
#legend { display: flex; gap: 4px 14px; align-items: center; flex-wrap: wrap; font-size: 12px; color: var(--muted); margin-left: auto; }
#legend .key { display: inline-flex; align-items: center; gap: 5px; }
#legend svg { width: 16px; height: 14px; }
.pill { display: inline-block; border-radius: 999px; padding: 1px 8px; font-size: 11px; color: #fff; }
#main { flex: 1; display: flex; min-height: 0; }
#tree { flex: 1; min-width: 0; position: relative; }
#tree svg { width: 100%; height: 100%; display: block; cursor: grab; }
#tree circle.node { cursor: pointer; }
#tree text { pointer-events: none; user-select: none; }
.panelTitle { font-size: 12px; font-weight: 600; color: var(--muted); margin: 0 0 4px 0; }
ol.path { margin: 0; padding: 0 0 0 26px; font-size: 12px; }
ol.path li { padding: 1px 0; cursor: pointer; }
ol.path li:hover .snp { text-decoration: underline; }
ol.path .snp { font-weight: 600; }
ol.path .bases { color: var(--faint); font-family: ui-monospace, Menlo, Consolas, monospace; font-size: 11px; margin-left: 4px; }
ol.path li.nodata .snp { font-weight: 400; color: var(--faint); }
ol.path li.final .snp { color: #1f6b43; }
.muted { color: var(--faint); }
#hint { position: absolute; bottom: 8px; right: 10px; color: var(--faint); font-size: 11px;
        background: rgba(255, 255, 255, 0.9); border: 1px solid #e5e7eb; border-radius: 999px; padding: 3px 10px; }
@media print {
  html, body { height: auto; display: block; }
  body, #head, #head .meta, #prov, #prov .k, #legend, .muted, ol.path .bases, td, th { color: #000 !important; }
  #bar .group, #bar .sep, #hint, #search, .noprint { display: none !important; }
  #bar { border: none; background: none; padding: 4px 0; }
  #legend { margin-left: 0; }
  #head { border: none; padding: 0 0 6px 0; }
  #notice { color: #000; border-color: #000; background: none; }
  #main { display: block; }
  #tree { height: 17cm; border: 1px solid #999; page-break-inside: avoid; }
}
'''

TREE_JS = r'''var C = DATA.colors, INFO = DATA.info || {};
var ROW = 22, MARGIN = 24, OVERVIEW_LEAVES = 60, CHAR = 6.6, MAXSCALE = 1.3;
var nodes = DATA.nodes;
var ANN = DATA.annotations || [];
var annByNode = {};
ANN.forEach(function(a){ (annByNode[a.node] = annByNode[a.node] || []).push(a); });
var nodeByName = {};
nodes.forEach(function(n){ if (!(n.name in nodeByName)) { nodeByName[n.name] = n; } });
var RANKCOLOR = {};
var RANKPALETTE = ["#8e44ad", "#d35400", "#2471a3", "#c0392b", "#16a085", "#7d6608", "#6c3483"];
function annColor(a){ if (a.color) { return a.color; } var r = a.rank || a.label; if (!RANKCOLOR[r]) { RANKCOLOR[r] = RANKPALETTE[Object.keys(RANKCOLOR).length % RANKPALETTE.length]; } return RANKCOLOR[r]; }
var children = {};
nodes.forEach(function(n){ if (n.parent !== null) { (children[n.parent] = children[n.parent] || []).push(n); } });
function isLabel(n){ return n.name === "ROOT" || (/^-.*-$/).test(n.name); }
/*the column width follows the longest SNP name, so a label never reaches the next column*/
var longest = 4;
nodes.forEach(function(n){ if (!isLabel(n)) { longest = Math.max(longest, String(n.name).length); } });
var COL = Math.max(90, Math.min(260, Math.round(longest * CHAR + 34)));
var DOT = {derived: [C.derived, "none", 0, 6], ancestral: [C.ancestral, "none", 0, 4.5], other_base: [C.other_base, "none", 0, 6],
	no_data: ["#ffffff", C.outline, 1.5, 4.5], untested: ["#ffffff", "#b8bec6", 1.2, 3.5]};
var onlyCalled = false, labels = true, vb = null, svgEl = null, lastSvg = "", overviewState = null;
function esc(s){ return String(s).replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;").replace(/"/g, "&quot;"); }
/*n.cb: the node or a node below it carries a derived SNP; such branches are never collapsed*/
function markCalled(){
	nodes.forEach(function(n){ n.cb = !!n.called; });
	for (var i = nodes.length - 1; i > 0; i--) { if (nodes[i].cb && nodes[i].parent !== null) { nodes[nodes[i].parent].cb = true; } }
}
function hasCalled(n){ return !!n.cb; }
function keep(n){ return !onlyCalled || hasCalled(n); }
function visKids(n){ if (n.collapsed) { return []; } return (children[n.id] || []).filter(keep); }
function hiddenKids(n){ return n.collapsed && (children[n.id] || []).some(keep); }
function canCollapse(n){ return (children[n.id] || []).some(keep) && !hasCalled(n); }
function subtreeSize(n){ var c = 0; (children[n.id] || []).forEach(function(k){ if (keep(k)) { c += 1 + subtreeSize(k); } }); return c; }
/*a collapsed branch takes one to three rows, by the size of the hidden subtree, and its
	triangle stays inside those rows*/
function triRows(n){ var c = subtreeSize(n); return c >= 50 ? 3 : (c >= 10 ? 2 : 1); }
function badges(n){ return (annByNode[n.name] || []).filter(function(a){ return a.kind !== "clade"; }); }
function layout(){
	nodes.forEach(function(n){ n.y = undefined; n.rows = 1; });
	var next = 0;
	(function place(n){
		n.x = n.depth * COL + MARGIN;
		var kids = visKids(n);
		if (!kids.length) {
			n.rows = hiddenKids(n) ? triRows(n) : 1;
			n.y = MARGIN + (next + (n.rows - 1) / 2) * ROW;
			next += n.rows + badges(n).length;
		} else { kids.forEach(place); n.y = (kids[0].y + kids[kids.length - 1].y) / 2; }
	})(nodes[0]);
	return next;
}
function edgeAttrs(n){
	if (n.onPath) { return ' stroke="' + C.derived + '" stroke-width="3"'; }
	if (n.kind === "derived" || n.kind === "ancestral" || n.kind === "other_base" || n.kind === "label") { return ' stroke="#6b7280" stroke-width="1.4"'; }
	return ' stroke="#c4c9cf" stroke-width="1.4" stroke-dasharray="4,3"';
}
function tipTitle(n){
	var t = esc(n.name);
	var words = {derived: "derived", ancestral: "ancestral", other_base: "other base", no_data: "no alignment data", untested: "not tested"};
	if (words[n.kind]) { t += ": " + words[n.kind]; }
	if (n.edge) { t += ", near an alignment edge"; }
	if (n.onPath) { t += ", on the SNP path"; }
	if (hasCalled(n)) { if (visKids(n).length) { t += " (branch with derived SNPs, stays open)"; } }
	else if (visKids(n).length || hiddenKids(n)) { t += " (click to " + (n.collapsed ? "expand" : "collapse") + ")"; }
	return t;
}
function triWidth(n){ return 26 + 14 * (n.rows - 1); }
function render(){
	layout();
	var parts = [];
	var maxX = 0, maxY = 0;
	nodes.forEach(function(n){ if (n.y !== undefined) { maxX = Math.max(maxX, n.x); maxY = Math.max(maxY, n.y); } });
	var W = maxX + COL + 260, H = maxY + 2 * MARGIN + ROW;
	parts.push('<rect x="-2000" y="-2000" width="' + (W + 4000) + '" height="' + (H + 4000) + '" fill="#ffffff"/>');
	[false, true].forEach(function(pass){
		nodes.forEach(function(n){
			if (n.parent === null || n.y === undefined || !keep(n) || !!n.onPath !== pass) { return; }
			var p = nodes[n.parent];
			parts.push('<path d="M ' + p.x + ' ' + p.y + ' V ' + n.y + ' H ' + n.x + '" fill="none"' + edgeAttrs(n) + '/>');
		});
	});
	ANN.forEach(function(a){
		if (a.kind !== "clade") { return; }
		var n = nodeByName[a.node];
		if (!n || n.y === undefined || !keep(n)) { return; }
		var col = annColor(a);
		var cx0 = n.x, cx1 = n.x + COL, cy0 = n.y, cy1 = n.y;
		(function walk(m){ (children[m.id] || []).forEach(function(k){ if (k.y === undefined || !keep(k)) { return; } cx1 = Math.max(cx1, k.x + COL); cy0 = Math.min(cy0, k.y); cy1 = Math.max(cy1, k.y); walk(k); }); })(n);
		parts.push('<rect x="' + (cx0 - 14) + '" y="' + (cy0 - ROW * 0.7) + '" width="' + (cx1 - cx0 + 28) + '" height="' + (cy1 - cy0 + ROW * 1.4) + '" rx="10" fill="' + col + '" fill-opacity="0.07" stroke="' + col + '" stroke-width="1.3" stroke-dasharray="6,4" pointer-events="none"><title>' + esc((a.rank ? a.rank + ": " : "") + a.label) + '</title></rect>');
		parts.push('<text x="' + (cx1 + 18) + '" y="' + (cy0 - ROW * 0.7 + 12) + '" fill="' + col + '" font-weight="600">' + esc(a.label) + '</text>');
	});
	nodes.forEach(function(n){
		if (n.y === undefined || !keep(n)) { return; }
		var triW = 0;
		if (hiddenKids(n)) {
			triW = triWidth(n);
			var h = n.rows * ROW / 2 - 3;
			parts.push('<path d="M ' + n.x + ' ' + n.y + ' L ' + (n.x + triW) + ' ' + (n.y - h) + ' L ' + (n.x + triW) + ' ' + (n.y + h) + ' Z" fill="#eceff2" stroke="#9ca3af" stroke-width="1"/>');
		}
		var r = 2.5;
		if (n.dot === false || n.kind === "label") {
			parts.push('<circle class="node" data-id="' + n.id + '" cx="' + n.x + '" cy="' + n.y + '" r="2.5" fill="#9ca3af"><title>' + esc(n.name) + '</title></circle>');
		} else {
			var d = DOT[n.kind] || DOT.untested;
			r = d[3];
			if (n.edge) { parts.push('<circle cx="' + n.x + '" cy="' + n.y + '" r="' + (r + 3) + '" fill="none" stroke="#d97706" stroke-width="2"/>'); }
			parts.push('<circle class="node" data-id="' + n.id + '" cx="' + n.x + '" cy="' + n.y + '" r="' + r + '" fill="' + d[0] + '" stroke="' + d[1] + '" stroke-width="' + d[2] + '"><title>' + tipTitle(n) + '</title></circle>');
		}
		if (labels && n.dot !== false && n.kind !== "label") {
			var weight = n.onPath ? ' font-weight="600"' : "";
			var fill = n.final ? C.derived_text : "#222222";
			if (visKids(n).length) {
				/*an inner SNP is labelled above the branch that leads to it, where no other line runs*/
				parts.push('<text x="' + (n.x - r - 3) + '" y="' + (n.y - 5) + '" text-anchor="end" fill="' + fill + '"' + weight + '>' + esc(n.name) + '</text>');
			} else {
				var extra = hiddenKids(n) ? " (" + subtreeSize(n) + " hidden)" : "";
				parts.push('<text x="' + (n.x + Math.max(r, triW) + 5) + '" y="' + (n.y + 4) + '" fill="' + fill + '"' + weight + '>' + esc(n.name) + '<tspan fill="#6b7280" font-weight="400">' + esc(extra) + '</tspan></text>');
			}
		}
		badges(n).forEach(function(a, k){
			var acol = annColor(a);
			var w = a.label.length * CHAR + 14;
			var cy = n.y + ROW * (k + 1);
			parts.push('<rect x="' + (n.x - w / 2) + '" y="' + (cy - 8) + '" width="' + w + '" height="16" rx="8" fill="' + acol + '" pointer-events="none"><title>' + esc((a.rank ? a.rank + ": " : "") + a.label) + '</title></rect>');
			parts.push('<text x="' + n.x + '" y="' + (cy + 4) + '" text-anchor="middle" fill="#ffffff">' + esc(a.label) + '</text>');
		});
	});
	lastSvg = parts.join("");
	document.getElementById("tree").innerHTML = '<svg xmlns="http://www.w3.org/2000/svg" preserveAspectRatio="xMidYMid meet" font-family="Helvetica, Arial, sans-serif" font-size="11">' + lastSvg + "</svg>";
	svgEl = document.querySelector("#tree svg");
	apply();
}
function labelEnd(n){
	if (!labels || n.kind === "label" || n.dot === false) { return n.x + 10; }
	var text = String(n.name).length + (hiddenKids(n) ? 12 : 0);
	return n.x + (hiddenKids(n) ? triWidth(n) : 6) + 8 + text * CHAR;
}
function bounds(list){
	var x0 = Infinity, y0 = Infinity, x1 = -Infinity, y1 = -Infinity;
	list.forEach(function(n){
		var ext = n.rows * ROW / 2 + 4;
		x0 = Math.min(x0, n.x); y0 = Math.min(y0, n.y - ext); x1 = Math.max(x1, labelEnd(n)); y1 = Math.max(y1, n.y + ext + ROW * badges(n).length);
	});
	return {x0: x0 - MARGIN, y0: y0 - MARGIN / 2, x1: x1 + MARGIN, y1: y1 + MARGIN / 2};
}
function shown(){ return nodes.filter(function(n){ return n.y !== undefined && keep(n); }); }
function apply(){ if (svgEl && vb) { svgEl.setAttribute("viewBox", vb.x + " " + vb.y + " " + vb.w + " " + vb.h); } }
function fitBox(box, most){
	if (!svgEl) { return 1; }
	most = most || MAXSCALE;
	var el = svgEl.getBoundingClientRect(), W = el.width || 900, H = el.height || 600;
	var bw = box.x1 - box.x0, bh = box.y1 - box.y0;
	var scale = Math.min(W / bw, H / bh, most);
	var w = W / scale, h = H / scale;
	/*a small tree or the SNP path is drawn from the top left at a readable size, not blown up*/
	var dx = scale >= most ? 0 : (w - bw) / 2, dy = scale >= most ? 0 : (h - bh) / 2;
	vb = {x: box.x0 - dx, y: box.y0 - dy, w: w, h: h};
	apply();
	return scale;
}
function fit(){ return fitBox(bounds(shown())); }
/*on load the whole tree is fitted; when that makes the labels too small to read, the view
	is fitted to the SNP path instead (the fit button shows the whole tree)*/
function fitInitial(){
	var scale = fit();
	var onPath = shown().filter(function(n){ return n.onPath; });
	if (scale >= 0.75 || !onPath.length) { return; }
	var box = bounds(onPath);
	box.x0 = Math.min(box.x0, MARGIN / 2);
	box.x1 += 2 * COL;
	var mid = (box.y0 + box.y1) / 2, half = Math.max((box.y1 - box.y0) / 2 + 3 * ROW, 12 * ROW);
	box.y0 = mid - half; box.y1 = mid + half;
	fitBox(box, 1);
}
function centerOn(n){
	if (!n || n.y === undefined || !vb) { return; }
	vb.x = n.x - vb.w / 3; vb.y = n.y - vb.h / 2; apply();
}
function zoomBy(f){ if (!vb) { return; } vb.x += vb.w * (1 - f) / 2; vb.y += vb.h * (1 - f) / 2; vb.w *= f; vb.h *= f; apply(); }
function leavesIfCollapsedAt(d){
	var c = 0;
	(function walk(n){
		var kids = (children[n.id] || []).filter(keep);
		if (!kids.length || (n.depth >= d && !hasCalled(n))) { c++; return; }
		kids.forEach(walk);
	})(nodes[0]);
	return c;
}
function autoCollapse(){
	nodes.forEach(function(n){ n.collapsed = false; });
	var maxDepth = Math.max.apply(null, nodes.map(function(n){ return n.depth; }));
	if (leavesIfCollapsedAt(maxDepth) <= OVERVIEW_LEAVES) { return; }
	var bestD = 0;
	for (var t = 1; t <= maxDepth; t++) { if (leavesIfCollapsedAt(t) <= OVERVIEW_LEAVES) { bestD = t; } else { break; } }
	if (bestD < 2) { return; }
	nodes.forEach(function(n){ if (n.depth >= bestD) { n.collapsed = canCollapse(n); } });
}
function relayoutFit(){ render(); fit(); }
function dotKey(kind){
	var d = DOT[kind];
	return '<svg viewBox="0 0 16 14"><circle cx="8" cy="7" r="' + Math.min(d[3], 5.5) + '" fill="' + d[0] + '" stroke="' + d[1] + '" stroke-width="' + d[2] + '"/></svg>';
}
function legendHtml(keys){
	var html = keys.map(function(k){
		if (k[0] === "path") { return '<span class="key"><svg viewBox="0 0 16 14"><line x1="1" y1="7" x2="15" y2="7" stroke="' + C.derived + '" stroke-width="3"/></svg>' + k[1] + "</span>"; }
		if (k[0] === "edge") { return '<span class="key"><svg viewBox="0 0 16 14"><circle cx="8" cy="7" r="5" fill="none" stroke="#d97706" stroke-width="2"/></svg>' + k[1] + "</span>"; }
		return '<span class="key">' + dotKey(k[0]) + k[1] + "</span>";
	}).join("");
	var seen = {};
	ANN.forEach(function(a){
		var key = (a.rank || a.label) + "|" + a.kind;
		if (seen[key]) { return; }
		seen[key] = 1;
		var col = annColor(a);
		if (a.kind === "clade") { html += '<span class="key"><span class="pill" style="background:transparent;border:1.5px dashed ' + col + ';color:' + col + '">' + esc(a.rank || a.label) + "</span>clade</span>"; }
		else { html += '<span class="key"><span class="pill" style="background:' + col + '">' + esc(a.rank || a.label) + "</span></span>"; }
	});
	document.getElementById("legend").innerHTML = html;
}
function shortSha(sha){ return sha ? sha.slice(0, 12) : ""; }
function provHtml(){
	var items = [];
	function add(k, v, title){ if (v === undefined || v === null || v === "") { return; } items.push('<span' + (title ? ' title="' + esc(title) + '"' : "") + '><span class="k">' + esc(k) + "</span> " + esc(v) + "</span>"); }
	add("Scheme", INFO.scheme_line);
	add("Reference genomes", (INFO.references || []).join(", "));
	add("Database", INFO.database ? INFO.database + (INFO.database_sha256 ? " (sha256 " + shortSha(INFO.database_sha256) + ")" : "") : "", INFO.database_sha256 ? "sha256:" + INFO.database_sha256 : "");
	add("CanSNPer2", INFO.version);
	add("Aligner", INFO.aligner);
	if (INFO.min_required_hits !== undefined) { add("Settings", "min_required_hits " + INFO.min_required_hits + ", strictness " + INFO.strictness + ", mask " + INFO.mask); }
	add("Run", INFO.date);
	document.getElementById("prov").innerHTML = items.join("");
	var notice = document.getElementById("notice");
	if (INFO.notice) { notice.textContent = INFO.notice; notice.style.display = "block"; }
}
function pathHtml(path, final){
	if (!path || !path.length) { return '<div class="muted">No SNP path.</div>'; }
	return '<ol class="path">' + path.map(function(p){
		var bases = p.ancestral ? p.ancestral + "&gt;" + p.derived + (p.target ? ", found " + esc(p.target) : "") : "";
		var cls = (p.name === final ? "final" : "") + (p.target === "N/A" ? " nodata" : "");
		return '<li class="' + cls + '" data-snp="' + esc(p.name) + '" title="' + esc(p.pos ? p.name + " at " + p.pos : p.name) + '"><span class="snp">' + esc(p.name) + '</span><span class="bases">' + bases + "</span></li>";
	}).join("") + "</ol>";
}
function wirePath(el){
	el.onclick = function(e){
		var t = e.target;
		while (t && t !== el && !(t.getAttribute && t.getAttribute("data-snp"))) { t = t.parentNode; }
		if (t && t !== el) { centerOn(nodeByName[t.getAttribute("data-snp")]); }
	};
}
function download(name, text, type){
	var blob = new Blob([text], {type: type});
	var a = document.createElement("a");
	a.href = URL.createObjectURL(blob);
	a.download = name;
	document.body.appendChild(a);
	a.click();
	document.body.removeChild(a);
	setTimeout(function(){ URL.revokeObjectURL(a.href); }, 1000);
}
function svgText(){
	var box = bounds(shown());
	var w = Math.ceil(box.x1 - box.x0), h = Math.ceil(box.y1 - box.y0);
	return '<?xml version="1.0" encoding="UTF-8"?>' + String.fromCharCode(10)
		+ '<svg xmlns="http://www.w3.org/2000/svg" width="' + w + '" height="' + h + '" viewBox="' + box.x0 + " " + box.y0 + " " + w + " " + h + '" font-family="Helvetica, Arial, sans-serif" font-size="11">'
		+ lastSvg + "</svg>" + String.fromCharCode(10);
}
document.getElementById("zoomIn").onclick = function(){ zoomBy(1 / 1.4); };
document.getElementById("zoomOut").onclick = function(){ zoomBy(1.4); };
document.getElementById("fit").onclick = fit;
document.getElementById("expandAll").onclick = function(){ nodes.forEach(function(n){ n.collapsed = false; }); relayoutFit(); };
document.getElementById("collapseAll").onclick = function(){ nodes.forEach(function(n){ if (n.depth >= 1) { n.collapsed = canCollapse(n); } }); relayoutFit(); };
document.getElementById("onlyCalled").onchange = function(){ onlyCalled = this.checked; relayoutFit(); };
document.getElementById("showLabels").onchange = function(){ labels = this.checked; render(); };
document.getElementById("overview").onclick = function(){
	if (!overviewState) { return; }
	nodes.forEach(function(n){ n.collapsed = overviewState[n.id]; });
	relayoutFit();
};
document.getElementById("saveSvg").onclick = function(){ download(svgName(), svgText(), "image/svg+xml"); };
var treeDiv = document.getElementById("tree");
var suppressClick = false;
treeDiv.addEventListener("click", function(e){
	if (suppressClick) { return; }
	var t = e.target;
	if (t.tagName !== "circle" || !t.hasAttribute("data-id")) { return; }
	var n = nodes[+t.getAttribute("data-id")];
	if (n.collapsed || canCollapse(n)) { n.collapsed = !n.collapsed; relayoutFit(); }
});
treeDiv.addEventListener("dblclick", function(e){ if (e.target.tagName !== "circle") { fit(); } });
treeDiv.addEventListener("wheel", function(e){
	e.preventDefault();
	if (!svgEl || !vb) { return; }
	var r = svgEl.getBoundingClientRect(), fx = (e.clientX - r.left) / r.width, fy = (e.clientY - r.top) / r.height;
	var f = e.deltaY > 0 ? 1.15 : 1 / 1.15;
	var px = vb.x + vb.w * fx, py = vb.y + vb.h * fy;
	vb.w *= f; vb.h *= f;
	vb.x = px - vb.w * fx; vb.y = py - vb.h * fy;
	apply();
}, {passive: false});
var drag = null, moved = 0;
treeDiv.addEventListener("mousedown", function(e){ drag = {x: e.clientX, y: e.clientY}; moved = 0; suppressClick = false; });
window.addEventListener("mousemove", function(e){
	if (!drag || !svgEl || !vb) { return; }
	var r = svgEl.getBoundingClientRect();
	var dx = e.clientX - drag.x, dy = e.clientY - drag.y;
	moved += Math.abs(dx) + Math.abs(dy);
	vb.x -= dx * vb.w / r.width; vb.y -= dy * vb.h / r.height;
	drag = {x: e.clientX, y: e.clientY};
	apply();
});
window.addEventListener("mouseup", function(){
	if (drag && moved > 5) { suppressClick = true; }
	drag = null;
});
'''


def render_template(template):
	'''a page template with the shared style sheet and tree drawing put in'''
	return template.replace("__BASECSS__", BASE_CSS).replace("__TREEJS__", TREE_JS)


'''The per-sample tree viewer written by NewickTree.write_tree_html: a header with the sample,
	the scheme, the final SNP and the run, the SNP path as an ordered list (root to final SNP,
	with the bases at each SNP) and the tree. The summary tree (--summary) shows the final SNPs
	of the run instead.'''
TREE_HTML_TEMPLATE = r'''<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>__TITLE__</title>
<style>
__BASECSS__
#pathPanel { width: 250px; flex: none; overflow-y: auto; padding: 10px 12px; border-right: 1px solid var(--line); background: #fdfefe; }
@media print { #pathPanel { width: auto; border: none; padding: 0 0 8px 0; } }
</style>
</head>
<body>
<div id="head">
<h1 id="title">CanSNPer2 tree</h1>
<div class="meta" id="meta"></div>
<div id="prov"></div>
<div id="notice"></div>
</div>
<div id="bar">
<div class="group">
<button id="overview" title="Back to the starting collapsed overview">Overview</button>
<button id="expandAll">Expand all</button>
<button id="collapseAll">Collapse</button>
<button id="fit" title="Fit the whole tree into the window">Fit</button>
<button id="zoomIn" title="Zoom in">+</button>
<button id="zoomOut" title="Zoom out">&minus;</button>
</div>
<span class="sep"></span>
<div class="group">
<label class="check" title="Show only the branches that lead to derived SNPs"><input type="checkbox" id="onlyCalled"><span>SNP paths only</span></label>
<label class="check"><input type="checkbox" id="showLabels" checked><span>Labels</span></label>
<button id="saveSvg" title="Save the tree as it is drawn now as an SVG file">Save SVG</button>
</div>
<div id="legend"></div>
</div>
<div id="main">
<div id="pathPanel"></div>
<div id="tree"></div>
</div>
<div id="hint" class="noprint">Click a dot to collapse or expand &middot; drag to pan &middot; wheel to zoom &middot; double click to fit &middot; branches with derived SNPs stay open</div>
<script>
var DATA = __PAYLOAD__;
window.onerror = function(message, source, line){
	var banner = document.createElement("div");
	banner.style.cssText = "background:#c0392b;color:#fff;padding:10px 14px;font-size:13px;";
	banner.textContent = "The CanSNPer2 tree page hit a script error: " + message + " (line " + line + "). Please report this message together with the CanSNPer2 version.";
	document.body.insertBefore(banner, document.body.firstChild);
};
__TREEJS__
var PAGE = DATA.page || {};
var pathNames = {};
if (PAGE.final && PAGE.final !== "NA") { (PAGE.path || []).forEach(function(p){ pathNames[p.name] = 1; }); }
nodes.forEach(function(n){
	if (isLabel(n)) { n.kind = "label"; n.dot = false; }
	n.kind = n.kind || "untested";
	n.onPath = !!pathNames[n.name];
	n.final = n.name === PAGE.final;
});
markCalled();
function svgName(){ return (PAGE.sample || "cansnper2") + "_tree.svg"; }
document.getElementById("title").textContent = "CanSNPer2 tree: " + (PAGE.sample || "") + (INFO.scheme ? " (" + INFO.scheme + ")" : "");
var meta = [];
if (INFO.scheme) { meta.push("Scheme <b>" + esc(INFO.scheme) + "</b>" + (INFO.status ? " (" + esc(INFO.status) + ")" : "")); }
if (PAGE.summary) {
	meta.push("Final SNPs of the run: <b>" + ((PAGE.finals || []).length ? esc(PAGE.finals.join(", ")) : "none") + "</b>");
} else {
	meta.push("Sample <b>" + esc(PAGE.sample || "") + "</b>");
	meta.push("Final SNP <b>" + esc(PAGE.final || "NA") + "</b>");
	if (PAGE.support !== undefined) { meta.push('<span title="Derived SNPs found on the SNP path / SNPs on the SNP path">Support ' + esc(PAGE.support) + "</span>"); }
	if (PAGE.no_data !== undefined && PAGE.no_data !== "NA") { meta.push("No data " + esc(PAGE.no_data)); }
}
document.getElementById("meta").innerHTML = meta.join(" &middot; ");
provHtml();
var panel = document.getElementById("pathPanel");
if (PAGE.summary) {
	panel.innerHTML = '<div class="panelTitle">Final SNPs</div>' + pathHtml((PAGE.finals || []).map(function(f){ return {name: f}; }), "");
	legendHtml([["derived", "Final SNP of a sample"], ["untested", "Other SNP"]]);
} else {
	var head = PAGE.final && PAGE.final !== "NA" ? "SNP path, root to final SNP" : "Derived SNPs (no SNP path confirmed)";
	panel.innerHTML = '<div class="panelTitle">' + head + "</div>" + pathHtml(PAGE.path || [], PAGE.final);
	legendHtml([["derived", "Derived"], ["ancestral", "Ancestral"], ["other_base", "Other base"], ["no_data", "No alignment data"], ["path", "SNP path"]]);
}
wirePath(panel);
autoCollapse();
overviewState = nodes.map(function(n){ return n.collapsed; });
render();
fitInitial();
if (!svgEl) { window.onerror("the tree did not draw at all", "", 0); }
</script>
</body>
</html>
'''


'''The results overview written by NewickTree.write_dashboard: one row per sample of the run
	in a side table, and the database tree embedded once next to it. The topology is the same
	for every sample, so the file stores only the per-sample SNP states and clicking a row
	repaints the tree with them. Samples that need a manual look (failed, not called, calls
	leaning on the strictness tolerance over positions without data, SNP paths over SNPs flagged
	at an alignment edge, derived SNPs outside the SNP path) sort to the top and the first one
	is selected on open. The header records the scheme and the run, a scheme that is not
	verified gets the notice of SchemeInfo.scheme_notice as a banner, and the table can be saved
	as TSV and the tree as SVG. Self-contained like the tree viewer: plain javascript, no
	server, no downloads'''
DASHBOARD_HTML_TEMPLATE = r'''<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>__TITLE__</title>
<style>
__BASECSS__
#side { width: 460px; flex: none; display: flex; flex-direction: column; min-height: 0;
        border-right: 1px solid var(--line); background: #fdfefe; }
#detail { padding: 10px 12px; border-bottom: 1px solid var(--line); font-size: 12px; min-height: 66px; max-height: 46%;
          overflow-y: auto; background: #fff; overflow-wrap: anywhere; word-break: break-word; }
#detail h2 { margin: 0 0 4px 0; font-size: 15px; font-weight: 600; }
#detail .line { margin: 2px 0; }
#detail .reason { color: #92400e; margin: 2px 0; }
#detail .links { margin: 4px 0 6px 0; }
#detail a { color: var(--accent); }
#detail .panelTitle { margin-top: 6px; }
#statusKey { padding: 6px 12px; border-bottom: 1px solid var(--line); font-size: 11.5px; color: var(--muted); }
#statusKey summary { cursor: pointer; font-weight: 600; font-size: 12px; }
#statusKey div { margin: 3px 0; }
#statusKey .badge { min-width: 66px; text-align: center; margin-right: 4px; }
#search { margin: 8px 12px; padding: 5px 9px; border: 1px solid #d4d9de; border-radius: 6px; font: inherit; font-size: 12px; }
#tableWrap { flex: 1; overflow-y: auto; }
table { border-collapse: collapse; width: 100%; font-size: 12px; table-layout: fixed; }
td.cut { overflow-wrap: anywhere; }
th { position: sticky; top: 0; background: #f4f6f8; text-align: left; padding: 5px 6px; font-weight: 600; border-bottom: 1px solid var(--line); }
td { padding: 4px 6px; border-bottom: 1px solid #eef1f3; vertical-align: top; }
td.num, th.num { text-align: right; font-variant-numeric: tabular-nums; }
tr.row { cursor: pointer; }
tr.row:hover { background: #f1f4f7; }
tr.selected { background: #e8f0fe; }
.badge { display: inline-block; border-radius: 999px; padding: 1px 7px; font-size: 11px; color: #fff; }
.sev-failed, .sev-notcalled { background: #b03a2e; }
.sev-review { background: #b45309; }
.sev-clean { background: #2e8b57; }
.sev-reused { background: #8a9099; }
.warn { color: #b45309; font-weight: 600; }
@media print {
  #side { width: auto; border: none; display: block; }
  #detail { max-height: none; border: none; padding: 0 0 6px 0; }
  #statusKey { border: none; padding: 0 0 6px 0; }
  #tableWrap { overflow: visible; }
  th { position: static; background: none; }
  tr.selected { background: none; }
  .badge { color: #000 !important; background: none !important; border: 1px solid #000; }
  #detail a { color: #000; text-decoration: none; }
}
</style>
</head>
<body>
<div id="head">
<h1 id="title">CanSNPer2 results</h1>
<div class="meta" id="counts"></div>
<div id="prov"></div>
<div id="notice"></div>
</div>
<div id="bar">
<div class="group">
<button id="overview" title="Back to the starting collapsed overview">Overview</button>
<button id="expandAll">Expand all</button>
<button id="collapseAll">Collapse</button>
<button id="fit" title="Fit the whole tree into the window">Fit</button>
<button id="zoomIn" title="Zoom in">+</button>
<button id="zoomOut" title="Zoom out">&minus;</button>
</div>
<span class="sep"></span>
<div class="group">
<label class="check" title="Show only the branches that lead to derived SNPs"><input type="checkbox" id="onlyCalled"><span>SNP paths only</span></label>
<label class="check"><input type="checkbox" id="showLabels" checked><span>Labels</span></label>
<button id="saveSvg" title="Save the tree as it is drawn now as an SVG file">Save SVG</button>
<button id="downloadTsv" title="Save the sample table, with the run information as # lines">Download table (TSV)</button>
</div>
<div id="legend"></div>
</div>
<div id="main">
<div id="side">
<div id="detail"></div>
<details id="statusKey" open><summary>Status key</summary>
<div><span class="badge sev-clean">clean</span>SNP path confirmed; every SNP on it has alignment data, none is near an alignment edge, and no derived SNP lies outside it.</div>
<div><span class="badge sev-review">review</span>SNP path confirmed, but SNPs on it have no alignment data or lie near an alignment edge, or derived SNPs lie outside it.</div>
<div><span class="badge sev-notcalled">not called</span>No SNP path met min_required_hits and strictness, no derived SNP was found, or two SNP paths competed at the same depth.</div>
<div><span class="badge sev-failed">failed</span>The alignment or the processing of the sample failed.</div>
<div><span class="badge sev-reused">reused</span>Result read from an earlier run; --rerun computes it again and shows it in the tree.</div>
</details>
<input id="search" type="text" placeholder="Filter samples by name">
<div id="tableWrap"><table><colgroup><col><col style="width: 25%"><col style="width: 60px"><col style="width: 58px"><col style="width: 100px"></colgroup><thead><tr><th>Sample</th><th>Final SNP</th><th class="num" title="Derived SNPs found on the SNP path / SNPs on the SNP path">Support</th><th class="num" title="SNP positions on the SNP path without alignment data">No data</th><th>Status</th></tr></thead><tbody id="rows"></tbody></table></div>
</div>
<div id="tree"></div>
</div>
<script>
var DATA = __PAYLOAD__;
window.onerror = function(message, source, line){
	var banner = document.createElement("div");
	banner.style.cssText = "background:#c0392b;color:#fff;padding:10px 14px;font-size:13px;";
	banner.textContent = "The CanSNPer2 results page hit a script error: " + message + " (line " + line + "). Please report this message together with the CanSNPer2 version.";
	document.body.insertBefore(banner, document.body.firstChild);
};
__TREEJS__
var samples = DATA.samples;
var SEVRANK = {"failed": 0, "not called": 1, "review": 2, "reused": 3, "clean": 4};
var SEVCLS = {"failed": "sev-failed", "not called": "sev-notcalled", "review": "sev-review", "clean": "sev-clean", "reused": "sev-reused"};
var TAB = String.fromCharCode(9), NL = String.fromCharCode(10);
var order = samples.map(function(s, i){ return i; });
order.sort(function(a, b){
	var sa = samples[a], sb = samples[b];
	var ra = SEVRANK[sa.severity], rb = SEVRANK[sb.severity];
	if (ra !== rb) { return ra - rb; }
	return sa.name < sb.name ? -1 : 1;
});
var counts = {};
samples.forEach(function(s){ counts[s.severity] = (counts[s.severity] || 0) + 1; });
document.getElementById("title").textContent = "CanSNPer2 results" + (INFO.scheme ? ": " + INFO.scheme : "");
document.getElementById("counts").innerHTML = samples.length + " samples: "
	+ "<b>" + (counts["clean"] || 0) + "</b> clean, "
	+ "<b>" + (counts["review"] || 0) + "</b> review, "
	+ "<b>" + (counts["not called"] || 0) + "</b> not called, "
	+ "<b>" + (counts["failed"] || 0) + "</b> failed, "
	+ "<b>" + (counts["reused"] || 0) + "</b> reused";
provHtml();
legendHtml([["derived", "Derived"], ["ancestral", "Ancestral"], ["other_base", "Other base"], ["no_data", "No alignment data"], ["untested", "Not tested"], ["path", "SNP path"], ["edge", "Near an alignment edge"]]);
function called(s){ return s.final && s.final !== "NA"; }
function support(s){ return called(s) && s.support ? s.support : "NA"; }
function noData(s){ return called(s) && s.path_no_data ? String(s.path_no_data.length) : "NA"; }
function pathOf(s){
	if (s.path_snps) { return s.path_snps; }
	return s.path && s.path !== "NA" ? s.path.split(";").map(function(p){ return {name: p}; }) : [];
}
function finalLabels(s){
	if (!called(s)) { return ""; }
	var onPath = {};
	pathOf(s).forEach(function(p){ onPath[p.name] = 1; });
	return ANN.filter(function(a){ return (a.kind === "node" && a.node === s.final) || (a.kind === "clade" && onPath[a.node]); }).map(function(a){
		var col = annColor(a);
		var style = a.kind === "clade" ? "background:transparent;border:1.5px dashed " + col + ";color:" + col : "background:" + col;
		return ' <span class="pill" style="' + style + '" title="' + esc((a.rank ? a.rank + ": " : "") + a.label) + '">' + esc(a.label) + "</span>";
	}).join("");
}
function fileLink(file, label){ return '<a href="' + encodeURIComponent(file) + '" title="' + esc(file) + '">' + esc(label) + "</a>"; }
var current = -1;
function svgName(){ return "cansnper2_" + (current >= 0 ? samples[current].name : "results") + "_tree.svg"; }
function paint(i){
	current = i;
	var s = samples[i];
	var states = s.states || {}, calledSet = {}, edge = {}, onPath = {};
	var hasStates = Object.keys(states).length > 0;
	(s.called || []).forEach(function(c){ calledSet[c] = 1; });
	(s.edge || []).forEach(function(c){ edge[c] = 1; });
	if (called(s)) { pathOf(s).forEach(function(p){ onPath[p.name] = 1; }); }
	var KIND = {1: "derived", 2: "ancestral", 3: "other_base"};
	nodes.forEach(function(n){
		n.called = !!calledSet[n.name];
		n.edge = !!edge[n.name];
		n.onPath = !!onPath[n.name];
		n.final = called(s) && n.name === s.final;
		if (isLabel(n)) { n.kind = "label"; n.dot = false; return; }
		n.dot = true;
		n.kind = KIND[states[n.name]] || (hasStates ? "no_data" : "untested");
	});
	markCalled();
	autoCollapse();
	overviewState = nodes.map(function(n){ return n.collapsed; });
	var reasons = (s.reasons || []).map(function(r){ return '<div class="reason">' + esc(r) + "</div>"; }).join("");
	var links = [];
	if (s.severity !== "failed" && INFO.snp_files !== false) {
		links.push(fileLink(s.name + "_snps.txt", "called SNPs"));
		links.push(fileLink(s.name + "_not_called.txt", "not called"));
	}
	if (s.severity !== "failed" && INFO.save_tree) { links.push(fileLink(s.name + "_tree.html", "tree")); }
	var line = "Final SNP <b>" + esc(s.final) + "</b>" + finalLabels(s);
	if (called(s)) { line += " &middot; Support " + esc(support(s)) + " &middot; No data " + esc(noData(s)); }
	var path = pathOf(s);
	var head = called(s) ? "SNP path, root to final SNP" : (path.length ? "Derived SNPs (no SNP path confirmed)" : "SNP path");
	document.getElementById("detail").innerHTML = "<h2>" + esc(s.name) + ' <span class="badge ' + SEVCLS[s.severity] + '">' + esc(s.severity) + "</span></h2>"
		+ '<div class="line">' + line + "</div>" + reasons
		+ (links.length ? '<div class="links">' + links.join(" &middot; ") + "</div>" : "")
		+ '<div class="panelTitle">' + head + "</div>" + pathHtml(path, s.final);
	var rows = document.querySelectorAll("#rows tr");
	for (var k = 0; k < rows.length; k++) { rows[k].className = rows[k].className.replace(" selected", ""); }
	var sel = document.querySelector("#rows tr[data-i='" + i + "']");
	if (sel) { sel.className += " selected"; }
	render();
	fitInitial();
}
wirePath(document.getElementById("detail"));
function buildTable(filter){
	var rows = document.getElementById("rows");
	rows.innerHTML = "";
	order.forEach(function(i){
		var s = samples[i];
		if (filter && s.name.toLowerCase().indexOf(filter) === -1) { return; }
		var tr = document.createElement("tr");
		tr.className = "row" + (i === current ? " selected" : "");
		tr.setAttribute("data-i", i);
		var note = (s.reasons || []).length ? ' <span class="warn" title="' + esc(s.reasons.join("; ")) + '">&#9888; ' + s.reasons.length + "</span>" : "";
		tr.innerHTML = '<td class="cut" title="' + esc(s.name) + '">' + esc(s.name) + '</td><td class="cut" title="' + esc(s.final) + '">' + esc(s.final) + "</td>"
			+ '<td class="num">' + esc(support(s)) + '</td><td class="num">' + esc(noData(s)) + "</td>"
			+ '<td><span class="badge ' + SEVCLS[s.severity] + '">' + esc(s.severity) + "</span>" + note + "</td>";
		tr.onclick = function(){ paint(i); };
		rows.appendChild(tr);
	});
}
function tsvText(){
	var clean = new RegExp("[" + TAB + NL + String.fromCharCode(13) + "]", "g");
	function cell(v){ return String(v === undefined || v === null ? "" : v).replace(clean, " "); }
	var out = (DATA.provenance || []).map(function(l){ return "# " + cell(l); });
	out.push(["sample", "final_snp", "support", "no_data", "status", "snp_path", "notes"].join(TAB));
	order.forEach(function(i){
		var s = samples[i];
		out.push([s.name, s.final, support(s), noData(s), s.severity, s.path, (s.reasons || []).join("; ")].map(cell).join(TAB));
	});
	return out.join(NL) + NL;
}
document.getElementById("downloadTsv").onclick = function(){ download("cansnper2_results" + (INFO.scheme ? "_" + INFO.scheme : "") + ".tsv", tsvText(), "text/tab-separated-values"); };
document.getElementById("search").oninput = function(){ buildTable(this.value.trim().toLowerCase()); };
buildTable("");
paint(order[0]);
if (document.querySelectorAll("#rows tr").length !== samples.length || !svgEl) {
	window.onerror("the sample table or the tree stayed empty", "", 0);
}
</script>
</body>
</html>
'''

NEWICK_SPECIAL = re.compile(r"[\s(),:;'\[\]]")


def quote_newick(name):
	'''A node name as it is written in a newick tree: a name with whitespace or one of ( ) , : ; [ ]
		or a single quote is put in single quotes, and a quote inside the name is doubled.'''
	if NEWICK_SPECIAL.search(name):
		return "'" + name.replace("'", "''") + "'"
	return name


class NewickNode(object):
	"""The NewickNode class stores the information of a taxonomy node
			ID
			name
			children
			parent
		The purpose of this class is to allow a fast printout of a newick tree.
		All nodes in a newick tree knows it's children and parent, so by selecting a
		print option for the class objects can inheritly print all its children or the lineage (parents).

		main functions
			add_child  ## Add a newick node object reference as a child of the node
			set_print  ## Set the print type  (name, id, lineage or newick <- default)
	"""

	def __init__(self, id, name, parent=False):
		super(NewickNode, self).__init__()
		self.id         = id            	## Node id
		self.name       = name          	## Node name
		self.parent     = parent        	## Parent id False for root
		self.children   = set()         	## Set of newick children
		self.__class__.print_opt = "newick" ## The default behaviour of this class is to print out a
											## 		newick tree from the given node (as root)


	def __str__(self):
		'''The print function of NewickNode allows any node to print all its children in newick format or the lineage to root
			There is also an option to print just the name of a node
		'''
		if self.__class__.print_opt == "name":
			return "{name}".format(name=self.name)
		elif self.__class__.print_opt == "lineage":
			if self.parent:
				return "{parent};{name}".format(name=self.name, parent=self.parent)
			else:
				return "root"
		elif self.__class__.print_opt == "newick":
			return self.newick()
		else:
			return "{id}: {name}; children: {nchildren} ".format(id=self.id, name=self.name, nchildren = len(self.children))

	def __repr__(self):
		return "NewickNode()"

	def newick(self):
		'''The subtree of the node in newick format, written without recursion: a tree of several
			thousand levels exceeded the recursion limit. Children are rendered in name order: the
			children container is a set, and its iteration order changes between processes with the
			python hash seed, which made the newick export (and with it the drawn trees and the
			release gate's check on them) different from one run to the next. Names with characters
			that are special in newick are written in single quotes.'''
		out = []
		stack = [self]
		while stack:
			item = stack.pop()
			if isinstance(item, str):
				out.append(item)
				continue
			children = sorted(item.children, key=lambda child: child.name)
			if not item.parent:  #Only the root node has this property
				tail = ")" + quote_newick(item.name) + ")ROOT;"
				head = "(("
			elif children:
				tail = ")" + quote_newick(item.name)
				head = "("
			else:
				out.append(quote_newick(item.name))
				continue
			stack.append(tail)
			for position, child in enumerate(reversed(children)):
				if position:
					stack.append(",")
				stack.append(child)
			out.append(head)
		return "".join(out)

	def add_child(self,child):
		'''Add a NewickNode object as child'''
		self.children.add(child)
		return

	def set_print(self,_type):
		'''Set the variable that controls the print style
			name
			newick - the complete subtree in newick format
			lineage
		'''
		self.__class__.print_opt = _type


class NewickTree(object):
	"""The NewickTree class parses a CanSNP or a FlexTaxD database and prints the database as a newick tree

		print styles
			name
			newick - the complete subtree in newick format
			lineage - prints all parents of a node up to root

	"""

	def __init__(self, database,name="newick",outdir="./",min_required_hits=3,strictness=0.7):
		super(NewickTree, self).__init__()
		self.database = CanSNPdbFunctions(database) ## Initiate database connection with CanSNPdbFunctions
		self.database_path = database
		self.nodeDict = {}							## Dictionary to store references to all newick nodes
		self.c_p_set = set()
		self.tree_file = "{outdir}/{name}_tree.html".format(outdir=outdir.rstrip("/"),name=name) ## output file
		self.min_required_hits = min_required_hits
		self.strictness = strictness
		## Build the newick tree
		self.newickTree = str(self.build_tree())
		## Tree colors
		self.snp_colors = {
			"derived": "#2e8b57", 			## green, also the line of the SNP path
			"ancestral": "#a3a3a3", 		## neutral grey
			"non_aligned": "#ffffff", 		## open circle (white fill) ...
			"outline": "#6b7280", 			## ... with a grey outline: no alignment data
			"other_base": "#2166ac",		## blue
			"derived_text": "#1f6b43"		## the label of the final SNP
		}

	def __repr__(self):
		return "NewickTree()"

	def close(self):
		'''close the database connection; one tree is made per sample, and the connections stayed
			open until the garbage collector freed them'''
		try:
			self.database.disconnect()
		except Exception:
			pass

	def print(self):
		print(self.newickTree)

	def get_tree(self,table="tree"):
		'''Function that returns the whole tree in the database the script expects
			the tree to be rooted at the lowest value and that the root has itself as parent'''
		SELECT = "SELECT parent,child FROM {table} ORDER BY child ASC".format(table=table)
		return self.database.query(SELECT).fetchall()

	def get_nodes(self,table="nodes"):
		'''Function that returns all nodes in the database'''
		nt = {}
		SELECT = "SELECT id,name,snp_id FROM {table} LEFT JOIN snp_annotation on (snp_annotation.node_id = {table}.id)".format(table=table)
		for id,name,snp_id in self.database.query(SELECT).fetchall():
			if snp_id != None:
				nt[id] = snp_id
			else:
				'''the tree is read once per sample and again for the dashboard: warn once per node'''
				if name not in _WARNED_NODES:
					_WARNED_NODES.add(name)
					logger.warning("#Node has no snp_annotation {node}".format(node=name))
				nt[id] = "-"+name+"-"
		return nt

	def get_parent(self,name):
		'''return parent'''
		#QUERY = '''SELECT parent,child,rank FROM tree LEFT JOIN rank on (tree.rank_i = rank.rank_i) WHERE child = "{node}"'''.format(node=name)
		QUERY = '''SELECT parent,child,rank_i FROM tree WHERE child = ?'''
		res = self.database.query(QUERY, insert_val=(name,), getres=True).fetchone()
		return res

	def get_child(self,name):
		'''return child'''
		#QUERY = '''SELECT child,child,rank FROM tree LEFT JOIN rank on (tree.rank_i = rank.rank_i) WHERE child = "{node}"'''.format(node=name)
		QUERY = '''SELECT parent,child,rank_i FROM tree WHERE parent = ?'''
		res = self.database.query(QUERY, insert_val=(name,), getres=True).fetchone()
		return res

	def _attach(self,child,nodes,parent):
		'''Add the node of child below the node of parent, raises KeyError when the parent has no node yet'''
		if len(self.c_p_set & set([(child,parent)])) == 0:
			node = NewickNode(child, nodes[child], self.nodeDict[parent])  		## this works also for root as root has itself as child
			'''Make sure link to parent was not made before'''
			self.c_p_set |= set([(child,parent)])
			self.c_p_set |= set([(parent,child)])
			self.nodeDict[child] = node
		else:
			logger.debug("Link {child}-{parent} already exists, retrieve node!".format(child=child,parent=parent))
			node = self.nodeDict[child]										## add a reference to the node so that children can be added
		pnode = self.nodeDict[parent]										## add child to the parent node
		pnode.add_child(node)
		logger.debug("NewickNode p:{parent} c: {child} added".format(parent=parent,child=child))

	def new_node(self,child,nodes,parent):
		'''Function that adds a new node to the newick tree. A parent that has no node yet is added
			first, walking up the tree table without recursion (a tree of thousands of levels whose
			ids run against the order of the levels exceeded the recursion limit)'''
		pending = [(child,parent)]
		waiting = set([child])
		while pending:
			current,current_parent = pending[-1]
			try:
				self._attach(current,nodes,current_parent)
				pending.pop()
				continue
			except KeyError:
				logger.debug("Error in adding NewickNode parent: {parent} does not exist, trying to add parent".format(parent=current_parent))
			parent_link = self.get_parent(current_parent)
			if parent_link is None:	## parent has no link in the tree table, node can not be attached
				logger.warning("Not adding node {child}: parent {parent} does not exist and has no parent link in the tree table".format(child=current,parent=current_parent))
				return
			t_parent,t_child,rank = parent_link
			if t_child in waiting:
				logger.warning("Not adding node {child}: the tree table links {node} to itself through its parents".format(child=current,node=t_child))
				return
			logger.debug("Adding parent: {parent} of parent: {child}".format(parent=t_parent,child=t_child))
			waiting.add(t_child)
			pending.append((t_child,t_parent))

	def build_tree(self):
		'''Build newick tree from database
			This function walks through a database of nodes and creates NewickNode objects
			self-aware of their decending newick tree or their parent lineage,
			Returns: The root of the tree, however all nodes are accesible from the
						NewickTree nodeDict by their node name
		'''
		tree = self.get_tree()
		nodes = self.get_nodes()
		logger.debug("Nodes: {n} Links: {l}".format(n=len(nodes),l=len(tree)))
		logger.debug([nodes,tree])
		newickTree = None
		for parent,child in tree:
			if parent == child:  ## root
				root = NewickNode(child, nodes[child], False)					## Create the root node
				self.nodeDict[child] = root										## Add the root node to the node dictionary
				self.nodeDict["root"] = root									## Also add this reference as "root"
				newickTree = root 												## The master parent will contain the full tree
				continue
			self.new_node(child,nodes,parent)
		if newickTree is None:
			'''without the initialization a database with an empty tree table died with an
				UnboundLocalError here; the clear error lets callers decide whether the tree is
				needed at all'''
			raise ValueError("the tree table of the database contains no root row (parent == child), the tree cannot be built")
		## The newickTree is the same as the master parent (containing the full tree, the root node is defined on row 135)
		self.root = newickTree ## kept as an object so draw_tree can walk the tree without re-parsing the string
		logger.debug("Tree complete, return newickTree")
		#logger.debug(newickTree)
		return newickTree

	def _tree_entries(self):
		'''Flatten the NewickNode tree into a list of dicts in preorder, shaped exactly like the
			ete3 tree used to be parsed from the newick string: a "ROOT" wrapper above the
			database root (NewickNode.__str__ renders the root as ")ROOT;"), every entry carrying
			id, name, parent id and depth. Children are walked in name order as in __str__, so the
			first entry of a name is the node search_nodes(name)[0] used to match (issue #78).'''
		entries = []
		def add(name, parent):
			entry = {"id": len(entries),
					 "name": name,
					 "parent": parent,
					 "depth": 0 if parent is None else entries[parent]["depth"] + 1}
			entries.append(entry)
			return entry["id"]
		root_id = add("ROOT", None)
		stack = [(self.root, root_id)]
		while stack:  ## preorder without recursion, children in name order
			node, parent = stack.pop()
			entry_id = add(node.name, parent)
			for child in sorted(node.children, key=lambda child: child.name, reverse=True):
				stack.append((child, entry_id))
		return entries

	def _parent_snp(self,entry,entries):
		'''name of the nearest SNP above the entry, skipping "-name-" label nodes; None at the top'''
		entry = entries[entry["parent"]]
		while entry["name"] != "ROOT":
			if not (entry["name"].startswith("-") and entry["name"].endswith("-")):
				return entry["name"]
			entry = entries[entry["parent"]]
		return None

	def _chain_names(self,entry,entries):
		'''the SNP names from the root of the tree down to the entry, the entry included; "-name-" label
			nodes and the ROOT wrapper are not SNPs and are left out'''
		names = []
		while entry["name"] != "ROOT":
			name = entry["name"]
			if not (name.startswith("-") and name.endswith("-")):
				names.append(name)
			entry = entries[entry["parent"]]
		return names[::-1]

	def _confirm_path(self,entry,entries,called_snps,snplist):
		'''Confirm the ACTUAL ancestor chain of the node, from the root of the tree down to the
			node itself (issue #25). Every SNP on the chain is one of:
			  called    derived in the query, counts as a hit;
			  no data   no alignment covered the position, or the query has a gap or N there
			            (state 0): there is no evidence either way, so it is a gap in the path;
			  against   ancestral (2) or a third base (3): the walk down the tree should have left
			            the branch here, so the chain is broken and the node cannot be called.
			"-name-" label nodes (a database root without snp_annotation) are not SNPs and are
			skipped. The node is confirmed when the chain is not broken, at least min_required_hits
			SNPs are called, the called fraction of the SNPs on the chain is at least strictness,
			and the SNP directly above the node is called: without that last rule one stray derived
			SNP on a neighbouring branch, whose parent had no data, outranked the real lineage.
			A single unaligned SNP near the root used to break every chain below it, so the sample
			was NA however many SNPs on its path were derived, and --strictness had no effect.
			Returns (confirmed, called, snps_on_chain, [SNPs without data]).'''
		logger.debug("Confirm path of {node}".format(node=entry["name"]))
		start = entry
		count = 0
		on_chain = 0
		missing = []
		called = set(called_snps)
		while entry["name"] != "ROOT":
			name = entry["name"]
			entry = entries[entry["parent"]]
			if name.startswith("-") and name.endswith("-"):
				continue
			on_chain += 1
			state = snplist.get(name,0)
			if name in called and state != 2:
				count +=1
			elif state in (2,3):
				logger.debug("SNP path broken at {node} ({state}), {count} called SNPs below it".format(
						node=name,state="ancestral" if state == 2 else "other base",count=count))
				return False,count,on_chain,missing
			else:
				missing.append(name)
		quota = float(count)/on_chain if on_chain else 0
		'''the SNP right above the node is the second SNP of the chain, counted from the node'''
		parent_missing = on_chain > 1 and self._parent_snp(start,entries) in missing
		if parent_missing:
			logger.debug("Path not confirmed: the SNP directly above {node} has no data".format(node=start["name"]))
			return False,count,on_chain,missing
		if count >= self.min_required_hits and quota >= self.strictness:
			logger.info("Confirm path: {count}/{n} SNPs called (strictness {strictness}, min required hits {minhit}){gaps}".format(
					count=count,n=on_chain,strictness=self.strictness,minhit=self.min_required_hits,
					gaps=", without data: "+", ".join(missing) if missing else ""))
			return True,count,on_chain,missing
		logger.debug("Path not confirmed: {count}/{n} SNPs called (strictness {strictness}, min required hits {minhit})".format(
				count=count,n=on_chain,strictness=self.strictness,minhit=self.min_required_hits))
		return False,count,on_chain,missing

	def draw_tree(self,snplist,called_snps=False,save_tree=True,summary=False,annotations=None,snp_info=None,run_info=None):
		'''Colors the SNP tree with the states found in the query alignments, calls the final SNP
			(_final_call) and writes the self-contained html viewer described in the module header.
		Keyword arguments:
		snplist -- a dict of the SNP names and their state
		snp_info -- the SNP_info rows of the sample, for the target base of each SNP on the path
		run_info -- what the page records about the run (CanSNPer2.run_info)
		'''
		logger.debug("Draw tree from snplist")
		entries = self._tree_entries()
		'''nodes of the tree that no alignment carried any base for: they have to be reported,
			not left out of the viewer and out of both result files (issue #27)'''
		no_data_nodes = set()
		called_set = set(called_snps) if called_snps else set()
		for entry in entries:
			n = entry["name"]
			entry["fill"] = self.snp_colors["non_aligned"]
			entry["r"] = 5           ## ete3 "size" was a diameter: 10 drew as r=5, keep that scale
			entry["line"] = "#DDDDDD"
			entry["dash"] = True
			entry["dot"] = True
			entry["kind"] = "untested"
			entry["called"] = n in called_set
			if n == "ROOT" or (n.startswith("-") and n.endswith("-")):
				'''ROOT and the "-name-" labels that get_nodes gives to nodes without any
					snp_annotation (the root of a database like that) are not SNPs: they are only
					meant to be line junctions and a label node must not be drawn as a grey SNP
					dot (issue #52, with #32 and #34)'''
				entry["dot"] = False
				entry["kind"] = "label"
				continue
			if len(snplist) > 0:
				snpvalue = snplist.get(n)
				entry["kind"] = SNP_KINDS.get(snpvalue, "no_data")
				if snpvalue in (1, 2, 3):
					entry["fill"] = self.snp_colors[entry["kind"]]
					entry["line"] = "#000000"
					entry["dash"] = False
					if snpvalue != 2:
						entry["r"] = 7.5
				elif n not in snplist:
					'''SNP not in list: no alignment covered this position at all. The node used
						to be skipped with a continue, which also skipped the style and left it at
						the ete default size, invisible in the pdf (issue #27)'''
					no_data_nodes.add(n)
			elif called_snps:
				### Only care about called SNPs (the final SNPs of the summary tree) and color them
				if entry["called"]:
					entry["fill"] = self.snp_colors["derived"]
					entry["r"] = 7.5
					entry["kind"] = "derived"
				entry["line"] = "#000000"
				entry["dash"] = False
		'''the nodes without data are drawn as open circles, and their database annotation is
			returned as well, so the per sample files can list them instead of hiding that they
			were never compared at all (issue #27)'''
		no_data_info = []
		if no_data_nodes:
			logger.warning("{n} SNP positions in the tree had no alignment data at all, no reference block covered them or --mask dropped them: {ids}".format(
					n=len(no_data_nodes), ids=", ".join(sorted(no_data_nodes)[:20])))
			no_data_info = self.get_no_data_info(no_data_nodes)
		if summary:
			if save_tree:
				self.write_tree_html(entries, annotations, page=self._summary_page(called_snps), run_info=run_info)
			return
		result = self._final_call(entries, snplist, called_snps, no_data_info)
		## Write the tree viewer
		if save_tree:
			self.write_tree_html(entries, annotations, page=self._sample_page(result, snplist, called_snps, snp_info), run_info=run_info)
		return result

	def _summary_page(self, final_snps):
		'''the header of the summary tree: the final SNPs of the run'''
		return {"sample": "summary", "summary": True, "finals": sorted(final_snps or [])}

	def _sample_page(self, result, snplist, called_snps, snp_info):
		'''the header of a sample tree: the final SNP, the support of the call and the SNP path
			(as in <sample>_snps.txt) with the bases of each SNP on it'''
		final_snp, _msg, called, _no_data = result
		called_names = [snp[1] for snp in called] if called else []
		path = called_names
		final, support, missing = "NA", "NA", []
		if final_snp and final_snp[2][0] and len(final_snp) > 5:
			path = [name for name in final_snp[5] if name in set(called_names)]
			final = final_snp[1]
			support = "{found}/{n}".format(found=final_snp[2][1], n=final_snp[3])
			missing = list(final_snp[4])
		return {"sample": os.path.basename(self.tree_file)[:-len("_tree.html")], "summary": False,
			"final": final, "support": support, "no_data": len(missing) if final != "NA" else "NA",
			"path": path_details(path, snp_info, snplist)}

	def _final_call(self, entries, snplist, called_snps, no_data_info):
		'''The final SNP call of the sample from the colored tree entries (see draw_tree)'''
		dlist = []
		logger.info("Called snps: {called}".format(called=called_snps))
		if called_snps:
			try:
				node_index = {} ## first entry of every name, preorder, like search_nodes()[0]
				for entry in entries:
					node_index.setdefault(entry["name"], entry["id"])
				'''Calculate the distance from all nodes to the root: every branch has unit length,
					so the ete distance to the root is just the depth of the node'''
				for tsnp in called_snps:
					if tsnp in node_index:
						dlist.append([entries[node_index[tsnp]]["depth"],tsnp])
					elif tsnp == "T/N.1":
						dlist.append([0,tsnp])
					else:
						logger.debug("Distance could non be calculated for {snp}".format(snp=tsnp))
				## Sort the list of distances from deepest
				dlist = sorted(dlist,key=lambda l:l[0], reverse=True)
				if len(dlist) < self.min_required_hits: ## the number of SNPs found is too small to be valid,
					msg = "The number of SNPs found is too small (min {min}), Cannot call SNPs.".format(min=self.min_required_hits)
					logger.debug(msg)
					return False,msg,[],no_data_info
				'''The final SNP is the deepest called SNP whose ACTUAL ancestor chain in the tree, from the
					root down to the node, is confirmed by _confirm_path: no ancestral or third base on it,
					and enough of its SNPs called (min_required_hits, strictness). The old walk accepted SNPs
					from branches it should have left: _check_start picked a start node by depth adjacency
					in the flat sorted list and _confirm_path counted any called ancestor, so a node on
					another branch could be final-called (issue #25).
					When confirmed paths compete at the same depth no path is auto-picked: the sample stays NA
					and the candidates are reported, leaving room for the calling policy of issue #47.'''
				logger.debug("Calling SNP")
				candidates = [] ## confirmed paths, all at the greatest confirmed depth
				deepest = None ## deepest called SNP found in the tree, kept for the "not confirmed" report
				for dist,tsnp in dlist: ## loop the called SNPs from deepest
					if tsnp not in node_index: ## the node name from the db is not in the tree, cannot call (issue #78)
						logger.debug("Called SNP {snp} was not found in the tree, cannot confirm it".format(snp=tsnp))
						continue
					ok,count,on_chain,missing = self._confirm_path(entries[node_index[tsnp]],entries,called_snps,snplist)
					if not deepest:
						deepest = [dist,tsnp,(False,count)]
					if not ok:
						continue
					if not candidates:
						candidates = [[dist,tsnp,count,on_chain,missing]]
					elif dist == candidates[0][0]:
						candidates.append([dist,tsnp,count,on_chain,missing])
					else: ## a shallower confirmed path cannot beat the deepest confirmed one
						break
				dlist = sorted(dlist,key=lambda l:l[0], reverse=False)
				if len(candidates) == 1:
					dist,tsnp,count,on_chain,missing = candidates[0]
					'''on_chain and missing are appended for the results dashboard: a call that
						needed the strictness tolerance over positions without data must stand
						out for manual review; the sixth item is the SNP chain from the root to the final SNP;
						consumers that index 0-2 are unaffected'''
					called = [dist,tsnp,(True,count),on_chain,missing,self._chain_names(entries[node_index[tsnp]],entries)]
					msg = "Final SNP called: {snp} ({count}/{n} SNPs on the path called)".format(snp=tsnp,count=count,n=on_chain)
					if missing:
						'''a call that relies on the strictness tolerance says which SNPs had no data'''
						msg += ", {k} without alignment data: {ids}".format(k=len(missing),ids=", ".join(missing))
					logger.info(msg)
				elif len(candidates) > 1:
					paths = ", ".join(["{snp} ({n} SNPs)".format(snp=c[1],n=c[2]) for c in candidates])
					msg = "Two or more confirmed SNP paths compete at the same depth: {paths}, no path is auto-picked, check {tree}".format(paths=paths,tree=self.tree_file)
					logger.warning(msg)
					called = deepest
				elif deepest:
					msg = "SNP could not be confirmed!"
					called = deepest
				else:
					msg = "None of the called SNPs was found in the tree, cannot call SNPs."
					logger.warning(msg)
					return False,msg,[],no_data_info
				logger.debug(dlist)
			except KeyError as e:
				logger.error("KeyError, {e}".format(e=e))
				msg = "called_snps function failed, KeyError {e}".format(e=e)
				logger.debug(msg)
				return False,msg,[],no_data_info
		else:
			msg = "No SNPs were called."
			logger.debug(msg)
			return False,msg,[],no_data_info
		return called,msg,dlist,no_data_info

	def _check_annotations(self, annotations, entries):
		'''drop annotations whose SNP name is not in the tree, loudly: an annotation file
			written against another database would otherwise silently draw nothing'''
		if not annotations:
			return []
		names = set(entry["name"] for entry in entries)
		kept = []
		unknown = []
		for annotation in annotations:
			if annotation["node"] in names:
				kept.append(annotation)
			else:
				unknown.append(annotation["node"])
		if unknown:
			logger.warning("{n} tree annotations name SNPs that are not in this database and were skipped: {ids}".format(
					n=len(unknown), ids=", ".join(sorted(set(unknown)))))
		return kept

	def provenance(self, run_info=None):
		'''what the pages record about the database (database_provenance) and the run (run_info,
			CanSNPer2.run_info; optional, so other callers can leave it out)'''
		info = database_provenance(self.database_path)
		info.update(run_info or {})
		return info

	def write_tree_html(self,entries,annotations=None,page=None,run_info=None):
		'''Write the self-contained viewer: the colored tree and its data are embedded in a single
			html file, no server, no javascript downloads, opens in any browser on any OS'''
		annotations = self._check_annotations(annotations or [], entries)
		page = dict(page or {"sample": os.path.basename(self.tree_file)[:-len("_tree.html")], "path": []})
		info = self.provenance(run_info)
		payload = json.dumps({"colors": self.snp_colors, "nodes": entries, "annotations": annotations,
			"page": page, "info": info, "provenance": provenance_lines(info)}, separators=(",", ":"))
		'''The payload sits inside a <script> block, and json.dumps leaves "<" as it is, so a node
			name from the database such as "B.1</script><script>..." closed the block and the rest
			ran as script when the tree was opened. "<" only occurs inside json strings, where
			\\u003c is the same character, so the data the viewer reads is unchanged.
			The title names the sample and the scheme and is escaped as html text. Both
			placeholders are filled in one pass, so a placeholder inside a name is not expanded.'''
		title = "CanSNPer2 tree: {sample}".format(sample=page.get("sample") or os.path.basename(self.tree_file))
		if info.get("scheme"):
			title += " ({scheme})".format(scheme=info["scheme"])
		values = {"TITLE": html.escape(title),
				  "PAYLOAD": payload.replace("<", "\\u003c")}
		page_html = re.sub(r"__(TITLE|PAYLOAD)__", lambda match: values[match.group(1)], render_template(TREE_HTML_TEMPLATE))
		with open(self.tree_file, "w") as out:
			out.write(page_html)

	def write_dashboard(self, reports, dashboard_file, annotations=None, run_info=None):
		'''Write the results overview (DASHBOARD_HTML_TEMPLATE): the sample table and the
			database tree embedded once, recolored per sample on click. Only the topology and
			the per-sample SNP states are stored, so the file stays small even for 100+ samples.
			run_info (CanSNPer2.run_info) is shown in the header; it is optional'''
		entries = self._tree_entries()
		annotations = self._check_annotations(annotations or [], entries)
		nodes = [{"id": entry["id"], "name": entry["name"], "parent": entry["parent"], "depth": entry["depth"]} for entry in entries]
		info = self.provenance(run_info)
		payload = json.dumps({"colors": self.snp_colors, "nodes": nodes, "samples": reports, "annotations": annotations,
			"info": info, "provenance": provenance_lines(info)}, separators=(",", ":"))
		'''same "<" escaping as write_tree_html: a database name must not close the script block'''
		title = "CanSNPer2 results"
		if info.get("scheme"):
			title += " ({scheme})".format(scheme=info["scheme"])
		values = {"TITLE": html.escape(title), "PAYLOAD": payload.replace("<", "\\u003c")}
		page = re.sub(r"__(TITLE|PAYLOAD)__", lambda match: values[match.group(1)], render_template(DASHBOARD_HTML_TEMPLATE))
		with open(dashboard_file, "w") as out:
			out.write(page)
		return dashboard_file

	def get_no_data_info(self,snp_ids):
		'''Look up the database annotation of the SNP nodes that no alignment carried a base for, in the
			same columns as the SNP_info rows (Name, Reference, Sequence, Pos, QuerySeq, QueryPosition, Ancestral
			base, Derived base, Target base, Edge flag), so the per sample output can list the positions
			that were never compared instead of leaving them out of both result files (issue #27). The
			target base is "N/A": there is no base there at all, which is different from a base that was
			read and did not confirm the SNP; the query columns are N/A for the same reason (issue #23).'''
		names = sorted(set(snp_ids))
		if not names:
			return []
		rows = []
		'''sqlite limits how many variables one statement may bind, so the lookup runs in chunks'''
		for i in range(0, len(names), 400):
			chunk = names[i:i+400]
			QUERY = '''SELECT snp_annotation.snp_id, COALESCE(snp_references.genome, snp_annotation.reference),
							snp_annotation.sequence, snp_annotation.sequence_position, snp_annotation.ancestral_base,
							snp_annotation.derived_base
						FROM snp_annotation
						LEFT JOIN snp_references ON (snp_references.id = snp_annotation.genome_i)
						WHERE snp_annotation.snp_id IN ({marks})'''.format(marks=",".join(["?"]*len(chunk)))
			for snp_id,genome,sequence,position,ancestral,derived in self.database.query(QUERY,insert_val=tuple(chunk),getres=True).fetchall():
				'''the empty last column keeps the rows in the shape of SNP_FILE_COLUMNS, whose
					edge flag cannot be computed for a position that no alignment carried; the
					query columns are N/A for the same reason (issue #23)'''
				rows.append([str(snp_id),str(genome),str(sequence),str(position),"N/A","N/A",str(ancestral),str(derived),"N/A",""])
		found = set([row[0] for row in rows])
		for name in names:
			if name not in found:
				'''a tree node without any annotation row: still a node the alignments never covered,
					list it rather than hide it again'''
				rows.append([name,"-","-","-","N/A","N/A","-","-","N/A",""])
		return sorted(rows, key=lambda row: row[0])
