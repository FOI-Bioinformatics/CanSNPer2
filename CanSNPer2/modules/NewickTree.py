'''
Module to read and write newick trees

'''
from CanSNPer2.modules.DatabaseConnection import CanSNPdbFunctions
import html
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


__version__ = "2.0.1"
__author__ = "David Sundell"
__credits__ = ["David Sundell"]
__license__ = "GPLv3"
__maintainer__ = "FOI bioinformatics group"
__email__ = ["bioinformatics@foi.se", "david.sundell@foi.se"]
__date__ = "2019-06-07"
__status__ = "Production"
__partof__ = "CanSNPer2"

'''The viewer written by NewickTree.write_tree_html: the colored tree arrives as JSON in
	__PAYLOAD__, everything else is plain javascript drawing an svg, no libraries, no network,
	so the html file works the same on any machine (issue #184 removed the ete3/Qt dependency
	that produced the old pdfs). The view is interactive: clicking a dot collapses or expands
	its subtree, collapse/expand/filters relayout the tree and refit it into the window, the
	viewBox pans by dragging and zooms with the wheel, and trees with more than 60 leaves open
	pre-collapsed to an overview - the full depth is only expanded when it is asked for.
	Branches that contain called SNPs are never collapsed (auto or manual): the called path
	and every green dot are always visible in full, and "overview" returns to the starting
	auto-collapsed depth'''
TREE_HTML_TEMPLATE = '''<!doctype html>
<html>
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>__TITLE__</title>
<style>
html, body { height: 100%; margin: 0; background: #fff; }
body { font-family: system-ui, -apple-system, "Segoe UI", Roboto, Helvetica, Arial, sans-serif; }
#bar { position: fixed; top: 0; left: 0; right: 0; z-index: 10; display: flex; align-items: center;
       gap: 8px; flex-wrap: wrap; padding: 8px 14px; background: #fbfcfd;
       border-bottom: 1px solid #e3e6e9; box-shadow: 0 1px 2px rgba(15, 23, 42, 0.05); }
#tools, #filters { display: flex; align-items: center; gap: 6px; }
button { appearance: none; border: 1px solid #d4d9de; background: #fff; color: #1f2937;
         border-radius: 7px; padding: 5px 12px; font: inherit; font-size: 13px; cursor: pointer; }
button:hover { background: #f1f4f7; border-color: #c2c8cf; }
button:active { background: #e8ecf0; }
#zoomIn, #zoomOut { width: 32px; padding: 5px 0; font-size: 15px; }
.sep { width: 1px; height: 24px; background: #e3e6e9; margin: 0 4px; }
label.check { display: inline-flex; align-items: center; gap: 6px; font-size: 13px; color: #374151;
              cursor: pointer; padding: 4px 9px; border-radius: 7px; }
label.check:hover { background: #f1f4f7; }
input[type="checkbox"] { width: 14px; height: 14px; accent-color: #3a7d5c; }
#legend { display: flex; gap: 14px; align-items: center; font-size: 12px; color: #4b5563; margin-left: auto; }
#legend .key { display: flex; align-items: center; gap: 5px; }
.dot { width: 11px; height: 11px; border-radius: 50%; }
#tree { position: fixed; top: 52px; bottom: 0; left: 0; right: 0; }
svg { width: 100%; height: 100%; display: block; cursor: grab; }
circle.node { cursor: pointer; }
text { fill: #333; font-size: 11px; pointer-events: none; user-select: none; }
#hint { position: fixed; bottom: 10px; right: 12px; z-index: 10; color: #6b7280; font-size: 12px;
        background: rgba(255, 255, 255, 0.85); border: 1px solid #e5e7eb; border-radius: 999px;
        padding: 4px 12px; }
</style>
</head>
<body>
<div id="bar">
<div id="tools">
<button id="overview" title="back to the starting collapsed overview">overview</button>
<button id="expandAll">expand all</button>
<button id="collapseAll">collapse</button>
<button id="fit" title="fit the tree into the window">fit</button>
<button id="zoomIn" title="zoom in">+</button>
<button id="zoomOut" title="zoom out">&minus;</button>
</div>
<span class="sep"></span>
<div id="filters">
<label class="check"><input type="checkbox" id="onlyCalled"><span>called paths only</span></label>
<label class="check"><input type="checkbox" id="showLabels" checked><span>labels</span></label>
</div>
<span class="sep"></span>
<div id="legend"></div>
</div>
<div id="tree"></div>
<div id="hint">click a dot: collapse / expand &middot; drag: pan &middot; wheel: zoom &middot; double click: fit &middot; called paths stay open</div>
<script>
var DATA = __PAYLOAD__;
var COL = 130, ROW = 22, MARGIN = 24, OVERVIEW_LEAVES = 60;
var nodes = DATA.nodes;
var children = {};
nodes.forEach(function(n){ if (n.parent !== null) { (children[n.parent] = children[n.parent] || []).push(n); } });
nodes.forEach(function(n){ n.collapsed = false; });
var LEGEND = [["derived","derived (present)"],["ancestral","ancestral (not present)"],["other_base","other base"],["non_aligned","no alignment data / untested"]];
document.getElementById("legend").innerHTML = LEGEND.map(function(l){ return '<span class="key"><span class=\\"dot\\" style=\\"background:' + DATA.colors[l[0]] + '"></span>' + l[1] + "</span>"; }).join("");
function hasCalled(n){ if (n.called) { return true; } return (children[n.id] || []).some(hasCalled); }
function keep(n){ return !onlyCalled || hasCalled(n); }
function visKids(n){ if (n.collapsed) { return []; } return (children[n.id] || []).filter(keep); }
function canCollapse(n){ return (children[n.id] || []).some(keep) && !hasCalled(n); }
var onlyCalled = false, vb = null, svgEl = null;
function layout(){
\tnodes.forEach(function(n){ n.y = undefined; });
\tvar nextY = 0;
\t(function place(n){
\t\tn.x = n.depth * COL + MARGIN;
\t\tvar kids = visKids(n);
\t\tif (!kids.length) { n.y = nextY * ROW + MARGIN; nextY += 1; }
\t\telse { kids.forEach(place); n.y = (kids[0].y + kids[kids.length-1].y) / 2; }
\t})(nodes[0]);
\treturn nextY;
}
function esc(s){ return String(s).replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;"); }
function render(){
\tlayout();
\tvar maxX = 0, maxY = 0;
\tnodes.forEach(function(n){ if (n.x !== undefined && n.y !== undefined) { maxX = Math.max(maxX, n.x); maxY = Math.max(maxY, n.y); } });
\tvar W = maxX + 320, H = maxY + 2 * MARGIN;
\tvar parts = ['<svg xmlns="http://www.w3.org/2000/svg" preserveAspectRatio="xMidYMid meet">'];
\tparts.push('<rect x="0" y="0" width="' + W + '" height="' + H + '" fill="white"/>');
\tnodes.forEach(function(n){
\t\tif (n.parent === null || n.y === undefined || !keep(n)) { return; }
\t\tvar p = nodes[n.parent];
\t\tvar dash = n.dash ? ' stroke-dasharray="4,3"' : "";
\t\tparts.push('<line x1="' + p.x + '" y1="' + p.y + '" x2="' + p.x + '" y2="' + n.y + '" stroke="' + n.line + '" stroke-width="2"' + dash + "/>");
\t\tparts.push('<line x1="' + p.x + '" y1="' + n.y + '" x2="' + n.x + '" y2="' + n.y + '" stroke="' + n.line + '" stroke-width="2"' + dash + "/>");
\t});
\tnodes.forEach(function(n){
\t\tif (n.y === undefined || !keep(n)) { return; }
\t\tvar hasKids = visKids(n).length > 0 || (n.collapsed && (children[n.id] || []).some(keep));
\t\tvar triW = 0;
\t\tif (n.collapsed && (children[n.id] || []).some(keep)) {
\t\t\tvar s = triScale(subtreeSize(n));
\t\t\ttriW = COL * 0.7 * s;
\t\t\tparts.push('<path d="M ' + n.x + ' ' + n.y + ' l ' + triW + ' ' + (-ROW * 1.5 * s) + ' l 0 ' + (ROW * 3 * s) + ' z" fill="' + n.fill + '" fill-opacity="0.35" stroke="' + n.line + '" stroke-width="1"/>');
\t\t}
\t\tif (n.dot !== false) {
\t\t\tparts.push('<circle class="node" data-id="' + n.id + '" cx="' + n.x + '" cy="' + n.y + '" r="' + Math.max(3, n.r) + '" fill="' + n.fill + '" stroke="' + (n.collapsed ? "#333" : "none") + '" stroke-width="' + (n.collapsed ? 2 : 0) + '"><title>' + esc(n.name) + (hasKids ? (hasCalled(n) ? " (called SNP path: always open)" : " (click to " + (n.collapsed ? "expand" : "collapse") + ")") : "") + "</title></circle>");
\t\t} else {
\t\t\tparts.push('<circle class="node" data-id="' + n.id + '" cx="' + n.x + '" cy="' + n.y + '" r="3" fill="#fff" stroke="#999"><title>' + esc(n.name) + "</title></circle>");
\t\t}
\t\tif (labels && n.name !== "ROOT") {
\t\t\tvar extra = n.collapsed ? " [" + subtreeSize(n) + "]" : "";
\t\t\tparts.push('<text x="' + (n.x + 6 + triW) + '" y="' + (n.y - 5) + '">' + esc(n.name) + extra + "</text>");
\t\t}
\t});
\tparts.push("</svg>");
\tdocument.getElementById("tree").innerHTML = parts.join("");
\tsvgEl = document.querySelector("#tree svg");
\tapply();
}
function subtreeSize(n){ var c = 0; (children[n.id] || []).forEach(function(k){ if (keep(k)) { c += 1 + subtreeSize(k); } }); return c; }
/*the collapsed triangle grows with the hidden subtree in steps: a 300-SNP branch should
\tlook like a clump and a 4-SNP one like a stub, without pretending exact proportion*/
function triScale(cnt){ if (cnt < 5) { return 0.6; } if (cnt < 10) { return 0.85; } if (cnt < 20) { return 1; } if (cnt < 50) { return 1.3; } if (cnt < 150) { return 1.6; } return 2; }
function bounds(){
\tvar x0 = Infinity, y0 = Infinity, x1 = -Infinity, y1 = -Infinity;
\tnodes.forEach(function(n){ if (n.y === undefined || !keep(n)) { return; }
\t\tvar ext = 0, tw = 0;
\t\tif (n.collapsed && (children[n.id] || []).some(keep)) {
\t\t\tvar s = triScale(subtreeSize(n));
\t\t\text = ROW * 1.5 * s; tw = COL * 0.7 * s;
\t\t}
\t\tx0 = Math.min(x0, n.x); y0 = Math.min(y0, n.y - ext); x1 = Math.max(x1, n.x + tw); y1 = Math.max(y1, n.y + ext); });
\treturn {x0: x0 - MARGIN, y0: y0 - MARGIN, x1: x1 + 300, y1: y1 + MARGIN};
}
function apply(){ if (svgEl && vb) { svgEl.setAttribute("viewBox", vb.x + " " + vb.y + " " + vb.w + " " + vb.h); } }
function fit(){
\tif (!svgEl) { return; }
\tvar box = bounds(), el = svgEl.getBoundingClientRect(), aspect = (el.width / el.height) || 1;
\tvar w = box.x1 - box.x0, h = box.y1 - box.y0;
\tif (w / h > aspect) { h = w / aspect; } else { w = h * aspect; }
\tvb = {x: (box.x0 + box.x1) / 2 - w / 2, y: (box.y0 + box.y1) / 2 - h / 2, w: w, h: h};
\tapply();
}
function zoomBy(f){ if (!vb) { return; } vb.x += vb.w * (1 - f) / 2; vb.y += vb.h * (1 - f) / 2; vb.w *= f; vb.h *= f; apply(); }
function leavesIfCollapsedAt(d){
\tvar c = 0;
\t(function walk(n){
\t\tvar kids = (children[n.id] || []).filter(keep);
\t\tif (!kids.length || (n.depth >= d && !hasCalled(n))) { c++; return; }
\t\tkids.forEach(walk);
\t})(nodes[0]);
\treturn c;
}
function autoCollapse(){
\tvar maxDepth = Math.max.apply(null, nodes.map(function(n){ return n.depth; }));
\tif (leavesIfCollapsedAt(maxDepth) <= OVERVIEW_LEAVES) { return; }
\tvar bestD = 0;
\tfor (var t = 1; t <= maxDepth; t++) { if (leavesIfCollapsedAt(t) <= OVERVIEW_LEAVES) { bestD = t; } else { break; } }
\tif (bestD < 2) { return; }
\tnodes.forEach(function(n){ if (n.depth >= bestD) { n.collapsed = canCollapse(n); } });
}
function relayoutFit(){ render(); fit(); }
document.getElementById("zoomIn").onclick = function(){ zoomBy(1 / 1.4); };
document.getElementById("zoomOut").onclick = function(){ zoomBy(1.4); };
document.getElementById("fit").onclick = fit;
document.getElementById("expandAll").onclick = function(){ nodes.forEach(function(n){ n.collapsed = false; }); relayoutFit(); };
document.getElementById("collapseAll").onclick = function(){ nodes.forEach(function(n){ if (n.depth >= 1) { n.collapsed = canCollapse(n); } }); relayoutFit(); };
document.getElementById("onlyCalled").onchange = function(){ onlyCalled = this.checked; relayoutFit(); };
document.getElementById("showLabels").onchange = function(){ labels = this.checked; render(); };
var labels = true;
var treeDiv = document.getElementById("tree");
var suppressClick = false;
treeDiv.addEventListener("click", function(e){
\tif (suppressClick) { return; }
\tvar t = e.target;
\tif (t.tagName !== "circle" || !t.hasAttribute("data-id")) { return; }
\tvar n = nodes[+t.getAttribute("data-id")];
\tif (n.collapsed || canCollapse(n)) { n.collapsed = !n.collapsed; relayoutFit(); }
});
treeDiv.addEventListener("dblclick", function(e){ if (e.target.tagName !== "circle") { fit(); } });
treeDiv.addEventListener("wheel", function(e){
\te.preventDefault();
\tif (!svgEl) { return; }
\tvar r = svgEl.getBoundingClientRect(), fx = (e.clientX - r.left) / r.width, fy = (e.clientY - r.top) / r.height;
\tvar f = e.deltaY > 0 ? 1.15 : 1 / 1.15;
\tvar px = vb.x + vb.w * fx, py = vb.y + vb.h * fy;
\tvb.w *= f; vb.h *= f;
\tvb.x = px - vb.w * fx; vb.y = py - vb.h * fy;
\tapply();
}, {passive: false});
var drag = null, moved = 0;
treeDiv.addEventListener("mousedown", function(e){ drag = {x: e.clientX, y: e.clientY}; moved = 0; suppressClick = false; });
window.addEventListener("mousemove", function(e){
\tif (!drag || !svgEl) { return; }
\tvar r = svgEl.getBoundingClientRect();
\tvar dx = e.clientX - drag.x, dy = e.clientY - drag.y;
\tmoved += Math.abs(dx) + Math.abs(dy);
\tvb.x -= dx * vb.w / r.width; vb.y -= dy * vb.h / r.height;
\tdrag = {x: e.clientX, y: e.clientY};
\tapply();
});
window.addEventListener("mouseup", function(){
\tif (drag && moved > 5) { suppressClick = true; }
\tdrag = null;
});
var overviewState = null;
document.getElementById("overview").onclick = function(){
\tif (!overviewState) { return; }
\tnodes.forEach(function(n){ n.collapsed = overviewState[n.id]; });
\trelayoutFit();
};
autoCollapse();
overviewState = nodes.map(function(n){ return n.collapsed; });
render();
fit();
</script>
</body>
</html>
'''

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
			'''children are rendered in name order: the children container is a set, and its
				iteration order changes between processes with the python hash seed, which made
				the newick export (and with it the drawn trees and the release gate's check on
				them) different from one run to the next'''
			children = sorted(self.children, key=lambda child: child.name)
			if len(self.children) > 0 and self.parent:
				return "({children}){name}".format(name=self.name,children=",".join(str(child) for child in children))
			elif not self.parent:  #Only the root node has this property
				return "(({children}){parent})ROOT;".format(children=",".join(str(child) for child in children),parent=self.name)
			else:
				return "{name}".format(name=self.name)
		else:
			return "{id}: {name}; children: {nchildren} ".format(id=self.id, name=self.name, nchildren = len(self.children))

	def __repr__(self):
		return "NewickNode()"

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
		self.nodeDict = {}							## Dictionary to store references to all newick nodes
		self.c_p_set = set()
		self.tree_file = "{outdir}/{name}_tree.html".format(outdir=outdir.rstrip("/"),name=name) ## output file
		self.min_required_hits = min_required_hits
		self.strictness = strictness
		## Build the newick tree
		self.newickTree = str(self.build_tree())
		## Tree colors
		self.snp_colors = {
			"derived": "#63e563", 			## Green
			"ancestral":"#984EA3" , 		## Purple
			"non_aligned": "#DDDDDD", 		## Grey
			"other_base":"#377EB8"			## Blue
		}

	def __repr__(self):
		return "NewickTree()"

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

	def new_node(self,child,nodes,parent):
		'''Function that adds a new node to the newick tree'''
		try:
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
		except KeyError:
			logger.debug("Error in adding NewickNode parent: {parent} does not exist, trying to add parent".format(parent=parent))
			parent_link = self.get_parent(parent)
			if parent_link is None:	## parent has no link in the tree table, node can not be attached
				logger.warning("Not adding node {child}: parent {parent} does not exist and has no parent link in the tree table".format(child=child,parent=parent))
				return
			t_parent,t_child,rank = parent_link
			logger.debug("Adding parent: {parent} of parent: {child}".format(parent=t_parent,child=t_child))
			self.new_node(t_child,nodes,t_parent)
			self.new_node(child,nodes,parent)
			return
		logger.debug("NewickNode p:{parent} c: {child} added".format(parent=parent,child=child))
		return

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
		for parent,child in tree:
			if parent == child:  ## root
				root = NewickNode(child, nodes[child], False)					## Create the root node
				self.nodeDict[child] = root										## Add the root node to the node dictionary
				self.nodeDict["root"] = root									## Also add this reference as "root"
				newickTree = root  												## The master parent will contain the full tree
				continue
			self.new_node(child,nodes,parent)
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
		def walk(node, parent):
			entry_id = add(node.name, parent)
			for child in sorted(node.children, key=lambda child: child.name):
				walk(child, entry_id)
		root_id = add("ROOT", None)
		walk(self.root, root_id)
		return entries

	def _confirm_path(self,entry,entries,called_snps,snplist):
		'''Confirm that the ACTUAL ancestor chain of the node, from the root of the tree down to
			the node itself, is entirely made up of called SNPs (issue #25).
			Returns (confirmed, number_of_snps_on_the_chain).
			A SNP that was not called, or that was called ancestral, breaks the chain: the walk down
			the tree should have left the branch at that point, so no node below it can be final-called.
			The quota of a complete chain is always 1, only min_required_hits gates the call now;
			competing complete paths are handled by the caller (policy of issue #47).'''
		logger.debug("Confirm path of {node}".format(node=entry["name"]))
		count = 0
		while entry["name"] != "ROOT":
			if entry["name"] in called_snps and snplist.get(entry["name"],0) != 2:
				count +=1
			else:
				'''not called (or missing from the alignments) here, or confirmed ancestral:
					the chain to this node is incomplete, the branch was left at this SNP'''
				logger.debug("Chain broken at {node}, {count} called SNPs below it".format(node=entry["name"],count=count))
				return False,count
			entry = entries[entry["parent"]]
		if count >= self.min_required_hits:
			logger.info("Confirm complete path of {count} SNPs (min required hits: {minhit})".format(count=count,minhit=self.min_required_hits))
			return True,count
		logger.debug("Complete chain too short: {count} < {minhit}".format(count=count,minhit=self.min_required_hits))
		return False,count

	def draw_tree(self,snplist,called_snps=False,save_tree=True,summary=False):
		'''Colors the SNP tree with the states found in the query alignments and writes the
			self-contained html viewer described in the module header.
		Keyword arguments:
		snplist -- a dict of the SNP names and their state
		'''
		logger.debug("Draw tree from snplist")
		entries = self._tree_entries()
		'''nodes of the tree that no alignment carried any base for: they have to be reported,
			not left out of the viewer and out of both result files (issue #27)'''
		no_data_nodes = set()
		for entry in entries:
			n = entry["name"]
			# The ancestral node is set to have a red colour, a node missing due to a gap is grey
			entry["fill"] = self.snp_colors["non_aligned"]
			entry["r"] = 5           ## ete3 "size" was a diameter: 10 drew as r=5, keep that scale
			entry["line"] = "#DDDDDD"
			entry["dash"] = True
			entry["dot"] = True
			entry["called"] = bool(called_snps) and n in set(called_snps)
			if n == "ROOT" or (n.startswith("-") and n.endswith("-")):
				'''ROOT and the "-name-" labels that get_nodes gives to nodes without any
					snp_annotation (the root of a database like that) are not SNPs: they are only
					meant to be line junctions and a label node must not be drawn as a grey SNP
					dot (issue #52, with #32 and #34)'''
				entry["dot"] = False
				continue
			if len(snplist) > 0:
				if n in snplist:
					snpvalue = snplist[n]
					if snpvalue == 2:
						## Ancestral: purple
						entry["fill"] = self.snp_colors["ancestral"]
						entry["line"] = "#000000"
						entry["dash"] = False
					elif snpvalue == 1:
						## If the SNP was derived make it green
						entry["fill"] = self.snp_colors["derived"]
						entry["r"] = 7.5
						entry["line"] = "#000000"
						entry["dash"] = False
					elif snpvalue == 3:
						## If the SNP was neither of ancestral or derived make it blue
						entry["fill"] = self.snp_colors["other_base"]
						entry["r"] = 7.5
						entry["line"] = "#000000"
						entry["dash"] = False
				else:
					'''SNP not in list: no alignment covered this position at all, keep the grey
						"non aligned" style set above. The node used to be skipped with a continue,
						which also skipped the style and left it at the ete default size, invisible
						in the pdf (issue #27)'''
					no_data_nodes.add(n)
			elif called_snps:
				### Only care about called SNPs and color all green
				if entry["called"]:
					entry["fill"] = self.snp_colors["derived"]
					entry["r"] = 7.5
				entry["line"] = "#000000"
				entry["dash"] = False
		'''the nodes without data are drawn grey now, and their database annotation is returned as
			well, so the per sample files can list them instead of hiding that they were never
			compared at all (issue #27)'''
		no_data_info = []
		if no_data_nodes:
			logger.warning("{n} SNP positions in the tree had no alignment data at all, no reference block covered them or --mask dropped them: {ids}".format(
					n=len(no_data_nodes), ids=", ".join(sorted(no_data_nodes)[:20])))
			no_data_info = self.get_no_data_info(no_data_nodes)
		## Write the tree viewer
		if save_tree:
			self.write_tree_html(entries)
			if summary:
				return
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
					root down to the node, is entirely made up of called SNPs. The old walk accepted SNPs from
					branches it should have left: _check_start picked a start node by depth adjacency in the
					flat sorted list and _confirm_path counted any called ancestor, so a node on another
					branch could be final-called (issue #25).
					When complete paths compete at the same depth no path is auto-picked: the sample stays NA
					and the candidates are reported, leaving room for the calling policy of issue #47.'''
				logger.debug("Calling SNP")
				candidates = [] ## complete called paths, all at the greatest confirmed depth
				deepest = None ## deepest called SNP found in the tree, kept for the "not confirmed" report
				for dist,tsnp in dlist: ## loop the called SNPs from deepest
					if tsnp not in node_index: ## the node name from the db is not in the tree, cannot call (issue #78)
						logger.debug("Called SNP {snp} was not found in the tree, cannot confirm it".format(snp=tsnp))
						continue
					ok,count = self._confirm_path(entries[node_index[tsnp]],entries,called_snps,snplist)
					if not deepest:
						deepest = [dist,tsnp,(False,count)]
					if not ok:
						continue
					if not candidates:
						candidates = [[dist,tsnp,count]]
					elif dist == candidates[0][0]:
						candidates.append([dist,tsnp,count])
					else: ## a shallower complete path cannot beat the deepest confirmed one
						break
				dlist = sorted(dlist,key=lambda l:l[0], reverse=False)
				if len(candidates) == 1:
					dist,tsnp,count = candidates[0]
					called = [dist,tsnp,(True,count)]
					msg = "Final SNP called: {snp}".format(snp=tsnp)
					logger.info(msg)
				elif len(candidates) > 1:
					paths = ", ".join(["{snp} ({n} SNPs)".format(snp=c[1],n=c[2]) for c in candidates])
					msg = "Two or more complete SNP paths compete at the same depth: {paths}, no path is auto-picked, check {tree}".format(paths=paths,tree=self.tree_file)
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

	def write_tree_html(self,entries):
		'''Write the self-contained viewer: the colored tree and its data are embedded in a single
			html file, no server, no javascript downloads, opens in any browser on any OS'''
		payload = json.dumps({"colors": self.snp_colors, "nodes": entries}, separators=(",", ":"))
		'''The payload sits inside a <script> block, and json.dumps leaves "<" as it is, so a node
			name from the database such as "B.1</script><script>..." closed the block and the rest
			ran as script when the tree was opened. "<" only occurs inside json strings, where
			\\u003c is the same character, so the data the viewer reads is unchanged.
			The title comes from the query file name and is escaped as html text. Both
			placeholders are filled in one pass, so a placeholder inside a name is not expanded.'''
		values = {"TITLE": html.escape(os.path.basename(self.tree_file)),
				  "PAYLOAD": payload.replace("<", "\\u003c")}
		page = re.sub(r"__(TITLE|PAYLOAD)__", lambda match: values[match.group(1)], TREE_HTML_TEMPLATE)
		with open(self.tree_file, "w") as out:
			out.write(page)

	def get_no_data_info(self,snp_ids):
		'''Look up the database annotation of the SNP nodes that no alignment carried a base for, in the
			same columns as the SNP_info rows (Name, Reference, Pos, Ancestral base, Derived base, Target
			base), so the per sample output can list the positions that were never compared instead of
			leaving them out of both result files (issue #27). The target base is "N/A": there is no base
			there at all, which is different from a base that was read and did not confirm the SNP.'''
		names = sorted(set(snp_ids))
		if not names:
			return []
		rows = []
		'''sqlite limits how many variables one statement may bind, so the lookup runs in chunks'''
		for i in range(0, len(names), 400):
			chunk = names[i:i+400]
			QUERY = '''SELECT snp_annotation.snp_id, COALESCE(snp_references.genome, snp_annotation.reference),
							snp_annotation.position, snp_annotation.ancestral_base, snp_annotation.derived_base
						FROM snp_annotation
						LEFT JOIN snp_references ON (snp_references.id = snp_annotation.genome_i)
						WHERE snp_annotation.snp_id IN ({marks})'''.format(marks=",".join(["?"]*len(chunk)))
			for snp_id,genome,position,ancestral,derived in self.database.query(QUERY,insert_val=tuple(chunk),getres=True).fetchall():
				'''the empty last column keeps the rows in the shape of SNP_FILE_COLUMNS, whose
					edge flag cannot be computed for a position that no alignment carried'''
				rows.append([str(snp_id),str(genome),str(position),str(ancestral),str(derived),"N/A",""])
		found = set([row[0] for row in rows])
		for name in names:
			if name not in found:
				'''a tree node without any annotation row: still a node the alignments never covered,
					list it rather than hide it again'''
				rows.append([name,"-","-","-","-","N/A",""])
		return sorted(rows, key=lambda row: row[0])
