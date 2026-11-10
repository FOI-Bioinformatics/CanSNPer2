'''
Module to read and write newick trees

'''
from CanSNPer2.modules.DatabaseConnection import CanSNPdbFunctions
import os
import sys
'''ete3 renders trees through Qt, which aborts the whole process when it cannot open a display
	(headless servers, containers), so select the offscreen platform unless the user chose one'''
if sys.platform.startswith("linux") and not os.environ.get("DISPLAY") and not os.environ.get("QT_QPA_PLATFORM"):
    os.environ["QT_QPA_PLATFORM"] = "offscreen"
import ete3
'''ete3 only exposes the drawing classes (NodeStyle, TreeStyle, AttrFace, faces) through its treeview
	subpackage, which fails to import when no Qt binding is installed, and ete3 swallows that error.
	An ete3 from the etetoolkit conda channel without pyqt is exactly that case: ete3 imported fine but
	ete3.NodeStyle() died with "module 'ete3' has no attribute 'NodeStyle'" in the middle of drawing
	(issue #19). Import the names explicitly from where they live, so the failure is one ImportError at
	import time, and raise a readable error when a tree is actually drawn (issue #122).
'''
try:
	'''the names live in different parts of the treeview subpackage: NodeStyle and TreeStyle in main,
		the face classes in faces'''
	from ete3.treeview.main import NodeStyle, TreeStyle
	from ete3.treeview.faces import AttrFace
	from ete3.treeview import faces
	ETENEWICK_OK = True
	ETE_TREEVIEW_ERROR = ""
except ImportError as e:
	NodeStyle = TreeStyle = AttrFace = faces = None
	ETENEWICK_OK = False
	'''the message ete3's own import produced, kept so the error at draw time names what is missing'''
	ETE_TREEVIEW_ERROR = str(e)

import logging
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
		self.tree_file = "{outdir}/{name}_tree.pdf".format(outdir=outdir.rstrip("/"),name=name) ## output file
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
		logger.debug("Tree complete, return newickTree")
		#logger.debug(newickTree)
		return newickTree

	def CanSNPer_tree_layout(self,node):
		'''Layout function for ETE3 trees.'''
		# Adds the name face to the image at the top side of the branch
		if not node.is_root():
			faces.add_face_to_node(AttrFace("name"), node, column=0, position="branch-top")

	def _confirm_path(self,node,called_snps,snplist):
		'''Confirm that the ACTUAL ancestor chain of the node, from the root of the tree down to
			the node itself, is entirely made up of called SNPs (issue #25).
			Returns (confirmed, number_of_snps_on_the_chain).
			A SNP that was not called, or that was called ancestral, breaks the chain: the walk down
			the tree should have left the branch at that point, so no node below it can be final-called.
			The quota of a complete chain is always 1, only min_required_hits gates the call now;
			competing complete paths are handled by the caller (policy of issue #47).'''
		logger.debug("Confirm path of {node}".format(node=node.name))
		count = 0
		while node.name != "ROOT":
			if node.name in called_snps and snplist.get(node.name,0) != 2:
				count +=1
			else:
				'''not called (or missing from the alignments) here, or confirmed ancestral:
					the chain to this node is incomplete, the branch was left at this SNP'''
				logger.debug("Chain broken at {node}, {count} called SNPs below it".format(node=node.name,count=count))
				return False,count
			node = node.up
		if count >= self.min_required_hits:
			logger.info("Confirm complete path of {count} SNPs (min required hits: {minhit})".format(count=count,minhit=self.min_required_hits))
			return True,count
		logger.debug("Complete chain too short: {count} < {minhit}".format(count=count,minhit=self.min_required_hits))
		return False,count

	def draw_ete3_tree(self,snplist,called_snps=False,save_tree=True,summary=False):
		'''Draws a phylogenetic tree using ETE3
		Keyword arguments:
		snplist -- a list of the SNP names, positions and state
		'''
		logger.debug("Draw tree from snplist")
		if not ETENEWICK_OK:
			raise ImportError("ete3 has no drawing classes (NodeStyle, TreeStyle, AttrFace), no tree can be drawn: {why}. Install a Qt binding for ete3: pip install PyQt5, or in conda conda install pyqt - an ete3 from the etetoolkit channel needs pyqt on its own (issues #19, #122)".format(why=ETE_TREEVIEW_ERROR))
		try:
			tree = ete3.Tree(self.newickTree, format=1)
		except Exception as e:
			logger.debug(self.newickTree)
			raise ValueError("The newick tree in the database could not be parsed by ete3, cannot draw the tree: {err}".format(err=e))
		farthest_leaf, tree_depth = tree.get_farthest_leaf()
		'''nodes of the tree that no alignment carried any base for: they have to be reported,
			not left out of the pdf and out of both result files (issue #27)'''
		no_data_nodes = set()
		for n in tree.traverse():
			# The ancestral node is set to have a red colour
			nstyle = NodeStyle()
			# If the SNP is missing due to a gap, make it grey
			nstyle["fgcolor"] = self.snp_colors["non_aligned"]
			nstyle["size"] = 10
			nstyle["vt_line_color"] = "#DDDDDD"
			nstyle["hz_line_color"] = "#DDDDDD"
			nstyle["vt_line_type"] = 1
			nstyle["hz_line_type"] = 1
			nstyle["vt_line_width"] = 2
			nstyle["hz_line_width"] = 2
			if len(snplist) > 0:
				try:
					snpvalue = snplist[n.name]
				except KeyError:
					if n.name == "ROOT" or (n.name.startswith("-") and n.name.endswith("-")):
						'''ROOT and the "-name-" labels that get_nodes gives to nodes without any
							snp_annotation (the root of a database like that) are not SNPs: they keep
							no style, the root is only meant to be a line and a label node must not be
							drawn as a grey SNP dot (issue #52, with #32 and #34)'''
						continue
					'''SNP not in list: no alignment covered this position at all, keep the grey
						"non aligned" style that was set above. The node used to be skipped with a
						continue, which also skipped set_style below and left the node at the ete3
						default size, so it was invisible in the pdf (issue #27)'''
					no_data_nodes.add(n.name)
				else:
					if snpvalue == 2:
						nstyle["fgcolor"] = self.snp_colors["ancestral"]
						nstyle["size"] = 10
						nstyle["vt_line_color"] = "#000000"
						nstyle["hz_line_color"] = "#000000"
						nstyle["vt_line_type"] = 0
						nstyle["hz_line_type"] = 0
						nstyle["vt_line_width"] = 2
						nstyle["hz_line_width"] = 2
					elif snpvalue == 1:
						## If the SNP was derived make it green
						nstyle["fgcolor"] = self.snp_colors["derived"]
						nstyle["size"] = 15
						nstyle["vt_line_color"] = "#000000"
						nstyle["hz_line_color"] = "#000000"
						nstyle["vt_line_type"] = 0
						nstyle["hz_line_type"] = 0
					elif snpvalue == 3:
						## If the SNP was neither of ancestral or derived make it blue
						nstyle["fgcolor"] = self.snp_colors["other_base"]
						nstyle["size"] = 15
						nstyle["vt_line_color"] = "#000000"
						nstyle["hz_line_color"] = "#000000"
						nstyle["vt_line_type"] = 0
						nstyle["hz_line_type"] = 0

			elif called_snps:
				if set([n.name]) & set(called_snps):
					### Only care about called SNPs and color all green
					nstyle["fgcolor"] = self.snp_colors["derived"]
					nstyle["size"] = 15
				nstyle["vt_line_color"] = "#000000"
				nstyle["hz_line_color"] = "#000000"
				nstyle["vt_line_type"] = 0
				nstyle["hz_line_type"] = 0

			if n.name != "ROOT": ## Root should be just a line not a false "ancenstral node"
				n.set_style(nstyle)
		'''the nodes without data are drawn grey now, and their database annotation is returned as
			well, so the per sample files can list them instead of hiding that they were never
			compared at all (issue #27)'''
		no_data_info = []
		if no_data_nodes:
			logger.warning("{n} SNP positions in the tree had no alignment data at all, no reference block covered them or --mask dropped them: {ids}".format(
					n=len(no_data_nodes), ids=", ".join(sorted(no_data_nodes)[:20])))
			no_data_info = self.get_no_data_info(no_data_nodes)
		ts = TreeStyle()
		ts.show_leaf_name = False  						# Do not print(leaf names, they are added in layout)
		ts.show_scale = False  							# Do not show the scale
		ts.layout_fn = self.CanSNPer_tree_layout  		# Use the custom layout
		ts.optimal_scale_level = 'full' 			 	# Fully expand the branches of the tree
		scale_factor = 500								# Tree scale factor, increases depth of tree to allow a higher resolution tree

		## Render tree and save to pdf
		if save_tree:
			tree.render(self.tree_file, tree_style=ts, w=tree_depth * scale_factor)
			if summary:
				return
		dlist = []
		logger.info("Called snps: {called}".format(called=called_snps))
		if called_snps:
			try:
				troot = tree.get_tree_root()
				'''Calculate the distance from all nodes to the root'''
				for tsnp in called_snps:
					try:
						dist = tree.get_distance(tsnp,troot)
						dlist.append([dist,tsnp])
					except ValueError:
						if tsnp == "T/N.1":
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
					matches = tree.search_nodes(name=tsnp)
					if not matches: ## the node name from the db is not in the ete3 tree, cannot call (issue #78)
						logger.debug("Called SNP {snp} was not found in the tree, cannot confirm it".format(snp=tsnp))
						continue
					ok,count = self._confirm_path(matches[0],called_snps,snplist)
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
					msg = "Two or more complete SNP paths compete at the same depth: {paths}, no path is auto-picked, check {pdf}".format(paths=paths,pdf=self.tree_file)
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
