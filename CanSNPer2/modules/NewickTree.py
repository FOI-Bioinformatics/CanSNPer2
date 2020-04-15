#!/usr/bin/env python3 -c

'''
Module to read and write newick trees

'''
from flextaxd.modules.database.DatabaseConnection import DatabaseFunctions
from CanSNPer2.modules.DatabaseConnection import CanSNPdbFunctions
from ete3 import Tree, faces, AttrFace, TreeStyle, NodeStyle

import sys
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
			if len(self.children) > 0 and self.parent:
				return "({children}){name}".format(name=self.name,children=",".join(str(child) for child in self.children))
			elif not self.parent:  #Only the root node has this property
				return "(({children}){parent})ROOT;".format(children=",".join(str(child) for child in self.children),parent=self.name)
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

	def __init__(self, database,name="newick",outdir="./",min_required_hits=3):
		super(NewickTree, self).__init__()
		self.database = CanSNPdbFunctions(database) ## Initiate database connection with CanSNPdbFunctions
		self.nodeDict = {}							## Dictionary to store references to all newick nodes
		self.c_p_set = set()
		self.tree_file = "{outdir}/{name}_tree.pdf".format(outdir=outdir.rstrip("/"),name=name) ## output file
		self.min_required_hits = min_required_hits
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
		QUERY = '''SELECT parent,child,rank_i FROM tree WHERE child = "{node}"'''.format(node=name)
		res = self.database.query(QUERY).fetchone()
		return res

	def get_child(self,name):
		'''return child'''
		#QUERY = '''SELECT child,child,rank FROM tree LEFT JOIN rank on (tree.rank_i = rank.rank_i) WHERE child = "{node}"'''.format(node=name)
		QUERY = '''SELECT parent,child,rank_i FROM tree WHERE parent = "{node}"'''.format(node=name)
		res = self.database.query(QUERY).fetchone()
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
			t_parent,t_child,rank = self.get_parent(parent)
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
		logger.info("Nodes: {n} Links: {l}".format(n=len(nodes),l=len(tree)))
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

	def _confirm_path(self,dist_list,called_snps,snplist):
		'''Confirm path of snps'''
		try:
			logger.debug("Confirm path")
			count = 0
			logger.info(dist_list)
			dist,node = dist_list
			while node:
				if node.name == "ROOT":
					pass
				elif node.name in called_snps and snplist[node.name] != 2:
					logger.debug("count1")
					count +=1
				else: ## If ancestral node is found in the path it means it is confirmed ancestral, therefore the current path is not correct!
					count = 0
				node = node.up
			logger.debug(count)
			quota = float(count)/dist
			logger.info("Confirm quota: {quota}".format(quota=quota))
			if quota > 0.7 and count >=self.min_required_hits:
				return True,count
			return False,count
		except KeyError as e:
			logger.error("KeyError, {e}".format(e=e))
			logger.critical("Could not confirm path!")

	def _check_start(self,dlist):
		'''Check so that the most distant node actually have support for at least three parents, else check the next deepest node'''
		n=0
		f_dist,f_node = dlist[n]
		_dist,_node = dlist[n]
		sublist = dlist[n:] ## remove first element is that first node is child of the first object in the list
		n+=1
		lcount = 0
		logger.debug("Check which start node is the deepest with at least three available parents")
		logger.debug(dlist)
		logger.debug(sublist)
		for dist,node in sublist:
			logger.debug("1-{dist}: {node}".format(dist=dist,node=node))
			logger.debug("2-{dist}: {node}".format(dist=_dist,node=_node))
			if (_dist != dist+1) and (_dist != dist):
				logger.debug("Leaf does not have parent support go to next leaf!")
				lcount = 0
				n+=1
				f_dist,f_node = dist,node ## Final node
				#sublist = dlist[n:] ## remove first n elements
			elif lcount >= self.min_required_hits:
				logger.debug("Final count ok return")
				if n >= 1:
					sublist = dlist[n-1:] ## remove elements -1
					logger.debug(sublist)
				else:
					sublist = dlist
				return f_dist,f_node,sublist
			else:
				logger.debug("Parent ok continue")
				lcount+=1
			_dist,_node = dist,node
			if lcount >= self.min_required_hits and len(dlist) <= self.min_required_hits:
				return f_dist,f_node,dlist
		return False,False,dlist

	def draw_ete3_tree(self,snplist,called_snps=False,save_tree=True):
		'''Draws a phylogenetic tree using ETE3
		Keyword arguments:
		snplist -- a list of the SNP names, positions and state
		'''
		logger.debug("Draw tree from snplist")
		try:
			tree = Tree(self.newickTree, format=1)
		except:
			logger.debug(self.newickTree)
		farthest_leaf, tree_depth = tree.get_farthest_leaf()
		for n in tree.traverse():
			# The ancestral node is set to have a red colour
			nstyle = NodeStyle()
			nstyle["fgcolor"] = self.snp_colors["ancestral"]
			nstyle["size"] = 10
			nstyle["vt_line_color"] = "#000000"
			nstyle["hz_line_color"] = "#000000"
			nstyle["vt_line_type"] = 0
			nstyle["hz_line_type"] = 0
			nstyle["vt_line_width"] = 2
			nstyle["hz_line_width"] = 2
			for snp in snplist.keys():
				if n.name == snp and snplist[snp] == 0:
					# If the SNP is missing due to a gap, make it grey
					nstyle["fgcolor"] = self.snp_colors["non_aligned"]
					nstyle["size"] = 10
					nstyle["vt_line_color"] = "#DDDDDD"
					nstyle["hz_line_color"] = "#DDDDDD"
					nstyle["vt_line_type"] = 1
					nstyle["hz_line_type"] = 1
				elif n.name == snp and snplist[snp] == 1:
					## If the SNP was derived make it green
					nstyle["fgcolor"] = self.snp_colors["derived"]
					nstyle["size"] = 15
					nstyle["vt_line_color"] = "#000000"
					nstyle["hz_line_color"] = "#000000"
					nstyle["vt_line_type"] = 0
					nstyle["hz_line_type"] = 0
				elif n.name == snp and snplist[snp] == 3:
					## If the SNP was neither of ancestral or derived make it blue
					nstyle["fgcolor"] = self.snp_colors["other_base"]
					nstyle["size"] = 15
					nstyle["vt_line_color"] = "#000000"
					nstyle["hz_line_color"] = "#000000"
					nstyle["vt_line_type"] = 0
					nstyle["hz_line_type"] = 0
			if n.name != "ROOT": ## Root should be just a line not a false "ancenstral node"
				n.set_style(nstyle)
		ts = TreeStyle()
		ts.show_leaf_name = False  						# Do not print(leaf names, they are added in layout)
		ts.show_scale = False  							# Do not show the scale
		ts.layout_fn = self.CanSNPer_tree_layout  		# Use the custom layout
		ts.optimal_scale_level = 'full' 			 	# Fully expand the branches of the tree
		scale_factor = 500								# Tree scale factor, increases depth of tree to allow a higher resolution tree

		## Render tree and save to pdf
		if save_tree:
			tree.render(self.tree_file, tree_style=ts, w=tree_depth * scale_factor)
		dlist = []
		confirmed = [False,0]
		final = "NA"
		logger.info("Called snps: {called}".format(called=called_snps))
		if called_snps:
			logger.info("Call SNPs!")
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
				logger.info(dlist)
				logger.info(len(dlist))
				logger.info(self.min_required_hits)
				if len(dlist) < self.min_required_hits: ## the number of SNPs found is too small to be valid,
					msg = "The number of SNPs found is too small (min {min}), Cannot call SNPs.".format(min=self.min_required_hits)
					logger.debug(msg)
					return False,msg
				## Check so that the starting node have at least tree consecutive nodes, otw look for another start
				dist,node,dlist = self._check_start(dlist)
				if not node:
					msg = "No valid start SNP found."
					logger.info(msg)
					return False,msg
				else:#dist,node = dlist[0]
					logger.debug("Calling SNP")
					node = tree.search_nodes(name=node)[0]
					### Loop different nodes
					confirmed = self._confirm_path([dist,node],called_snps,snplist)
					dlist[0].append(confirmed)
					if confirmed[0]:
						final = node.name.strip()
						msg = "Final SNP called: {snp}".format(snp=final)
						logger.info(msg)
					else:
						msg = "SNP could not be confirmed!"
				logger.debug(dlist)
			except KeyError as e:
				logger.error("KeyError, {e}".format(e=e))
				msg = "called_snps function failed, KeyError {e}".format(e=e)
				logger.debug(msg)
				return False,msg
		else:
			msg = "No SNPs were called."
			logger.debug(msg)
			return False,msg
		return dlist[0],msg
