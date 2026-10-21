#!/usr/bin/env python3 -c

'''
Read CanSNPer taxonomy holds a dictionary with taxonomy tree and name translation

Ported from flextaxd.modules.ReadTaxonomy(CanSNPer) to remove the flextaxd
dependency (ISSUES #22/#60), using CanSNPer2's own database functions.
'''

__version__ = "1.0"
__author__ = "David Sundell"
__credits__ = ["David Sundell"]
__license__ = "GPLv3"
__maintainer__ = "FOI bioinformatics group"
__email__ = "bioinformatics@foi.se"
__date__ = "2020-01-17"
__status__ = "Production"

from CanSNPer2.modules.DatabaseConnection import DatabaseFunctions
import logging
import gzip
import sys
logger = logging.getLogger(__name__)

class InputError(Exception):
	"""Exception raised for errors in the input."""
	def __init__(self, message):
		#self.expression = expression
		self.message = message

class ImportFormatError(Exception):
	def __init__(self, value):
		self.value = value
	def __str__(self):
		return repr(self.value)

class ReadTaxonomyBase(object):
	"""docstring for ReadTaxonomy."""
	def __init__(self, taxonomy_file=False, taxonomy_name=False, database=False,verbose=False,**kwargs):
		super(ReadTaxonomyBase, self).__init__()
		### Connect to or create database
		if database:
			self.database = DatabaseFunctions(database,verbose=verbose)
		else:
			raise InputError("No database was provided to ReadTaxonomy, abort!")
		self.taxonomy_file = taxonomy_file
		self.verbose = verbose
		self.taxonomy = {}
		self.rank = {}
		self.levelDict = {}
		self.sep = "\t"
		self.length = 0
		self.ids = 0
		self.qiime = 0
		self.skipped = []

		## Add base node
		self.root = self.add_node("root")
		self.add_rank("no rank")
		## Add basic link
		self.add_link(child=1, parent=1,rank="no rank")

	def set_qiime(self,opt):
		self.qiime = opt

	def zopen(self, path,*args, **kwargs):
		'''Redefine open to handle zipped files automatically'''
		if path.endswith(".gz"):
			if str(*args) in ["r","w"] and sys.version_info.major >= 3:
				## Python version three and above interprets binary formats as binary, add t (rt) to get text returned
				args = (str(*args)+"t",)
			else:
				args = ("rt",)
			return gzip.open(path,*args,**kwargs)
		return open(path,*args,**kwargs)

	def parse_taxonomy(self,treefile=False):
		'''Parse taxonomy information'''
		self.read_nodes(treefile=treefile)

	def set_separator(self,sep):
		'''Change separator from default (\t)'''
		self.sep = sep

	def add_rank(self, rank,qiime=False):
		'''Insert rank into database'''
		try:
			return self.rank[rank]
		except KeyError:
			pass
		if rank:
			if qiime:
				rank_i = self.database.add_rank(self.levelDict.get(rank))
			else:
				rank_i = self.database.add_rank(rank)
			self.rank[rank] = rank_i
		else:
			return False
		return rank_i

	def add_link(self, child=None, parent=None,rank="no rank"):
		'''Add relationship in tree'''
		if not child and parent:
			raise InputError("A link requires both a child and a parent!")
		self.database.add_link(child=child,parent=parent,rank=self.rank[rank])
		self.ids+=1

	def add_node(self, description,id=False):
		'''Add node to tree
			extend databaseFunction add_node function using self.taxonomy
			dict to keep track of which nodes were already added
		'''
		logger.debug("Node desc: {desc} id: {id}".format(desc=description,id=id))
		self.taxid_base = self.database.add_node(description,id)
		self.taxonomy[description] = self.taxid_base
		return self.taxid_base

	def read_nodes(self, treefile=False):
		'''Read a tab separated node file and store a dictionary with names in object'''
		logger.info("Read nodes in taxonomy file {}".format(treefile))
		swap = False
		rank = "no rank"  #Base rank if rank is not used
		if not treefile:
			treefile = self.taxonomy_file
		with self.zopen(treefile, "r") as _treefile:
			headers = _treefile.readline().strip().split(self.sep)
			if "parent" not in headers or "child" not in headers:
				logger.debug("Headers:  {h} separator [{sep}]".format(h=headers,sep=self.sep))
				raise InputError("Your input tree file does not contain the headers to specify child and parent!")
			if headers[0] == "parent":
				swap = True
			logger.debug("Swap: {swap}".format(swap=swap))
			for tree_row in _treefile:
				data = tree_row.strip().split(self.sep)
				logger.debug(data)
				if len(data) > 2:
					rank = data.pop().strip()
					if rank.strip() != "":
						self.add_rank(rank)
				if swap:
					data[0],data[1] = data[1],data[0]
				'''Add nodes into database'''
				if data[0] == "":
					'''Check for empty rows'''
					continue
				for node in data:
					logger.debug(node)
					node = node.strip()
					try:
						self.taxonomy[node]
					except KeyError:
						logger.debug("Add node: {node}".format(node=node))
						self.taxonomy[node] = self.add_node(node)
						self.ids +=1
				self.add_link(child=self.taxonomy[data[0].strip()],parent=self.taxonomy[data[1].strip()],rank=rank)
				self.length +=1
			self.database.commit()


class ReadTaxonomyCanSNPer(ReadTaxonomyBase):
	"""docstring for ReadTaxonomyCanSNPer."""
	def __init__(self, taxonomy_file=False, database=".canSNPdb",  taxid_base=1,root_name=False,rank="family", cnp=False, verbose=False,**kwargs):
		super(ReadTaxonomyCanSNPer, self).__init__(taxonomy_file=taxonomy_file, database=database,verbose=verbose,**kwargs)
		self.input = taxonomy_file
		self.taxonomy = {}
		self.taxid_num = taxid_base
		## Initiate database
		logger.info(taxonomy_file)
		logger.debug(root_name)
		if not root_name:
			logger.info("Fetching root name from file")
			root_name = self.get_root_name(taxonomy_file)
		if not cnp:
			logger.info("Adding, cellular organism node")
			c_i = self.add_node("cellular organisms")
			logger.info("Adding root node {node}!".format(node=root_name))
		root_i = self.add_node(root_name)
		self.taxid_num =taxid_base ## reset
		logger.debug("Adding ranks!")
		if rank != "no rank":
			self.add_rank("no rank")
		self.add_rank(rank)
		if root_name != "cellular organisms" and not cnp:
			self.add_link(c_i,1)
		if root_name != "root" and not cnp:
			self.add_link(root_i,c_i,rank=rank) ## root is always added as top parent if not one force 1
			self.taxonomy[root_i] = 1  ## Add local link to root
		self.names = {}
		self.root = root_i
		self.length = 0
		self.ids = 0

		self.database.commit()
		logger.debug("Root: {root} link: [{base}]".format(root=self.root, base=self.taxid_num))

	def get_root_name(self,fin):
		'''Get the root name of the CanSNP file unless root is specified'''
		with self.zopen(fin) as f:
			firstline = f.readline()
			firstline = firstline.strip().replace("\t",";")
			root = firstline.split(";") ## Root should be either the leftmost annotation in the tree or the only annotation on that row
			if isinstance(root,list):
				root = root[0]
		return root

	def add_SNP(self,nodes,i):
		'''Add name of node to database'''
		name = nodes[i].strip()
		logger.debug("Parent did not exit Add parent: {name}".format(name=name))

		#i-=1 ## Check next parent
		try:
			parent_i = self.taxonomy[nodes[i-1]]  ## check if current nodes parent exists
		except KeyError:                          ## Parent did not exist, keep walking up the tree and add all parents until root!
			if i < -len(nodes):
				return self.root                ## The root has been reached return index of root
			else:                                ## Parent didn't exist again, add parent to this node and then add the link to that parent
				parent_i = self.add_SNP(nodes,i-1)## Add node of parent
		node_i = self.add_node(name)  ## Add node
		if node_i:
			self.taxid_num += 1
		self.add_link(node_i, parent_i)            ## Add link to next parent
		return node_i     ##  index of child

	def parse_taxonomy(self):
		'''Retrieve node description from CanSNPer formatted tree'''
		logger.info("Parse CanSNP tree file")
		with self.zopen(self.taxonomy_file,"r") as f:
			for row in f:
				row = row.strip().replace("\t",";")  ## Also accept tab separated tree files
				nodes = row.strip().split(";")  ## get node and all its parents in a list
				logger.debug(nodes)
				child = nodes[-1]  ## get name of child node
				'''If the tree was not properly formatted an a parent is missing make sure that function works anyway by adding any parent node above child'''
				try:
					pname = nodes[-2].strip()
					parent_i = self.taxonomy[pname]  ## Check if parent of child node exists
				except KeyError: ## parent node does not exist, add parent
					parent_i = self.add_SNP(nodes,-2)
				except IndexError: ## Should be first row with only one node (parent)
					self.taxonomy[child] = self.root
					parent_i = self.root
					continue
				'''Now all parents exists, add new child node and add the link'''
				child_i = self.add_node(child)
				self.database.add_link(child_i,parent_i)
		self.database.commit()
		self.length = self.taxid_num - self.root                ## Check number of new nodes added
		logger.info("New taxonomy ids assigned {taxidnr}".format(taxidnr=self.length))
