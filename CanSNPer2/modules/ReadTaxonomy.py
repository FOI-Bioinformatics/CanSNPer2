'''
Read NCBI taxonomy dmp files (nodes or names) and holds a dictionary
'''

from CanSNPer2.modules.DatabaseConnection import CanSNPdbFunctions
import logging
logger = logging.getLogger(__name__)

class InputError(Exception):
	"""Exception raised for errors in the input."""
	def __init__(self, message):
		#self.expression = expression
		self.message = message

class ReadTaxonomy(object):
	"""docstring for ReadTaxonomy."""
	def __init__(self, database=".taxonomydb",verbose=False):
		super(ReadTaxonomy, self).__init__()
		### Connect to or create database
		logger.info("Read Taxonomy init")
		self.database = CanSNPdbFunctions(database,verbose=verbose)
		self.verbose = verbose
		self.taxonomy = {}
		self.rank = {}
		self.levelDict = {}
		self.names = {}
		self.sep = "\t"
		self.length = 0
		self.ids = 0

	def set_separator(self,sep):
		'''Change separator from default (\t)'''
		self.sep = sep

	def add_rank(self, rank,ncbi=False):
		'''Insert rank into database'''
		if not ncbi:
			rank_i = self.database.add_rank(self.levelDict.get(rank))
		else:
			rank_i = self.database.add_rank(rank)
		self.rank[rank] = rank_i
		return rank_i

	def add_link(self, child, parent,rank="no rank"):
		'''Add relationship in tree'''
		#print(child,parent,rank)
		try:
			rank = self.rank[rank]
		except KeyError:
			rank = False
		self.database.add_link(child,parent,rank)
		self.ids+=1

	def add_node(self, description,id=False,table="Node"):
		'''Add node to tree
			extend databaseFunction add_node function using self.taxonomy
			dict to keep track of which nodes were already added
		'''
		self.taxid_base = self.database.add_node(description,id,table=table)
		self.taxonomy[description] = self.taxid_base
		return self.taxid_base

	def parse_taxonomy(self,treefile=False):
		'''Parse taxonomy information'''
		self.read_nodes(treefile=treefile)

	def read_nodes(self, treefile=False):
		'''Read a tab separated node file and store a dictionary with names in object'''
		if self.verbose: print("Read nodes in taxonomy file {}".format(treefile))
		swap = False
		if not treefile:
			treefile = self.taxonomy_file
		with open(treefile, "r") as _treefile:
			headers = _treefile.readline().strip().split(self.sep)
			if "parent" not in headers or "child" not in headers:
				raise InputError("Your input tree file does not contain the headers to specify child and parent!")
			if headers[0] == "child":
				swap = True
			for tree_row in _treefile:
				data = tree_row.strip().split(self.sep)
				if swap:
					data[0],data[1] = data[1],data[0]
				'''Add nodes into database'''
				if data[0] == "":
					'''Check for empty rows'''
					continue
				for node in data:
					try:
						self.names[node]
					except KeyError:
						self.names[node] = self.database.add_node(node)
						self.ids +=1
				self.database.add_link(child=data[0],parent=data[1])
				self.length +=1
			self.database.commit()
