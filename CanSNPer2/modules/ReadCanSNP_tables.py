'''
Read CanSNPer1 tables takes CanSNPer version 1 input files and parses them to a database
'''

__author__ = "David Sundell"
__credits__ = ["David Sundell"]
__license__ = "GPLv3"
__maintainer__ = "FOI bioinformatics group"
__email__ = "bioinformatics@foi.se"
__date__ = "2019-07-08"
__status__ = "Production"

from CanSNPer2.modules.ReadTaxonomy import ReadTaxonomy

class ReadCanSNPer(ReadTaxonomy):
	"""docstring for ReadCanSNPer."""
	def __init__(self, database=".canSNPdb",  taxid_base=1,root_name="Francisella", version=2, verbose=False):
		super(ReadCanSNPer, self).__init__(database=database,verbose=verbose)
		self.annotation = {}
		self.taxid_base = taxid_base

		'''Annotations'''
		self.genomes = {}		## path to genomes

		self.names = {}
		self.root = 1
		self.length = 0
		self.ids = 0

		self.version = version

	def add_SNP(self,nodes,i):
		'''Add name of node to database'''
		name = nodes[i].strip()
		#i-=1 ## Check next parent
		try:
			parent_i = self.annotation[nodes[i-1]]  ## check if current nodes parent exists
		except KeyError:  						## Parent did not exist, keep walking up the tree and add all parents until root!
			if i < -len(nodes):
				return self.root				## The root has been reached return index of root
			else:								## Parent didn't exist again, add parent to this node and then add the link to that parent
				parent_i = self.add_SNP(nodes,i-1)## Add node of parent
				#self.annotation[] = parent_i
		node_i = self.add_node_annotation([name],["snp_id"],table="nodes")  ## Add node
		self.add_link(node_i, parent_i)			## Add link to next parent
		return node_i 	##  index of child

	def add_node_annotation(self,data,headers = [], table="nodes"):
		'''Add full annotation of a CanSNP'''
		#print(self.taxid_base, data)
		info = {}
		for i in range(len(headers)):
			info[headers[i]] = data[i]  ## the
		try:
			taxid_base = self.annotation[info[headers[0]]]
			update_info = { "set_value": taxid_base,
					"where_column": "node_i",
					"set_column": info
				}
			self.database.update(update_info,table)
		except KeyError:
			taxid_base = self.database.insert(info, table=table)
			self.annotation[info[headers[0]]] = taxid_base
			self.taxid_base+=1
		return taxid_base

	def add_annotation(self,data,headers = [], table="nodes"):
		'''Add full annotation of a CanSNP'''
		info = {}
		for i in range(len(headers)):
			info[headers[i]] = data[i]  ## the
		annotation_i = self.database.insert(info, table=table)
		return annotation_i

	def parse_genomes(self,genomes_file):
		'''Parse genomes file and add information about genome references in the database'''
		if self.verbose: print("Parse genome files")
		self.genomeDict = {}
		with open(genomes_file, "r") as f:
			headers = f.readline().strip().split("\t")
			for row in f:
				genomeinfo = row.strip().split("\t")
				genome_i = self.add_annotation(genomeinfo, headers, table="genomes")
		self.database.commit()
		res = self.database.get_all(database=self.database, table="genomes",columns="genome_i,genome_id")
		for genome_i,genome_id in res:
			self.genomeDict[genome_id] = genome_i
		return

	def parse_snp_annotations(self,annotation_file):
		'''This functions parses the canSNP node table and adds defined SNPs to the database
			#SNP	Reference	Genome	Genome_Position	Derived_base	Ancestral_base
			translate headers into -> snp_id organism        reference       genome  position        derived_base    ancestral_base
		'''
		if self.verbose: print("Parse CanSNP node file")
		#print(self.annotation)
		with open(annotation_file,"r") as f:
			headers = f.readline().strip("#").split(self.sep)
			if self.version == 1:
				headers = "snp_id,organism,reference,genome,position,derived_base,ancestral_base".split(",")
				del headers[1]
			for row in f:
				data = row.strip().split(self.sep)  ## split up all the information in a node
				if self.version == 1:
					'''Remove "organism" column as databases will be separated on organim'''
					del data[1]
				headers[2] = "genome_i"
				data[2] = self.genomeDict[data[2]]  ## Get genome index third column must be genome
				self.add_node_annotation(data,headers)
		self.database.commit()  ## Commit changes to database
		return

	def parse_taxonomy(self,taxonomy_file):
		'''Retrieve node description from CanSNPer formatted tree'''
		if self.verbose: print("Parse CanSNP tree file")
		with open(taxonomy_file,"r") as f:
			for row in f:
				nodes = row.strip().split(";")  ## get node and all its parents in a list
				child = nodes[-1]  ## get name of child node
				if child == "":
					continue
				'''If the tree was not properly formatted an a parent is missing make sure that function works anyway by adding any parent node above child'''
				try:
					parent_i = self.annotation[nodes[-2].strip()]  ## Check if parent of child node exists
				except KeyError: ## parent node does not exist, add parent
					#print(self.taxid_base)
					parent_i = self.add_SNP(nodes,-2)
				except IndexError: ## Should be first row with only one node (parent)
					child_i = self.add_node_annotation([child],["snp_id"],table="nodes")
					self.annotation[child] = self.root
					parent_i = self.root
					self.database.add_link(parent_i,parent_i)

					#self.taxid_base = self.root+1
					continue
				'''Now all parents exists, add new child node and add the link'''
				try:
					child_i = self.annotation[child]
				except KeyError:
					child_i = self.add_node_annotation([child],["snp_id"],table="nodes")
				self.database.add_link(child_i,parent_i)
		self.database.commit()									## Commit changes to database
		self.length = self.taxid_base - self.root				## Check number of new nodes added
		print("New taxonomy ids assigned {taxidnr}".format(taxidnr=self.length))
