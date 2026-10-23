'''
	Modify and update the CanSNPer2 database
'''

from CanSNPer2.modules.ModifyTree import ModifyTree, InputError
from CanSNPer2.modules.DatabaseConnection import NameNotFoundError, check_source_headers, GENOME_SOURCE_COLUMNS, ANNOTATION_SOURCE_COLUMNS
from datetime import date
import logging
logger = logging.getLogger(__name__)

class ModifyCanSNPer2Database(ModifyTree):
	"""ModifyCanSNPer2Database adds a few important functions only nessesary when modifying a database"""
	def __init__(self, database, mod_database=False, mod_file=False, separator="\t",verbose=False,parent=False,replace=False,remove=False,**kwargs):
		logger.info("Load ModifyFunctions")
		super().__init__(database, mod_database=mod_database, mod_file=mod_file, separator=separator,verbose=verbose,parent=parent,replace=replace,remove=remove)
		self.genomes = {}
		self.today = date.today()

	def get_nodes(self, database=False,col=False):
		'''Retrieve the whole node info table of the database to decrease the number of database calls!'''
		return self.taxonomydb.get_nodes()

	def get_genomes(self, table):
		return self.taxonomydb.get_genomes(table="snp_references")

	def add_annotation(self, data, id=False,genome_i=False):
		'''Add snp annotation to annotation table'''
		logger.debug("Add annotation")

		try:
			data["date"]
		except KeyError:
			data["date"] = self.today.strftime("%d/%m/%y")
		info = {
				"snp_id":            data["snp_id"],
				"position":         data["position"],
				"ancestral_base":     data["ancestral_base"],
				"derived_base":        data["derived_base"],
				"reference":         data["reference"],
				"date":             data["date"],
				"genome_i":         genome_i
		}
		### If ID is supplied skip autoincrement and add specific ID
		if id:
			info["node_id"] = id
		logger.debug(info)
		taxid_base = self.taxonomydb.insert(info, table="snp_annotation")
		self.taxonomydb.commit()
		return taxid_base

	def add_genome(self,data):
		'''Add genome annotation'''
		logger.debug("Add genome")
		info = {
				"genome":            data["genome"],
				"strain":         data["strain"],
				"genbank_id":     data["genbank_id"],
				"refseq_id":        data["refseq_id"],
				"assembly_name":         data["assembly_name"],
				#"comment":             data["comment"]
		}
		logger.debug(info)
		return self.taxonomydb.insert(info, table="snp_references")

	def load_annotation_file(self,annotation_file):
		'''Read the annotation source file'''
		nodes = self.get_nodes()
		if len(self.genomes) == 0:
			self.genomes = self.get_genomes(table="snp_references")
		with open(annotation_file) as af:
			headers = af.readline().lstrip("#").strip().split("\t")
			logging.debug(headers)
			check_source_headers(headers,ANNOTATION_SOURCE_COLUMNS,annotation_file)
			dataArr = []
			for row in af:
				data = row.strip().split("\t")
				if len(data) < len(headers): ## short row would run out of values in the header zip below (issue #88)
					logging.warning("Skipped short row in annotation file: '{row}'".format(row=row.strip()))
					continue
				ddict = {}
				for i in range(len(headers)):
					head,val = headers[i],data[i]
					ddict[head] = val
				dataArr.append(ddict)
				try:
					_id = nodes[ddict["snp_id"]]
				except KeyError:
					logging.warning("WARNING: {id} was not found in the tree, check your source files!".format(id=ddict["snp_id"]))
					_id = False
				try:
					genome_i = self.genomes[ddict["genome"]]
				except KeyError:
					logger.warning("No genome annotation could be found for given reference, please add {genome}!".format(genome=ddict["genome"]))

					_id = False
				logger.debug(ddict)
				if _id:
					self.add_annotation(ddict, id=_id, genome_i=genome_i)
		self.taxonomydb.commit()
		return dataArr

	def load_genome_reference_file(self,reference_file):
		'''Add new reference genomes from file!'''
		with open(reference_file) as af:
			headers = af.readline().lstrip("#").strip().split("\t")
			logging.debug(headers)
			check_source_headers(headers,GENOME_SOURCE_COLUMNS,reference_file)
			dataArr = []
			for row in af:
				data = row.strip().split("\t")
				if len(data) < len(headers): ## short row would run out of values in the header zip below (issue #89)
					logging.warning("Skipped short row in reference file: '{row}'".format(row=row.strip()))
					continue
				ddict = {}
				for i in range(len(headers)):
					head,val = headers[i],data[i]
					ddict[head] = val
				dataArr.append(ddict)
				genome_i = self.add_genome(ddict)
				self.genomes[ddict["genome"]] = genome_i
				self.genomes[ddict["genbank_id"]] = genome_i
				self.genomes[ddict["refseq_id"]] = genome_i
		self.taxonomydb.commit()
		return dataArr


	def update_database(self):
		'''Update the database file'''
		'''This override carried the old flextaxd attribute names (modified_links), which
			ModifyTree never sets, so every --mod_file run died with an AttributeError (issue #133),
			the names below follow ModifyTree.parse_modification
		'''
		removed_links = (self.old_links | self.non_overlapping_old_links) - set(self.parent_link)
		if self.replace:
			self.taxonomydb.delete_links(removed_links)
			self.taxonomydb.delete_nodes(self.old_nodes)
			'''the snp annotations of removed nodes have to go with them, leftovers would still be
				picked up by the annotation lookups (issue #134)'''
			for node_id in self.old_nodes:
				self.taxonomydb.query("DELETE FROM snp_annotation WHERE node_id = ?",insert_val=(node_id,),error=True)
			self.taxonomydb.commit()
			logger.info("Deleted nodes {nodes}".format(nodes=self.old_nodes))
		links,nodes = self.taxonomydb.add_links(self.new_links)
		logger.debug("Added nodes {nodes}, added links {links}".format(nodes=nodes,links=links))

		if len(links) + len(nodes) + len(removed_links) > 0:
			if self.replace:
				logger.info("Deleting {n} links and {n2} nodes that are no longer valid".format(n=len(removed_links),n2=len(self.old_nodes)))
			logger.info("Adding {n} new nodes".format(n=len(nodes)))
			logger.info("Adding {n} updated and/or new links".format(n=len(links)))

			''' Commit changes (only commit once both deletion and addition of new nodes and links are completed!)'''
			self.taxonomydb.commit()
			self.nodeDict = self.taxonomydb.get_nodes()
			if self.mod_genomes:
				logger.info("Transfering genomeid2taxid annotation from incoming database")
				self.update_genomes()
		else:
			logger.info("All updates already found in database, nothing has been changed!")
		return


	def remove_branch(self,node):
		'''--remove was documented but no code ever read the flag, so CanSNPer2-database silently did
			nothing (issue #147). Remove the branch starting at the given node: the node itself, all of
			its descendants, every link to or from them and their SNP annotations.
		'''
		try:
			node_id = self.taxonomydb.get_id(node)
		except NameNotFoundError:
			raise InputError("The node {node} could not be found in the database {database}!".format(node=node,database=self.taxonomydb.database))
		branch = self.taxonomydb.get_children(set([node_id])) | set([node_id])
		links = self.taxonomydb.get_links(branch)
		'''CanSNPer marks the root of a tree with a link from the node to itself, deleting that node
			would wipe the whole database'''
		if (node_id,node_id,node_id) in links or (node_id,node_id) in [(link[0],link[1]) for link in links]:
			raise InputError("{node} is a root of the tree, removing it would delete the whole database, remove a branch below it instead!".format(node=node))
		logger.info("Removing the branch from {node}: {n} nodes and {m} links".format(node=node,n=len(branch),m=len(links)))
		self.taxonomydb.delete_links(links,hold=True)
		'''the annotations of the removed nodes have to go with them, leftover rows would still be
			picked up by the annotation lookups (issue #134)'''
		for removed_id in branch:
			self.taxonomydb.query("DELETE FROM snp_annotation WHERE node_id = ?",insert_val=(removed_id,),error=True)
		self.taxonomydb.delete_nodes(branch,hold=True)
		self.taxonomydb.commit()
		logger.info("Branch from {node} removed!".format(node=node))
		return branch
