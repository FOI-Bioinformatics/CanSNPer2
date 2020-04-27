'''
	Modify and update the CanSNPer2 database
'''

from flextaxd.modules.ModifyTree import ModifyTree
from flextaxd.modules.database.DatabaseConnection import DatabaseFunctions
from datetime import date
import logging
logger = logging.getLogger(__name__)

class ModifyCanSNPer2Database(ModifyTree):
	"""ModifyCanSNPer2Database adds a few important functions only nessesary when modifying a database"""
	def __init__(self, database, mod_database=False, mod_file=False, separator="\t",verbose=False,parent=False,replace=False,snp_annotation=False,**kwargs):
		logger.info("Load ModifyFunctions")
		super().__init__(database, mod_database, mod_file, separator,verbose,parent,replace)
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
		return self.insert(info, table="snp_references")

	def load_annotation_file(self,annotation_file):
		'''Read the annotation source file'''
		nodes = self.get_nodes()
		if len(self.genomes) == 0:
			self.genomes = self.get_genomes(table="snp_references")
		with open(annotation_file) as af:
			headers = af.readline().lstrip("#").strip().split("\t")
			logging.debug(headers)
			dataArr = []
			for row in af:
				data = row.strip().split("\t")
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
		logger.critical("The load reference modification function is not yet completed!")
		with open(annotation_file) as af:
			headers = af.readline().lstrip("#").strip().split("\t")
			logging.debug(headers)
			dataArr = []
			for row in af:
				data = row.strip().split("\t")
				ddict = {}
				for i in range(len(headers)):
					head,val = headers[i],data[i]
					ddict[head] = val
				dataArr.append(ddict)
		return


	def update_database(self):
		'''Update the database file'''
		if self.replace:
			self.taxonomydb.delete_links(self.modified_links)
			self.taxonomydb.delete_links(self.existing_links-set(self.parent_link))  ## delete existing links from old nodes, except parent
			self.taxonomydb.delete_nodes(self.old_nodes)
			logger.info("Deleted nodes {nodes}".format(nodes=self.old_nodes))
		links,nodes = self.taxonomydb.add_links(self.new_links)
		logger.debug("Added nodes {nodes}, added links {links}".format(nodes=nodes,links=links))

		if len(links) + len(nodes) + len(self.modified_links) > 0:
			if self.replace:
				logger.info("Deleting {n} links and {n2} nodes that are no longer valid".format(n=len(self.modified_links | self.existing_links),n2=len(self.old_nodes)))
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


#
# class AddReference(object):
#     """docstring for AddReference."""
#
#     def __init__(self, database, reference):
#         super(AddReference, self).__init__()
#
#     def add_genome(self,data):
#         '''Add genome annotation'''
#         logger.debug("Add genome")
#         logger.debug(data)
#         info = {
#                 "genome":            data["genome"],
#                 "strain":         data["strain"],
#                 "genbank_id":     data["genbank_id"],
#                 "refseq_id":        data["refseq_id"],
#                 "assembly_name":         data["assembly_name"],
#                 #"comment":             data["comment"]
#         }
#         return self.insert(info, table="snp_references")
#
#     def load_genome_annotation_file(self,annotation_file):
#         self.genomes = {}
#         with open(annotation_file) as f:
#             headers = f.readline().lstrip("#").strip().split("\t")
#             #dataArr = []
#             logging.debug(headers)
#             for row in f:
#                 data = row.strip().split("\t")
#                 ddict = {}
#                 for i in range(len(headers)):
#                     head,val = headers[i],data[i]
#                     ddict[head] = val
#                 genome_i = self.add_genome(ddict)
#                 self.genomes[ddict["genome"]] = genome_i
#                 self.genomes[ddict["genbank_id"]] = genome_i
#                 self.genomes[ddict["refseq_id"]] = genome_i
#         self.conn.commit()
