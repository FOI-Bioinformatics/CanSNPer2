import sys, os, inspect
import sqlite3

#from CanSNPer2.modules.DatabaseConnection import DatabaseConnection


import logging
logger = logging.getLogger(__name__)
from flextaxd.modules.database.CreateDatabase import CreateDatabase
from flextaxd.modules.database.DatabaseConnection import DatabaseFunctions
from importlib import import_module
from datetime import date
#script_path = os.path.dirname(os.path.abspath(__file__))  ## Retrieve the abspath of the location of this script

def dynamic_import(abs_module_path, class_name):
	module = import_module(".".join([abs_module_path,class_name]))
	target_class = getattr(module, class_name)
	return target_class

class CreateSNPDatabase(CreateDatabase):
	'''
		InitiateDatabase creates an empty CanSNPer2 database with all nessesary tables
		Does not include any functionality for modifying or filling the database and is
		written to be independend of other database object!
		Extends CreateDatabase from FlexTaxD
	'''
	def __init__(self, database,create=False,verbose=False):
		super().__init__(verbose=verbose)  ## This should create the batabase and all its tablestable
		#self.database = database
		'''If database does not exist create database, else nothing will happen'''
		if create:
			self.create_database(database)
			logger.info("Create CanSNPer2 table")
			self.create_cansnp_table()

	def create_cansnp_table(self):
		sql_create_SNP_table =      """CREATE TABLE IF NOT EXISTS snp_annotation (
									node_id INTEGER PRIMARY KEY,
									snp_id VARCHAR(6) ,
									position INTEGER,
									ancestral_base VARCHAR(1),
									derived_base VARCHAR(1),
									reference VARCHAR(20),
									date DATETIME,
									genome_i INTEGER,
									FOREIGN KEY (genome_i) REFERENCES genomes (genome_i)
								);
							"""
		sql_create_snp_references_table =  """CREATE TABLE IF NOT EXISTS snp_references (
									id INTEGER PRIMARY KEY,
									genome VARCHAR(6),
									strain VARCHAR(200),
									genbank_id VARCHAR(20),
									refseq_id VARCHAR(20),
									assembly_name VARCHAR(200)
								);
							"""
		logger.info("Add SNP table")
		self.add_table(sql_create_SNP_table)
		logger.info("Add SNP reference table")
		self.add_table(sql_create_snp_references_table)
		return

class CanSNPDatabase(DatabaseFunctions):
	"""Class containing functions for Loading data into the database, inherits functions from flextaxd database"""
	def __init__(self,*args,**kwargs):
		super().__init__(kwargs["database"],verbose=kwargs["verbose"])
		logger.debug(kwargs)
		'''Connect to database'''
		if kwargs["create"]:
			logger.debug("Create: {create}".format(create=kwargs["create"]))
			CanSNPdatabase_obj = CreateSNPDatabase(kwargs["database"],create=kwargs["create"],verbose=kwargs["verbose"])
		else:
			logger.debug("Connect to database {db}".format(db=kwargs["database"]))
		#CanSNPdatabase_obj.create_cansnptable(kwargs["create"])
		source_type = kwargs["source_type"]
		self.today = date.today()
		if kwargs["tree"]:
			logger.debug("Loading module ReadTaxonomy{type}".format(type=source_type))
			read_module = dynamic_import("flextaxd.modules", "ReadTaxonomy{type}".format(type=source_type))
			read_obj = read_module(kwargs["tree"], database=kwargs["database"],root_name=False,rank="family",verbose=kwargs["verbose"])
			read_obj.verbose = kwargs["verbose"]                                                          ## Set verbose mode
			logger.info("Parse taxonomy")
			read_obj.parse_taxonomy()                                                           ## Parse taxonomy file
		if kwargs["references"]:
			logger.info("Load genome reference file!")
			self.load_genome_annotation(kwargs["references"])
		if kwargs["annotation"]:
			logger.info("Load annotation file!")
			self.load_cansnp_annotation(kwargs["annotation"])

	def add_annotation(self, data, id=False,genome_i=False):
		'''Add snp annotation to annotation table'''
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
		taxid_base = self.insert(info, table="snp_annotation")
		return taxid_base

	def add_genome(self,data):
		'''Add genome annotation'''
		logger.debug("Add genome")
		logger.debug(data)
		info = {
				"genome":            data["genome"],
				"strain":         data["strain"],
				"genbank_id":     data["genbank_id"],
				"refseq_id":        data["refseq_id"],
				"assembly_name":         data["assembly_name"],
				#"comment":             data["comment"]
		}
		return self.insert(info, table="snp_references")

	def load_genome_annotation_file(self,annotation_file):
		self.genomes = {}
		with open(annotation_file) as f:
			headers = f.readline().lstrip("#").strip().split("\t")
			#dataArr = []
			logging.debug(headers)
			for row in f:
				data = row.strip().split("\t")
				ddict = {}
				for i in range(len(headers)):
					head,val = headers[i],data[i]
					ddict[head] = val
				genome_i = self.add_genome(ddict)
				self.genomes[ddict["genome"]] = genome_i
				self.genomes[ddict["genbank_id"]] = genome_i
				self.genomes[ddict["refseq_id"]] = genome_i
		self.conn.commit()
		return self.genomes

	def load_annotation_file(self,annotation_file):
		'''Read the annotation source file'''
		nodes = self.get_nodes()
		logger.info(nodes)
		logger.debug("nodes: {nnodes}".format(nnodes=len(nodes)/2))
		if len(self.genomes) == 0:
			self.genomes = self.get_genomes(table="snp_references")
		logger.debug(self.genomes)
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
					logger.warning("No genome annotation could be found add a genome annotation file to command!")
					logger.debug(ddict)
					_id = False
				if _id:
					logger.debug("Add annotation: {dic}, {_id}, {gid}".format(dic=ddict,_id=_id,gid=genome_i))
					self.add_annotation(ddict, id=_id, genome_i=genome_i)
		### commit changes to database
		self.conn.commit()
		return dataArr

	def load_cansnp_annotation(self,annotation_file):
		'''Function that loads data into the cansnp annotation table'''
		dbConn = dynamic_import("flextaxd.modules.database", "DatabaseConnection")
		data = self.load_annotation_file(annotation_file)
		logging.info("CanSNPer annotation loaded!")

	def load_genome_annotation(self,annotation_file):
		'''Function that loads genome data into the cansnp genome annotation table'''
		dbConn = dynamic_import("flextaxd.modules.database", "DatabaseConnection")
		data = self.load_genome_annotation_file(annotation_file)
		logging.debug(data)
		logging.info("Genome annotation loaded!")
