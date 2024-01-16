import sys, os
import sqlite3

'''
DatabaseConnections and CanSNPdbFunctions are classes written to simplify database work using sqlite3
This particular script is specially adapted to suit the use for CanSNPer web database!
'''

#__name__="DatabaseConnection"
__version__ = "0.1.5"
__author__ = "David Sundell"
__credits__ = ["David Sundell"]
__license__ = "GPLv3"
__maintainer__ = "FOI bioinformatics group"
__email__ = ["bioinformatics@foi.se", "david.sundell@foi.se"]
__date__ = "2019-04-17"
__status__ = "Production"
__partof__ = "CanSNPer2"

import logging
logger = logging.getLogger(__name__)

class ConnectionError(Exception):
	def __init__(self, value):
		self.value = value
	def __str__(self):
		return repr(self.value)

class DatabaseConnection(object):
	"""docstring for DatabaseConnection"""
	def __init__(self, database, verbose=False):
		super().__init__()
		self.verbose = verbose
		self.database = database
		logging.debug("DatabaseConnection variables database: {database}, verbose: {verbose}".format(database=database,verbose=verbose))
		self.conn = self.connect(self.database)
		self.cursor = self.create_cursor(self.conn)
		logging.debug("Database connected!")

	def __str__(self):
		return "Object of class DatabaseConnection, connected to {database}".format(database=self.database)

	def __repr__(self):
		return "DatabaseConnection()"

	def set_verbose(self,val):
		self.verbose = val

	def connect(self,database):
		'''Create database connection'''
		try:
			conn = sqlite3.connect(database)
			logging.info("{database} opened successfully... ".format(database=database))
			return conn
		except Exception as e:
			sys.stderr.write(str(e))
		raise ConnectionError("Count not connect to the database {database} see above message for details!".format(database=database))

	def disconnect(self):
		self.conn.close()

	def create_cursor(self,conn):
		'''Create a db cursor'''
		try:
			cursor = conn.cursor()
			logging.info("cursor created")
			return cursor
		except Exception as e:
			sys.stderr.write(str(e))
		return None

	def commit(self):
		self.conn.commit()

	def query(self,query,insert_val = False,getres=False, cursor=False,error=False):
		'''    The query function controls  the error handling of sqlite3 execute'''
		res = False
		if not cursor:
			cursor = self.cursor
		try:
			if insert_val:
				res = cursor.execute(query,insert_val)
				if error or getres:
					return res
				return cursor.lastrowid
			else:
				return cursor.execute(query)
		except Exception as e:
			if "UNIQUE constraint failed" not in str(e):
				## UNIQUE constraint is an accepted error as it keeps multiple edges from being added
				logging.error("Error in DatabaseConnection query")
				logging.error(query)
				logging.error("Insert val", insert_val)
				logging.debug(str(e))
				raise
			return(e)

	def insert(self,data,table):
		'''Insert function
				data is a dictionary with keys matching
				table columns with respective value to be inserted
		'''
		INSERT_QUERY = '''
			INSERT INTO {table}({columns})
			VALUES ({values})
		'''

		columns = data.keys()
		values = tuple([data[col] for col in columns])
		insertStr = INSERT_QUERY.format(
				columns=",".join(columns),
				values=','.join(["?" for x in values]),
				table=table
		)
		return self.query(insertStr,insert_val=values)

	def update(self,data,table):
		'''Update function requires table column which column to identify row with and value to replace'''
		columns = []
		colkeys = list(data["set_column"].keys())
		for col in colkeys:
			if col != "snp_id":
				columns.append("{set_column}='{val}'".format(set_column=col, val=data["set_column"][col]))
		set_query = ", ".join(columns)
		UPDATE_QUERY = '''
			UPDATE {table}
			SET {columns}
			WHERE {where_column} = ?
		'''.format(table=table,columns=set_query,where_column=data["where_column"])
		###  UPDATE genomes SET id = newnode WHERE genome = oldname
		logging.debug(UPDATE_QUERY,data)
		datasend = tuple([data["set_value"]])
		return self.query(UPDATE_QUERY,datasend,getres=True)

class LoadCanSNPAnnotation(DatabaseConnection):
	def __init__(self, database, verbose=False):
		super().__init__(database,verbose)
		logger.info.debug(self.database)
		if not os.path.exists(self.database):
			logging.error(ConnectionError("The database {database} does not exist, please create a FlexTaxDatabase or check your database path!".format(self.database)))

	def getSNPdict(self, table=False,columns="*"):
		'''Get full table from table'''
		logging.debug(self.database)
		QUERY = '''SELECT {columns} FROM {table}'''.format(table=table,columns=columns)
		logging.info(QUERY)
		data = self.database.query(QUERY).fetchall()
		dic = {}
		for id,name in data:
			dic[name] = id
		return dic

	def load_annotations(self,annotation):
		'''Function to load data into CanSNPer2 database table'''
		self.snpDict = self.getSNPdict(table="nodes",columns="id,name")
		if not annotation:
			logging.error(ConnectionError("The annotation file does not exist!".format(self.database)))

		with open(annotation) as f:
			headers = f.readline().strip().lstrip("#").split("\t")
			'''Define new headers as they must match the database'''
			headers = ["node_id", "snp_id", "position", "ancestral_base","derived_base","reference","date","genome_i"]
			for row in f:
				data = row.strip().split("\t")
				self.add_node(data, headers)
				print(headers, data)
				exit()

	def add_node(self, data, headers, table="snp_annotation"):
		'''Add node to tree'''
		'''Must get node ID from FlexTaxDatabase!!'''

		info = { "node_id", self.snpDict[data.pop(0)]}
		### If ID is supplied skip autoincrement and add specific ID
		for head in headers:
			info[head] = data.pop(0)
		logger.debug(info, table)
		taxid_base = self.insert(info, table=table)
		return taxid_base


class CanSNPdbFunctions(DatabaseConnection):
	"""CanSNPerdb database function class contains multiple additional database
		functions to simplify data access related to the website, its a subclass of DatabaseConnection"""
	def __init__(self, database,  verbose=False):
		super().__init__(database,verbose)
		logging.info("Load CanSNPdbFunctions")
		### store DatabaseConnection object reference

	def __str__(self):
		return "Object of class CanSNPdbFunctions, connected to {database}".format(database=self.database)

	'''
		Set functions
	'''
	def set_database(self, database):
		'''Change the database object default in class'''
		self.database = database
		self.conn = self.connect(self.database)
		self.cursor = self.create_cursor(self.conn)

	'''
		Get functions of class
	'''
	def get_taxid_base(self):
		'''Fetch the next incremental node from the current database'''
		QUERY = "SELECT MAX(id) AS max_id FROM nodes"
		return self.query(QUERY).fetchone()[0]+1

	def get_genomes(self, database=False,limit=0):
		'''Get the list of genomes in the database'''
		## This is a many to many relation, so all genomes has to be put in a set for each taxonomy id
		genomeDict = {}
		QUERY = '''SELECT id,genome FROM {table}'''.format(table="genomes")
		if limit > 0:
			QUERY += " LIMIT {limit}".format(limit=limit)
		#print(QUERY)
		for id,genome in database.query(QUERY).fetchall():
			genomeDict[genome] = id
		return genomeDict

	def get_all(self, database=False, table=False,columns="*"):
		'''Get full table from table'''
		if not database:
			database = self.database
		QUERY = '''SELECT {columns} FROM {table}'''.format(table=table,columns=columns)
		logging.info(QUERY)
		return database.query(QUERY).fetchall()

	def get_nodes(self, database=False,col=False):
		'''Retrieve the whole node info table of the database to decrease the number of database calls!'''
		nodeDict = {}
		QUERY = '''SELECT id,name FROM nodes'''
		if not database:
			database = self
		for node in database.query(QUERY).fetchall():
			nodeDict[node[0]] = node[1]
			if col == 1:
				continue
			nodeDict[node[1]] = node[0]
		return nodeDict

	def get_links(self, nodes,database=False):
		'''This function returns all links in the given database'''
		QUERY = '''SELECT parent,child FROM tree WHERE parent in ({nodes}) OR child in ({nodes})'''.format(nodes=",".join(map(str,nodes)))
		if not database:
			database = self
		links = database.query(QUERY).fetchall()
		return links

	'''
		Add functions of class
	'''
	def add_node(self, description, id=False,table="nodes"):
		'''Add node to tree'''
		info = { "name": description }
		### If ID is supplied skip autoincrement and add specific ID
		if id:
			info["node_i"] = id
		taxid_base = self.insert(info, table=table)
		return taxid_base

	def add_link(self, child, parent,rank=False, table="tree"):
		'''Add relationship in tree'''
		info = {
			"child": child,
			"parent": parent
		}
		return self.insert(info, table=table)

	def add_genome(self, id, table="genomes"):
		'''Add genome annotation to nodes'''
		info = {
			"node_i": id,
			"genome": genome
		}
		return self.insert(info, table=table)

	def add_rank(self, rank,id=False):
		'''Add node to tree'''
		info = { "rank": rank }
		### If ID is supplied skip autoincrement and add specific ID
		if id:
			info["id"] = id
		rank_id = self.insert(info, table="rank")
		if self.verbose: print("rank added: ",info)
		return rank_id

	def add_links(self,links, table="tree",hold=False):
		'''Add links from a list to tree'''
		added_links = 0
		for parent,child in links:
			res = self.add_link(parent,child,table=table)
			### Check if the link already exist in the database, this overlap may occur when a large new branch is added
			if "UNIQUE constraint failed" not in str(res):
				added_links +=1
		## Commit changes
		if not hold:
			self.commit()
		return added_links

	def add_nodes(self,nodes, table="Tree",hold=False):
		'''Add nodes from a list of nodes'''
		added_nodes = 0
		for node in nodes:
			res = self.add_node(node)
			### Check if the link already exist in the database, this overlap may occur when a large new branch is added
			if "UNIQUE constraint failed" not in str(res):
				added_nodes +=1
		## Commit changes
		if not hold:
			self.commit()
		return added_nodes

	'''
		Delete functions of class
	'''
	def delete_links(self,links, table="Tree",hold=False):
		'''This function deletes all links given in links'''
		QUERY = "DELETE FROM {table} WHERE parent = {parent} AND child = {child}"
		for parent,child in links:
			res = self.query(QUERY.format(table=table, parent=parent, child=child))
		## Commit changes
		if not hold:
			self.commit()

	def delete_nodes(self, nodes, table="nodes",hold=False):
		'''This function deletes all nodes given in nodes'''
		QUERY = "DELETE FROM {table} WHERE node_i = {node}"
		for node in nodes:
			res = self.query(QUERY.format(table=table, node=node))
		## Commit changes
		if not hold:
			self.commit()

	def num_rows(self,table):
		'''Return the number of rows in a table'''
		QUERY = '''SELECT Count(*) FROM {table}'''
		return self.query(QUERY.format(table=table)).fetchall()[0][0]

class XMFAFunctions(DatabaseConnection):
	"""CanSNPerdb database function class contains multiple additional database
		functions to simplify data access related to the website, its a subclass of DatabaseConnection"""
	def __init__(self, database, verbose=False):
		import time
		super().__init__(database,verbose)
		logging.info("Load XMFAFunctions")
		### store DatabaseConnection object reference

	'''
		XMFAfunctions
	'''
	def get_results(self,snp_string,reference):
		SNPs = {}
		res = self.query(snp_string, (reference,),getres=True)
		for strain, pos,tbase,rbase,SNP in res.fetchall():
			# check if this position was already added and if so throw an error (cansnper does not save SNPs by their ID but rather by their position in each reference)
			if pos in SNPs:
				logger.error("Error: there was a duplicate CanSNP in your CanSNPer database at position "+str(pos)+" in reference "+reference)
				logger.error("Please check your database input files")
				raise Exception("Multiple SNPs with same position")
			#/
			SNPs[pos] = tuple([pos,rbase, tbase,SNP])
		return SNPs

	def get_snps(self, reference="SCHUS4.2"):
		'''Returns a list of all SNPs and their positions.
		Keyword arguments:
		returns: results as a dictionary with tuple SNP for each position {pos: (pos, refBase, TargetBase, SNPid)}
				 and a list of positions sorted ASC
		'''
		snp_string = """SELECT genome, position, derived_base, ancestral_base, snp_id
										FROM snp_annotation
										LEFT JOIN snp_references on (snp_references.id = snp_annotation.genome_i)
										WHERE genome = ?
									"""
		# res = self.query(snp_string, (reference,),getres=True)
		# for strain, pos,tbase,rbase,SNP in res.fetchall():
		# 	SNPs[pos] = tuple([pos,rbase, tbase,SNP])
		SNPs = self.get_results(snp_string,reference)
		if len(SNPs) == 0:
			logger.error("No SNPs could be found in database check your database connection!")
			logger.debug(snp_string)
			raise ValueError("No SNP's was found in the database for reference: {reference}".format(reference=reference))
		snp_positions = list(SNPs.keys())
		snp_positions.sort()
		return SNPs,snp_positions

	# def get_references(self):
	# 	'''Finds all references in the database and returns a list'''
	# 	query = """SELECT DISTINCT(Strain) FROM Sequences"""
	# 	return [x[0] for x in self.query(query).fetchall()]
