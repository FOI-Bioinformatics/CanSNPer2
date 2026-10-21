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
				logging.error("Insert val: {vals}".format(vals=insert_val))
				logging.debug(str(e))
				raise
			if getres or error or not insert_val:
				## UNIQUE is only swallowed for the plain inserts where callers test str(res) for it,
				## when a cursor was requested the exception must propagate instead of breaking fetch* chains
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
		setvals = []
		colkeys = list(data["set_column"].keys())
		for col in colkeys:
			if col != "snp_id":
				columns.append("{set_column}=?".format(set_column=col))
				setvals.append(data["set_column"][col])
		set_query = ", ".join(columns)
		UPDATE_QUERY = '''
			UPDATE {table}
			SET {columns}
			WHERE {where_column} = ?
		'''.format(table=table,columns=set_query,where_column=data["where_column"])
		###  UPDATE genomes SET id = newnode WHERE genome = oldname
		logging.debug("{query} {data}".format(query=UPDATE_QUERY,data=data))
		### the WHERE value is data["where"] when supplied, data["set_value"] only as legacy fallback (issue #65)
		where_val = data["where"] if "where" in data else data["set_value"]
		datasend = tuple(setvals + [where_val])
		return self.query(UPDATE_QUERY,datasend,getres=True)

class LoadCanSNPAnnotation(DatabaseConnection):
	def __init__(self, database, verbose=False):
		super().__init__(database,verbose)
		logger.debug("LoadCanSNPAnnotation connected to {database}".format(database=self.database))
		if not os.path.exists(self.database):
			logging.error(ConnectionError("The database {database} does not exist, please create a FlexTaxDatabase or check your database path!".format(database=self.database)))

	def getSNPdict(self, table=False,columns="*"):
		'''Get full table from table'''
		if not table:
			raise Exception("Table was not supplied to LoadCanSNPAnnotation getSNPdict()")
		logging.debug(self.database)
		QUERY = '''SELECT {columns} FROM {table}'''.format(table=table,columns=columns)
		logging.info(QUERY)
		data = self.query(QUERY).fetchall()
		if data and len(data[0]) != 2: ## the result is unpacked as (id,name), more or fewer columns cannot be turned into a dict (issue #75)
			raise ValueError("getSNPdict expects exactly two columns (id, name) from table {table}, got {n}".format(table=table,n=len(data[0])))
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

		info = { "node_id": self.snpDict[data.pop(0)]}
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
		res = self.query(QUERY).fetchone()
		if res is None or res[0] is None:
			return 1
		return res[0]+1

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
			info["id"] = id
		taxid_base = self.insert(info, table=table)
		return taxid_base

	def add_link(self, child, parent,rank=False, table="tree"):
		'''Add relationship in tree'''
		info = {
			"child": child,
			"parent": parent
		}
		return self.insert(info, table=table)

	def add_genome(self, id, genome, table="genomes"):
		'''Add genome annotation to nodes'''
		info = {
			"id": id,
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
		QUERY = "DELETE FROM {table} WHERE parent = ? AND child = ?"
		for parent,child in links:
			res = self.query(QUERY.format(table=table),insert_val=(parent,child),getres=True)
		## Commit changes
		if not hold:
			self.commit()

	def delete_nodes(self, nodes, table="nodes",hold=False):
		'''This function deletes all nodes given in nodes'''
		QUERY = "DELETE FROM {table} WHERE id = ?"
		for node in nodes:
			res = self.query(QUERY.format(table=table),insert_val=(node,),getres=True)
		## Commit changes
		if not hold:
			self.commit()

	def num_rows(self,table):
		'''Return the number of rows in a table'''
		QUERY = '''SELECT Count(*) FROM {table}'''
		return self.query(QUERY.format(table=table)).fetchall()[0][0]

class NameNotFoundError(Exception):
	def __init__(self, value):
		self.value = value
	def __str__(self):
		return repr(self.value)

class TreeError(Exception):
	def __init__(self, value):
		self.value = value
	def __str__(self):
		return repr(self.value)

'''
	CreateDatabase, TaxonomyDBConnection, DatabaseFunctions and ModifyFunctions are
	portable replacements for the former flextaxd.modules.database dependencies
	(ISSUES #59), CanSNPer2 now carries its own taxonomy database schema and functions.
'''

class CreateDatabase(object):
	"""Creates an empty taxonomy database with the nodes, tree, rank and genomes tables"""
	def __init__(self, verbose=False):
		super().__init__()
		self.verbose=verbose
		self.sql_create_nodes_table = """ CREATE TABLE IF NOT EXISTS nodes (
											id integer PRIMARY KEY,
											name text NOT NULL
										); """

		self.sql_create_tree_table = """CREATE TABLE IF NOT EXISTS tree (
										parent integer NOT NULL,
										child integer NOT NULL,
										rank_i integer,
										FOREIGN KEY (parent) REFERENCES nodes (id),
										FOREIGN KEY (child) REFERENCES nodes (id),
										FOREIGN KEY (rank_i) REFERENCES rank (rank_i),
										unique (parent, child,rank_i)
									);"""

		self.sql_create_rank_table = """CREATE TABLE IF NOT EXISTS rank (
										rank_i integer PRIMARY KEY,
										rank VARCHAR(15)
									);"""

		self.sql_create_genomes_table = """ CREATE TABLE IF NOT EXISTS genomes (
											id integer NOT NULL,
											genome VARCHAR(255) NOT NULL UNIQUE,
											reference VARCHAR(255),
											FOREIGN KEY (id) REFERENCES nodes (id)
										); """

	def create_connection(self,db_file):
		'''Create a database connection to the sqlite database specified by db_file'''
		try:
			conn = sqlite3.connect(db_file)
			return conn
		except Exception as e:
			logger.debug(e)
		return None

	def create_table(self,create_table_sql):
		'''Create a table from the create_table_sql statement'''
		try:
			c = self.conn.cursor()
			c.execute(create_table_sql)
		except Exception as e:
			logger.debug(e)

	def add_table(self,table):
		self.create_table(table)
		self.conn.commit()
		return

	def create_database(self,database=False):
		'''Connect to the database and create all taxonomy tables'''
		self.conn = self.create_connection(database)
		if self.conn is not None:
			self.create_table(self.sql_create_nodes_table)
			self.create_table(self.sql_create_tree_table)
			self.create_table(self.sql_create_genomes_table)
			self.create_table(self.sql_create_rank_table)
			self.conn.commit()
		else:
			raise ConnectionError("Error! cannot connet to the database {db}".format(db=database))
		return

class TaxonomyDBConnection(CreateDatabase):
	'''Connection to a taxonomy database, replaces the flextaxd DatabaseConnection class,
		creates a missing database with CreateDatabase instead of calling it as a script'''
	def __init__(self, database, verbose=False):
		super().__init__(verbose=verbose)
		self.verbose = verbose
		self.database = database
		if not os.path.exists(self.database):
			logger.info("Database {database} does not exist, creating it!".format(database=self.database))
			self.create_database(self.database)
		try: ## If database connection already exists
			self.conn
		except AttributeError:
			logger.debug("Connecting to {database}".format(database=self.database))
			self.conn = self.create_connection(self.database)
		if self.conn is None:
			raise ConnectionError("Count not connect to the database {database} see above message for details!".format(database=self.database))
		self.cursor = self.create_cursor(self.conn)

	def __str__(self):
		return "Object of class TaxonomyDBConnection, connected to {database}".format(database=self.database)

	def __repr__(self):
		return "TaxonomyDBConnection()"

	def set_verbose(self,val):
		'''Set the verbose level of the connection object'''
		self.verbose = val
		return self.verbose

	def create_cursor(self,conn):
		'''Create a db cursor'''
		try:
			self.cursor = conn.cursor()
			logger.debug("cursor created.")
			return self.cursor
		except Exception as e:
			sys.stderr.write(str(e))
		raise ConnectionError("Count not create db cursor to {database}".format(database=self.database))

	def commit(self):
		self.conn.commit()

	def query(self,query,insert_val = False, cursor=False,error=False):
		'''The query function is a wrapper around sqlite3 execute to form responses related to the request

		Returns
			int - 	if insert: lastrowid
			list -	if request: sqlite3 execute result
			str - 	if Exception: returs error message
		'''
		res = False
		if not cursor:
			cursor = self.cursor
		try:
			if insert_val:
				res = cursor.execute(query,insert_val)
				if error:
					return res
				return cursor.lastrowid
			else:
				return cursor.execute(query)
		except Exception as e:
			if "UNIQUE constraint failed" not in str(e):
				## UNIQUE constraint is an accepted error as it keeps multiple edges from being added
				logger.warning("Error in TaxonomyDBConnection query")
				logger.warning(query)
				if insert_val:
					logger.debug("Insert val: {vals}".format(vals=insert_val))
				sys.stderr.write(str(e)+"\n")
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
		UPDATE_QUERY = '''
			UPDATE {table}
			SET {set_column} = ?
			WHERE {where_column} = ?
		'''.format(table=table,set_column=data["set_column"],where_column=data["where_column"])
		udata = tuple([data["set_value"],data["where"]])
		logger.debug("{q} {d}".format(q=UPDATE_QUERY,d=udata))
		res = self.query(UPDATE_QUERY,udata,error=True)
		if self.rowcount() == 0:
			res = self.insert({data["set_column"]:data["set_value"],data["where_column"]:data["where"]}, table)
			return False
		return True

	def multi_update(self,data,table):
		'''Update function from a list of (set_value, where_value) pairs.
			Updates existing rows and inserts the missing ones, values are bound
			as parameters so no UNIQUE constraint on the where_column is required (issue #66)'''
		set_column = data["set_column"]
		where_column = data["where_column"]
		UPDATE_QUERY = '''
			UPDATE {table}
			SET {set_column} = ?
			WHERE {where_column} = ?
		'''.format(table=table,set_column=set_column,where_column=where_column)
		logger.debug("{q}".format(q=UPDATE_QUERY))
		updated = 0
		added = 0
		for set_val,where_val in data["data"]:
			self.query(UPDATE_QUERY,insert_val=(set_val,where_val))
			if self.rowcount() == 0:
				self.insert({set_column:set_val,where_column:where_val}, table)
				added += 1
			else:
				updated += 1
		logger.info("{updated} updated and {added} added in {table}".format(updated=updated,added=added,table=table))
		return True

	def delete(self,nodes,table):
		'''Deleting a node should make sure all related data is also removed'''
		DELETE_QUERY = '''
				DELETE FROM {table}
					WHERE id in({values})
		'''.format(table=table,
					values=','.join(["?" for x in nodes])
		)
		return self.query(DELETE_QUERY,insert_val=tuple(nodes))

	def rowcount(self):
		'''Get the number row number the cursor is currently at'''
		return self.cursor.rowcount

class DatabaseFunctions(TaxonomyDBConnection):
	"""DatabaseFunctions class defines additional functions to the TaxonomyDBConnection class
	"""
	def __init__(self, database, verbose=False):
		super().__init__(database, verbose)
		logger.debug("Load DatabaseFunctions")

	'''Validate tree function'''
	def validate_tree(self):
		'''This function validates the tree structure in the databases
			1. All nodes must only have one parent (check_parent)
			2. All edges must be attatched to the tree (match tree links and edges)
			3. All nodes must be attatched to the tree (match tree nodes and nodes)
			4. All nodes from links must have a node description (match the number of nodes from links with number of annotated nodes)

		------
		Returns
			True: if all edges has one and only one parent
		'''
		self.taxonomy = self.get_nodes()

		logger.info("Get all database nodes")
		self.nodes = self.get_nodes(col=1)					## get all nodes
		logger.info("Get all database edges")
		self.edges = self.get_links()					## Get all links
		logger.info("Get all children from root node")
		self.tree_connections = self.get_children([1])  ## Get all chilren from root
		logger.info("Get tree edges from children")
		self.tree_links = self.get_links(self.tree_connections) ## get all links from rooted tree
		logger.info("Get nodes from tree edges")
		self.tree_nodes = set()
		for edge in self.edges:
			self.tree_nodes.add(edge[0])
			self.tree_nodes.add(edge[1])
		logger.info("Validate parents")
		res = self.check_parent()
		extra = ""
		if len(res) > 0:
			failed_nodes = [self.nodes[x[0]] for x in list(res)]
			logger.debug(failed_nodes)
			if len(failed_nodes) > 10:
				logger.info(failed_nodes[:10])
				failed_nodes = len(failed_nodes)
				extra = "total failed nodes: "
			raise TreeError("There are nodes with multiple parents {extra}{nodes}".format(nodes=failed_nodes,extra=extra))
		stats = """Tree statistics
					Nodes: {nodes}
					Links: {links}
					Tree: n({tnodes}), l({tlinks})
					LinkNodes: {link_nodes}
					Parent_ok: True
					""".format(nodes = len(self.nodes),
								links=len(self.edges),
								tnodes=len(self.tree_connections),
								tlinks = len(self.tree_links),
								link_nodes = len(self.tree_nodes)
		)
		logger.info(stats)
		nodeset = set(self.nodes.keys())
		if len(nodeset - self.tree_connections) != 0:
			lset = nodeset - self.tree_connections
			if len(lset) > 10:
				lset = len(lset)
			logger.info("{}".format(lset))
			raise TreeError("The number of nodes and the number nodes under root does not match!")
		if len(set(self.edges) - set(self.tree_links)) != 0:
			lset = set(self.edges) - set(self.tree_links)
			if len(lset) > 10:
				lset = len(lset)
			logger.info("{}".format(lset))
			raise TreeError("The number of edges and the number edges under root does not match!")
		if len(self.tree_nodes - nodeset) != 0:
			lset = self.tree_nodes - nodeset
			if len(lset) > 10:
				lset = len(lset)
			logger.info("{}".format(lset))
			raise TreeError("The number of annotated nodes does not match with the number of nodes connected to edges!")
		if len(self.tree_connections - nodeset) != 0:
			lset = self.tree_connections - nodeset
			if len(lset) > 10:
				lset = len(lset)
			logger.info("{}".format(lset))
			raise TreeError("There are nodes in the database not connected to the tree")
		logger.info("Validation OK!")
		return True

	def statistics(self):
		'''Print statistics of the database'''
		logger.info("Get all database nodes")
		self.nodes = self.get_nodes(col=1)					## get all nodes
		logger.info("Get all database edges")
		self.edges = self.get_links()					## Get all links
		logger.info("Get annotated genomes")
		self.genomes = self.get_genomes()

		stats = """Tree statistics
					Nodes: {nodes}
					Links: {links}
					Genomes: {ngenomes}
					""".format(nodes = len(self.nodes),
								links=len(self.edges),
								ngenomes=len(self.genomes)
		)
		print(stats)
		return True

	def check_parent(self):
		'''check tree structure parent, make sure no tree node has multiple parents'''
		nrank = 0
		try:
			QUERY = "SELECT child FROM tree GROUP BY child HAVING count(parent) > 1"  ## Thanks to andrewjmc@github for this suggestion
			logger.debug(QUERY)
			child_w_dp = self.query(QUERY).fetchall()
			logger.debug(child_w_dp)
			## Bugfix  only allowed one parent to be different, now multiple is ok as long as ranks are different
			for child in map(str,[x[0] for x in child_w_dp]):
				QUERY = "SELECT parent,rank_i FROM tree WHERE child in ({children})".format(children=child)
				p_ranks = self.query(QUERY).fetchall()
				if len(p_ranks) != len(set(p_ranks)):
					logger.error("Nodes with two parents have identical ranks, Fatal Error")
				else:
					nrank += 1
			if nrank > 0:
				logger.info("Found: {n} nodes with multiple parents, however in different lineages, OK.".format(n=nrank))
				return []
			return child_w_dp
		except AttributeError as ae:
			logger.info("AttributeError occured")
			logger.info(QUERY)
			logger.info(child_w_dp)
			logger.debug(ae)
		except TypeError as ae:
			logger.info("TypeError occured")
			logger.info(QUERY)
			logger.info(child_w_dp)
			logger.debug(ae)
		raise TreeError("Following Nodes: {node} has more than one parent!".format(node=child_w_dp))

	'''Get functions of class'''
	def get_all(self, database=False, table=False,sort=False):
		'''Get full table from table'''
		if not table:
			raise Exception("Table was not supplied to DatabaseFuction get_all()")
		QUERY = '''SELECT * FROM {table}'''.format(table=table)
		if sort:
			QUERY += " ORDER BY {col} ASC".format(col=sort)
		return self.query(QUERY).fetchall()

	def get_taxid_base(self):
		'''Fetch the next incremental node from the current database'''
		QUERY = "SELECT MAX(id) AS max_id FROM nodes"
		res = self.query(QUERY).fetchone()[0]
		if res == None:
			res = 1
		else:
			res += 1
		return res

	def get_genomes(self, database=False,limit=0,table="genomes",cols="id,genome"):
		'''Get the list of genomes in the database'''
		## This is a many to many relation, so all genomes has to be put in a set for each taxonomy id
		genomeDict = {}
		QUERY = '''SELECT {cols} FROM {table}'''.format(cols=cols,table=table)
		if limit > 0:
			QUERY += " LIMIT {limit}".format(limit=limit)
		for id,genome in self.query(QUERY).fetchall():
			genomeDict[genome] = id
		return genomeDict

	def get_node(self, name, database=False):
		'''Retrieve one node by name'''
		QUERY = '''SELECT id FROM nodes where name = ?'''
		if not database:
			database = self.database
		return [n[0] for n in self.query(QUERY,insert_val=(name,),error=True).fetchall()]

	def get_nodes(self, database=False,col=5):
		'''Retrieve the whole node info table of the database to decrease the number of database calls!'''
		nodeDict = {}
		QUERY = '''SELECT id,name FROM nodes'''
		logger.debug(QUERY)
		if not database:
			database = self.database
		for node in self.query(QUERY).fetchall():
			if col == 1 or col == 5:
				nodeDict[node[0]] = node[1]
			if col == 2 or col == 5:
				nodeDict[node[1]] = node[0]
		return nodeDict

	def get_links(self, nodes=False,database=False,swap=False,only_parents=False,order=False,simple=False):
		'''This function returns all links in the given database'''
		order = ["parent","child"]
		if not simple:
			order += ["rank_i"]
		logger.debug(nodes)
		if swap:
			order[1],order[0] = order[0],order[1]

		if only_parents and nodes:
			QUERY = '''SELECT {order} FROM tree WHERE child in ({nodes})'''.format(nodes=",".join(map(str,nodes)),order=",".join(order))
		elif nodes and order:
			QUERY = '''SELECT {order} FROM tree WHERE child in ({nodes}) ORDER BY parent ASC'''.format(nodes=",".join(map(str,nodes)),order=",".join(order))
		elif nodes:
			QUERY = '''SELECT {order} FROM tree WHERE parent in ({nodes}) OR child in ({nodes})'''.format(nodes=",".join(map(str,nodes)),order=",".join(order))
		else:
			QUERY = '''SELECT {order} FROM tree'''.format(order=",".join(order))
		logger.debug(QUERY)
		if not database:
			database = self.database
		links = self.query(QUERY).fetchall()
		return links

	'''Add functions of class'''
	def add_node(self, description, id=False, table="nodes"):
		'''Add node to tree'''
		if description.strip() == "":  ## for some reason an empty node has been given, do not add empty nodes
			return False
		info = { "name": description }
		### If ID is supplied skip autoincrement and add specific ID
		if id:
			info["id"] = id
		taxid_base = self.insert(info, table=table)
		logger.debug("node added [description {}, taxid base {}]: ".format(info, taxid_base))
		return taxid_base

	def add_rank(self, rank ,id=False):
		'''Add node to tree'''
		if rank == None:
			raise AttributeError("Rank cannot be None")
		info = { "rank": rank }
		### If ID is supplied skip autoincrement and add specific ID
		if id:
			info["id"] = id
		rank_id = self.insert(info, table="rank")
		return rank_id

	def add_link(self, child=None, parent=None, rank=1, table="tree"):
		'''Add relationship in tree'''
		info = {
			"child": child,
			"parent": parent,
			"rank_i": rank
		}
		logger.debug("link added:  child {}, parent {}, rank {} ".format(child,parent,rank))
		return self.insert(info, table="tree")

	def add_genome(self, genome, _id=False,reference=False):
		'''Add genome annotation to nodes'''
		info = {
			"genome": genome
		}

		if _id:
			info["id"] = _id
		if reference:
			info["reference"] = reference
		logger.debug(info)
		return self.insert(info, table="genomes")

	def add_links(self,links, table="tree",hold=False):
		'''Add links from a list to tree'''
		added_links = []
		nodes = set()
		for parent,child,rank in links:
			res = self.add_link(child,parent,rank,table=table)
			### Check if the link already exist in the database, this overlap may occur when a large new branch is added
			if "UNIQUE constraint failed" not in str(res):
				added_links.append([parent,child,rank])
				nodes.add(parent)
				nodes.add(child)
		## Commit changes
		if not hold:
			self.commit()
		return added_links,nodes

	def add_nodes(self,nodes, table="tree",hold=False):
		'''Add nodes from a list of nodes'''
		added_nodes = 0
		for node in nodes:
			res = self.add_node(node)
			### If the node already exist in the database
			if "UNIQUE constraint failed" not in str(res):
				added_nodes +=1
		## Commit changes
		if not hold:
			self.commit()
		return added_nodes

	'''Delete functions of class'''
	def delete_links(self,links, table="tree",hold=False):
		'''This function deletes all links given in links'''
		logger.info("Slow clean")
		QUERY = "DELETE FROM {table} WHERE parent = {parent} AND child = {child} AND rank_i = {rank_i}"
		logger.debug("Deleting {nlinks} links!".format(nlinks=len(links)))
		logger.debug(QUERY.format(table=table,parent="",child="",rank_i=""))
		logger.info("{links}".format(links=len(links)))
		for parent,child,rank in links:
			logger.debug("{}-{} rank: {} deleted!".format(parent,child,rank))
			res = self.query(QUERY.format(table=table, parent=parent, child=child, rank_i=rank))
		## Commit changes
		if not hold:
			logger.debug("Commit changes!")
			self.commit()
		return True

	def fast_delete_links(self,links,table="tree",hold=False):
		'''This function is used when general clean function is executed'''
		logger.info("Fast clean")
		parents = set()
		children = set()
		QUERY = "DELETE FROM {table} WHERE parent in ({parents}) AND child in ({children})"
		logger.debug("Deleting {nlinks} links!".format(nlinks=len(links)))
		for parent,child,rank in links:
			parents |= set([parent])
			children |= set([child])
		logger.debug(QUERY.format(table=table,parents=len(parents),children=len(children)))
		res = self.query(QUERY.format(table=table, parents=",".join(list(map(str,parents))), children=",".join(list(map(str,children)))))
		## Commit changes
		if not hold:
			logger.debug("Commit changes!")
			self.commit()
		return True

	def ambigious_delete_links(self,nodes,hold=False):
		'''This funtion takes as input a list of nodes. It will search the "tree" table for rows that return on "parent in list_of_nodes" and "child in list_of_nodes"'''
		logger.info("Entered ambigious_delete_links")
		QUERY = "DELETE FROM {table} WHERE parent in ({parents}) AND child in ({children})"
		logger.info("Attempting delete where parent&child in: ")
		logger.debug( str(nodes))
		res = self.query(QUERY.format(table='tree', parents=",".join(list(map(str,nodes))), children=",".join(list(map(str,nodes)))))
		## Commit changes
		if not hold:
			logger.debug("Commit changes!")
			self.commit()
		return True

	def delete_nodes(self, nodes, table="nodes",hold=False):
		'''This function deletes all nodes given in nodes'''
		QUERY = "DELETE FROM {table} WHERE id in ({nodes})"
		logger.debug("Deleting {nnodes} nodes!".format(nnodes=len(nodes)))
		logger.debug(QUERY.format(table=table,nodes=""))
		res = self.query(QUERY.format(table=table, nodes=",".join(list(map(str,nodes)))))
		## Commit changes
		if not hold:
			logger.debug("Commit changes!")
			self.commit()
		return True

	def delete_genomes(self,nodes,genomes='',match_genome_only=False,table="genomes",hold=False):
		'''This function clean up genomes from given nodes'''
		if genomes: # check if genomes were specified for targeted deletion (sometimes a node can hold multiple genomes)
			if not match_genome_only:
				QUERY = "DELETE FROM {table} WHERE id in({nodes}) AND genome in({genomes})"
			else:
				QUERY = "DELETE FROM {table} WHERE genome in({genomes})"
		else:
			QUERY = "DELETE FROM {table} WHERE id in({nodes})"

		logger.info("Deleting {nnodes} annotations!".format(nnodes=len(nodes)))
		if genomes:
			if not match_genome_only:
				logger.debug(QUERY.format(table=table,nodes="",genomes=""))
				res = self.query(QUERY.format(table=table, nodes=",".join(list(map(str,nodes))), genomes='"'+'","'.join(list(map(str,genomes)))+'"'))
			else:
				res = self.query(QUERY.format(table=table, genomes='"'+'","'.join(list(map(str,genomes)))+'"'))
		else:
			logger.debug(QUERY.format(table=table,nodes=""))
			res = self.query(QUERY.format(table=table, nodes=",".join(list(map(str,nodes)))))

		## Commit changes
		if not hold:
			logger.debug("Commit changes!")
			self.commit()
		return True

	def num_rows(self,table):
		'''Find out the number of rows in a table'''
		QUERY = '''SELECT Count(*) FROM {table}'''
		logger.debug(QUERY)
		return self.query(QUERY.format(table=table)).fetchall()[0][0]


class ModifyFunctions(DatabaseFunctions):
	"""ModifyFunctions adds nessesary functions when modifying a database"""
	def __init__(self, database, verbose=False):
		super().__init__(database, verbose)
		logger.debug("Load ModifyFunctions")

	def get_rank(self,col=1):
		'''Get rank index from database'''
		QUERY = "SELECT rank_i,rank FROM rank"
		logger.debug(QUERY)
		rankDict = {}
		for rank in self.query(QUERY).fetchall():
			rankDict[rank[0]] = rank[1]
			if col == 1:
				continue
			rankDict[rank[1]] = rank[0]
		return rankDict

	def get_children(self,parents,children=set(),level=0,maxdepth=50,selected=False):
		'''Get all children from a parent'''
		if level >= maxdepth:
			return children
		if not parents:
			return children
		QUERY = '''SELECT child FROM tree WHERE parent in({nodes})'''.format(nodes=",".join(map(str,list(parents))))

		if selected:
			QUERY = '''SELECT child FROM tree WHERE parent in({nodes}) and rank_i = {rank}'''.format(nodes=",".join(map(str,list(parents))), rank=selected)
		logger.debug(QUERY)
		res = self.query(QUERY).fetchall()
		if (len(res) + len(children)) != 0:
			children = set([child_i[0] for child_i in res])
			if level < maxdepth:
				try:
					children |= self.get_children(parents=children,level=level+1,maxdepth=maxdepth)
				except RecursionError:
					logger.debug('Maximum depth reached, halting tree traversal')
					return children
		return children

	def get_parent(self,name,all=False):
		'''Get parent from node id parent'''
		QUERY = '''SELECT parent,child,rank_i FROM tree WHERE child = "{node}"'''.format(node=name)
		logger.debug(QUERY)
		if all:
			res = self.query(QUERY).fetchall()
		else:
			res = self.query(QUERY).fetchone()
		return res

	def parse_parents(self,parents):
		pret = set()
		for pc in parents:
			pret |= set([int(pc[0])])
			pret |= set([int(pc[1])])
		return pret

	def get_parents(self,name,parents=set(),depth=0,find_all=False,simple=False,ncbi=False):
		'''Get all parents until root'''
		ld = depth+1
		lmax = 12
		if ncbi and depth == 0:
			logger.info("NCBI mode on, minimum level to check 50 (since euk has red many levels) Will take time!")
		lmax = 50
		if simple:
			simple=""
		else:
			simple=",rank_i"
		try:
			name = set(list(map(int,name))) ## Make sure all names are int
		except TypeError:
			name = set(list(map(int,[x[0] for x in name]))) ## Make sure all names are int
		if isinstance(name, int):
			name = [name]
		if not name:
			return parents
		if len(name & set([1])) > 0 and not ncbi:
			return parents
		QUERY = '''SELECT parent,child{simple} FROM tree WHERE child in ({node})'''.format(simple=simple, node=",".join(map(str,name)))
		res = self.query(QUERY).fetchall()
		if not find_all:
			try:
				res = int(res[0])
			except IndexError:
				logger.warning("WARNING: parent could not be found for node {res} \n{query}".format(query=QUERY, res=name))
				return set()
			except TypeError:
				logger.error("ERROR: cannot fetch parents when parents are duplicated!")
				raise TreeError("cannot fetch parents when parents are duplicated!")
			if res[0] == res[1]:  ## Special for root
				return set([res[0]])
			else:
				## Add all parents
				try:
					parents |= self.get_parents([int(res[0])],depth=ld,ncbi=ncbi)

				except RecursionError:
					print(res, name)
			## Add current node
			parents |= set([int(res[0])])
		else:
			parents = self.parse_parents(res) #set([(int(pc[0]),int(pc[-1])) for pc in res])
			if len(parents-name) == 1 and ld > lmax:  ## There are no more parents to fetch
				return parents
			else:
				## Add all parents
				if ld > lmax:
					return parents
				try:
					parents |= self.get_parents(parents,find_all=find_all,depth=ld,ncbi=ncbi)
				except RecursionError:
					logger.info("RecursionError, node requested")
					print(res, name)
		return parents

	def get_id(self,name):
		'''get node id from name'''
		QUERY = '''SELECT id FROM nodes WHERE name = ? COLLATE NOCASE'''
		try:
			res = self.query(QUERY,insert_val=(name,),error=True).fetchone()[0]
		except TypeError:
			logger.debug(QUERY)
			raise NameNotFoundError("Name not found in the database! {name}".format(name=name))
		return res

	def get_name(self,id):
		'''get node name from id'''
		QUERY = '''SELECT name FROM nodes WHERE id = ?'''
		try:
			res = self.query(QUERY,insert_val=(id,),error=True).fetchone()[0]
		except TypeError:
			logger.debug(QUERY)
			raise NameNotFoundError("Name not found in the database! {id}".format(id=id))
		return res

	def update_table(self,data,table):
		'''Add annotation to table'''
		return self.update(data, table=table)

	def update_genome(self,data):
		'''Add genome annotation to nodes'''
		return self.update_table(data, table="genomes")

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
