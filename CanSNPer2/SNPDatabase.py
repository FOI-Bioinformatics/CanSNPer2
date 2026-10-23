#!/usr/bin/env python3
'''
CreateDatabase is a simple script to simplify the creation of a new database
'''

#**** Imports ****
import sys, os
import sqlite3

import logging

from CanSNPer2.modules.InitiateDatabase import CanSNPDatabase
from CanSNPer2.modules.ModifyDatabase import ModifyCanSNPer2Database
from CanSNPer2.modules.ModifyTree import InputError
from CanSNPer2.modules.DatabaseConnection import check_database, ConnectionError, SourceFileError
from CanSNPer2.modules.WriteTaxonomy import WriteTaxonomy
from CanSNPer2.modules.NewickTree import NewickTree
from CanSNPer2 import __version__  ## used by the --version flag (issue #119)
import CanSNPer2.modules

def get_read_modules():
	'''Find (ReadTaxonomy) modules and return options'''
	read_modules = []
	modulepath = os.path.dirname(os.path.abspath(CanSNPer2.modules.__file__))
	for script in os.listdir(modulepath):
		if script.startswith("ReadTaxonomy") and script.endswith(".py"):
			modname = script[len("ReadTaxonomy"):-len(".py")]
			if modname != "":
				read_modules.append(modname)
	return read_modules

supported_input = set(get_read_modules())
supported_output = set(["tab", "newick"])# set(["cansnper","ncbi","tab","newick"])


import argparse
parser = argparse.ArgumentParser(description='CanSNPer2-database ')

required_args = parser.add_argument_group('Required')
required_args.add_argument('-db',  '--database',     metavar='',                            help='CanSNPer2 database name')


create_database = parser.add_argument_group("Load database")
create_database.add_argument('--tree',         metavar='', default=False,                	help="CanSNPer tree source file")
create_database.add_argument('--annotation',   metavar='', default=False,                	help="CanSNPer snp source file")
create_database.add_argument('--references',   metavar='', default=False,                 	help="File containing all reference genomes listed")
create_database.add_argument('--source_type',  metavar='', default="CanSNPer",
													choices=supported_input, 				help="Select source file type")
create_database.add_argument('--create',             action='store_true',                   help="Create new database!")


modify_database = parser.add_argument_group("Database modifications")
modify_database.add_argument('--mod_file',     metavar='',                                 	help="File with modifications/update to the tree")
modify_database.add_argument('--parent',         metavar='',                                 	help="Node (or nodes matching tree file) from which to update/replace/remove")
modify_database.add_argument('--remove',     action='store_true',                         	help="Remove the whole branch rooted at the --parent node, with its nodes, links and SNP annotations")
modify_database.add_argument('--replace',     action='store_true',                          help="replace node")


export_database = parser.add_argument_group("Export database")
export_database.add_argument('--export',  action='store_true',					help="Export database to text format (exports tree and annotation file)")
export_database.add_argument('--export_format', metavar='', default='tab',  				choices=supported_output, help="Select output format [{format}]".format(format=", ".join(supported_output)))
export_database.add_argument('-o', '--outdir', metavar='',									help="outdir for database export!")


debugopts = parser.add_argument_group("Logging and debug options")
debugopts.add_argument('--tmpdir', 			metavar='', default="/tmp",						help="Specify tmp directory default (/tmp)")
debugopts.add_argument('--logs', 				metavar='', default="logs", 				help="Specify log directory")
debugopts.add_argument('--verbose',			action='store_const', const=logging.INFO,		help="print process info, default no output")
debugopts.add_argument('--debug',				action='store_const', const=logging.DEBUG,	help="print debug info")
debugopts.add_argument('--supress',				action='store_const', const=logging.ERROR,
																default=logging.WARNING,	help="supress warnings")

parser.add_argument("--version", action='store_true', help=argparse.SUPPRESS)

args = parser.parse_args()
if len(sys.argv)==1:
	parser.print_help()
	parser.exit()


if args.version:
	print("CanSNPer2 - version {version}".format(version=__version__))
	exit()


logval = args.supress
if args.debug:
	logval = args.debug
elif args.verbose:
	logval = args.verbose

from datetime import date,time
t = time()
today = date.today()
args.logs = args.logs.rstrip("/")+"/"
logpath = args.logs+"CanSNPer2-database"+today.strftime("%b-%d-%Y")+"{}.log"
if not os.path.exists(args.logs):
	os.makedirs(args.logs, exist_ok=True)
if os.path.exists(logpath):
	logpath=logpath.format("-{:%H:%M}".format(t))
else: logpath = logpath.format("")

logging.basicConfig(
		#filename=logpath,
		level=logval,
		format="%(asctime)s %(module)s [%(levelname)-5.5s]  %(message)s",
	    handlers=[
	        logging.FileHandler(logpath),
	        logging.StreamHandler()
	    ])
logger = logging.getLogger(__name__)

def main():
	'''Modify the CanSNPer2 database, if it doesn´t exist (create is added create database from files)'''
	logger.debug(args)
	if not args.create:
		'''a wrong path let sqlite create an empty database and the modification or export only failed
			later with "no such table" errors (issue #125)'''
		if not args.database:
			logger.error("A database is required, use --database, or --create together with a source file!")
			exit(1)
		try:
			check_database(args.database)
		except ConnectionError as e:
			logger.error(str(e))
			exit(1)
	if args.remove:
		'''--remove was documented but no code ever read the flag, so it silently did nothing (issue #147)'''
		if args.create:
			logger.error("--remove removes a branch from an existing database and cannot be combined with --create!")
			exit(1)
		if args.mod_file:
			logger.error("--remove only needs the branch root given with --parent, it does not read --mod_file!")
			exit(1)
		if not args.parent:
			logger.error("--remove needs the branch to remove, give the branch root with --parent!")
			exit(1)
		try:
			ModifyCanSNPer2Database(**vars(args)).remove_branch(args.parent)
		except InputError as e:
			logger.error(str(e))
			exit(1)
		exit()
	if args.export: ## Dump database to file
		outdir = args.outdir or "." ## without -o the export goes to the current directory (issue #132)
		if args.export_format == "newick":
			logger.info("Export Newick tree!")
			NT = NewickTree(args.database, outdir=outdir)
			if args.outdir:
				'''-o is documented as the export folder, so write the tree to a file instead of only
					printing it to stdout when it is given (issue #131)'''
				if not os.path.exists(args.outdir):
					os.makedirs(args.outdir, exist_ok=True)
				newick_file = os.path.join(args.outdir,"{db}_tree.nwk".format(db=os.path.basename(args.database).rsplit(".",1)[0]))
				with open(newick_file,"w") as newick_out:
					print(NT.newickTree,file=newick_out)
				logger.info("Newick tree written to {file}".format(file=newick_file))
			else:
				NT.print()
		else:
			logger.info("Export CanSNPer2 database")
			WriteObj = WriteTaxonomy( outdir,database=args.database, minimal=True, separator="\t", prefix="annotation,snptree",desc=True)
			WriteObj.set_order(True)
			WriteObj.names()
			WriteObj.nodes()
		exit()
	if not args.create:
		'''This should be reprogrammed so that one can add annotations without supplying a mod file!'''
		try:
			CanSNPer2Mod_DB = ModifyCanSNPer2Database(**vars(args))
			CanSNPer2Mod_DB.update_database()
		except InputError as e:
			'''a modification without a mod_file or mod_database is a usage error, not a crash (issue #150)'''
			logger.error(str(e))
			exit(1)
		if args.annotation:
			'''Update for snp annotations supplied, run load annotation file function'''
			logger.info("Load new SNP annotations!")
			try:
				CanSNPer2Mod_DB.load_annotation_file(args.annotation)
			except SourceFileError as e:
				'''a source file without the expected columns would otherwise crash the loader (issue #154)'''
				logger.error(str(e))
				exit(1)
		if args.references:
			logger.info("Load new genomes!")
			try:
				CanSNPer2Mod_DB.load_genome_reference_file(args.references)
			except SourceFileError as e:
				logger.error(str(e))
				exit(1)
	elif args.annotation or args.tree or args.references:
		'''a database can also be created from only a tree or only a reference list (issue #129)'''
		try:
			CanSNPDB = CanSNPDatabase(**vars(args))
		except (InputError,SourceFileError) as e:
			logger.error(str(e))
			exit(1)
	else:
		logger.error("No datafile supplied, nothing to process!")
		exit(1)

if __name__ == "__main__":
	'''without this guard "python -m CanSNPer2.SNPDatabase" parsed its arguments and then did nothing (issue #128)'''
	main()
