#!/usr/bin/env python3 -c
'''
CreateDatabase is a simple script to simplify the creation of a new database
'''

#**** Imports ****
import sys, os
import sqlite3

import logging
#logger = logging.getLogger(__name__)

#from CanSNPer2.modules.ReadCanSNP_tables import ReadCanSNPer
from CanSNPer2.modules.InitiateDatabase import CanSNPDatabase
from CanSNPer2.modules.ModifyDatabase import ModifyCanSNPer2Database
import flextaxd

def get_read_modules():
	'''Find (ReadTaxonomy) modules and return options'''
	read_modules = []
	for script in os.listdir(os.path.dirname(os.path.abspath(flextaxd.__file__))+"/modules"):
		if script.startswith("ReadTaxonomy"):
			modname = script.lstrip("ReadTaxonomy").rstrip(".py")
			if modname != "":
				read_modules.append(modname)
	return read_modules

supported_input = set(get_read_modules())
supported_output = set(["cansnper","ncbi","tab","newick"])


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
modify_database.add_argument('--remove',     action='store_true',                         	help="If node is given, instead of replace/update remove branch from node")
modify_database.add_argument('--replace',     action='store_true',                          help="replace node")
#modify_database.add_argument('--add_node',     metavar='',                                 	help="Add a single node, parent,node,children")
#modify_database.add_argument('--add_snp',     metavar='',                                 	help="Add a single snp,  parent,node,children")


export_database = parser.add_argument_group("Export tree")
export_database.add_argument('--export_tree', metavar='', choices=supported_output, 		help="Export tree to text format")

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
	os.system("mkdir -p {path}".format(path=args.logs))
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
		'''This should be reprogrammed so that one can add annotations without supplying a mod file!'''
		CanSNPer2Mod_DB = ModifyCanSNPer2Database(**vars(args))
		CanSNPer2Mod_DB.update_database()
		if args.annotation:
			'''Update for snp annotations supplied, run load annotation file function'''
			logger.info("Load new SNP annotations!")
			CanSNPer2Mod_DB.load_annotation_file(args.annotation)
		if args.references:
			logger.info("Load new genomes!")
			CanSNPer2Mod_DB.load_genome_reference_file(args.references)
	elif args.annotation:
		CanSNPDB = CanSNPDatabase(**vars(args))
	else:
		logger.error("No datafile supplied, nothing to process!")
		exit()
