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
from CanSNPer2.modules.ExportDatabase import export_source, ExportError
from CanSNPer2.modules.ImportDatabase import import_from_sources, DatabaseImportError
from CanSNPer2.modules.NewickTree import NewickTree
from CanSNPer2.modules.LogSetup import setup_logging
from CanSNPer2 import __version__, COMMANDS_HELP  ## __version__ used by the --version flag (issue #119)
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
## "source" writes the tree, SNP and reference files that --create reads; "tab" is its old name,
## which wrote node and link tables that could not rebuild the database
supported_output = ["source", "newick", "tab"]


import argparse


def get_arguments():
	'''The command used to parse its arguments and open its log file at import time, so
		importing CanSNPer2.SNPDatabase consumed sys.argv and wrote a log file as a side
		effect (issue #174)'''
	parser = argparse.ArgumentParser(
		description='CanSNPer2-database - create, modify and export CanSNP databases',
		epilog=COMMANDS_HELP,
		formatter_class=argparse.RawDescriptionHelpFormatter)

	required_args = parser.add_argument_group('Required')
	required_args.add_argument('-db',  '--database',     metavar='',                            help='CanSNPer2 database name')


	create_database = parser.add_argument_group("Load database")
	create_database.add_argument('--tree',         metavar='', default=False,                	help="CanSNPer tree source file")
	create_database.add_argument('--annotation',   metavar='', default=False,                	help="CanSNPer snp source file")
	create_database.add_argument('--references',   metavar='', default=False,                 	help="File containing all reference genomes listed")
	create_database.add_argument('--source_type',  metavar='', default="CanSNPer",
														choices=supported_input, 				help="Select source file type")
	create_database.add_argument('--create',             action='store_true',                   help="Create new database!")
	create_database.add_argument('--import',     dest='do_import', action='store_true',         help="build a database from the source files written by --export (--tree --annotation --references), written to --database")


	convert_database = parser.add_argument_group("Convert from the format before 2.3.0")
	convert_database.add_argument('--convert',   action='store_true',                           help="convert a SNP source file (--annotation) or a database (--database) from positions along the joined reference sequences to positions on each reference sequence, written to -o")
	convert_database.add_argument('--refdir',    metavar='',                                    help="folder with the reference files <genome>.fna used by --convert; missing files are downloaded")

	modify_database = parser.add_argument_group("Database modifications")
	modify_database.add_argument('--mod_file',     metavar='',                                 	help="File with modifications/update to the tree")
	modify_database.add_argument('--parent',         metavar='',                                 	help="Node (or nodes matching tree file) from which to update/replace/remove")
	modify_database.add_argument('--remove',     action='store_true',                         	help="Remove the whole branch rooted at the --parent node, with its nodes, links and SNP annotations")
	modify_database.add_argument('--replace',     action='store_true',                          help="replace node")


	export_database = parser.add_argument_group("Export database")
	export_database.add_argument('--export',  action='store_true',					help="Export the database to text files")
	export_database.add_argument('--export_format', metavar='', default='source',  			choices=supported_output, help="source: the tree, SNP and reference files that --create reads; newick: the tree; tab is the old name of source (default source)")
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
	return args


def get_logger(args):
	'''Set the log level and open the log file of this command'''
	logval = args.supress
	if args.debug:
		logval = args.debug
	elif args.verbose:
		logval = args.verbose

	args.logs = args.logs.rstrip("/")+"/"
	logpath = setup_logging(args.logs, "CanSNPer2-database", level=logval) ## shared log naming for all commands (issue #217)
	return logging.getLogger(__name__)


def convert(args, logger):
	'''--convert: a SNP source file or a database with positions along the joined reference
		sequences (the format before 2.3.0) is converted to positions on each reference sequence.
		It is handled before the database checks, which refuse the old format. Returns the exit code.'''
	from CanSNPer2.modules.ConvertPositions import (ConversionError, convert_table, convert_database,
		genome_files, table_genomes, database_genomes, references_from_file, references_from_database, check_target)
	if bool(args.annotation) == bool(args.database):
		logger.error("--convert converts either a SNP source file (--annotation) or a database (--database), give one of them!")
		return 1
	if not args.outdir:
		logger.error("--convert writes the converted file to a folder, give it with -o/--outdir!")
		return 1
	if not args.refdir:
		logger.error("--convert checks each SNP against its reference, give the folder with the reference files <genome>.fna with --refdir!")
		return 1
	source = args.annotation or args.database
	target = os.path.join(args.outdir, os.path.basename(source))
	sequences = set()
	try:
		check_target(source, target)
		if args.annotation:
			references = references_from_file(args.references) if args.references else None
			files = genome_files(args.refdir, table_genomes(source), references)
			count = convert_table(source, target, files, sequences=sequences)
		else:
			try:
				check_database(source, current_format=False)
			except ConnectionError as e:
				raise ConversionError(str(e))
			references = references_from_database(source)
			if args.references:
				references.update(references_from_file(args.references))
			files = genome_files(args.refdir, database_genomes(source), references)
			count = convert_database(source, target, files, sequences=sequences)
	except (ConversionError, ConnectionError, sqlite3.Error, OSError) as e:
		'''every failure of the conversion is one error message and exit code 1, not a traceback'''
		logger.error(str(e))
		return 1
	print("converted {n} SNPs of {name} to positions on {m} reference sequences, written to {target}".format(
		n=count, name=os.path.basename(source), m=len(sequences), target=target))
	return 0


def main():
	'''Modify the CanSNPer2 database, if it doesn´t exist (create is added create database from files)'''
	args = get_arguments()
	logger = get_logger(args)
	logger.debug(args)
	if args.convert:
		exit(convert(args, logger))
	if args.do_import:
		'''import is the counterpart of --export: it reads the source triple and writes a
			current-format database, so it is handled before the database checks below,
			which expect --database to be an existing database or a --create target'''
		missing = [name for name, value in (("--tree", args.tree), ("--annotation", args.annotation),
			("--references", args.references)) if not value]
		if missing:
			logger.error("--import builds a database from the files --export writes, it needs all "
				"three source files, missing: {0}!".format(", ".join(missing)))
			exit(1)
		if not args.database:
			logger.error("--import needs the new database to write, give it with --database!")
			exit(1)
		try:
			counts = import_from_sources(args.tree, args.annotation, args.references,
				args.database, tmpdir=args.tmpdir, verbose=args.verbose)
		except DatabaseImportError as e:
			logger.error(str(e))
			exit(1)
		print("Imported {tree}, {snp}, {references} into {dst} (current format)".format(
			tree=args.tree, snp=args.annotation, references=args.references, dst=args.database))
		print("{nodes} nodes, {tree} links, {snp_annotation} SNPs, "
			"{snp_references} reference genomes".format(**counts))
		exit()
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
			try:
				files = export_source(args.database, outdir)
			except ExportError as e:
				logger.error("The database cannot be exported: {error}".format(error=e))
				exit(1)
			print("Exported to {files}".format(files=", ".join(files.values())))
			print("Import back with: CanSNPer2-database --import --database <new.db> --tree {tree} --annotation {snp} --references {references}".format(**files))
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
