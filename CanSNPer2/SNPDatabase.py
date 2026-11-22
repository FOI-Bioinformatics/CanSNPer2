#!/usr/bin/env python3
'''
CreateDatabase is a simple script to simplify the creation of a new database
'''

#**** Imports ****
import sys, os
import shutil
import sqlite3
import tempfile

import logging

from CanSNPer2.modules.TopLevel import run_main
from CanSNPer2.modules.InitiateDatabase import CanSNPDatabase
from CanSNPer2.modules.ModifyDatabase import ModifyCanSNPer2Database
from CanSNPer2.modules.ModifyTree import InputError, TreeError
from CanSNPer2.modules.DatabaseConnection import check_database, ConnectionError, SourceFileError, DatabaseWriteError
from CanSNPer2.modules.DatabaseConnection import TreeError as DatabaseTreeError
from CanSNPer2.modules.SourceCheck import check_modification
from CanSNPer2.modules.ReadTaxonomyCanSNPer import InputError as TreeFileError
from CanSNPer2.modules.ExportDatabase import export_source, tree_paths, ExportError, write_outputs
from CanSNPer2.modules.ImportDatabase import import_from_sources, DatabaseImportError
from CanSNPer2.modules.NewickTree import NewickTree
from CanSNPer2.modules.SchemeInfo import (STATUSES, MISSING_TEXT, database_counts, read_scheme_info, ordered_keys,
	log_scheme_notice, build_info, modified_info, write_scheme_info)
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
	create_database.add_argument('--annotation',   metavar='', default=False,                	help="CanSNPer canSNP source file (<name>_canSNPs.txt; a 2.3.0 <name>_snp.txt is read as well)")
	create_database.add_argument('--references',   metavar='', default=False,                 	help="File containing all reference genomes listed")
	create_database.add_argument('--source_type',  metavar='', default="CanSNPer",
														choices=supported_input, 				help="Select source file type")
	create_database.add_argument('--create',             action='store_true',                   help="Create new database!")
	create_database.add_argument('--import',     dest='do_import', action='store_true',         help="build a database from the source files written by --export (--tree --annotation --references), written to --database")
	create_database.add_argument('--scheme_name',   metavar='',                                 help="with --create or --import: the name of the scheme, recorded in the database")
	create_database.add_argument('--scheme_status', metavar='',                                 help="with --create or --import: the status of the scheme, verified, draft or unknown, recorded in the database")
	create_database.add_argument('--scheme_source', metavar='',                                 help="with --create or --import: where the scheme files came from (a repository or folder), recorded in the database")

	info_database = parser.add_argument_group("Database information")
	info_database.add_argument('--info',         action='store_true',                           help="print the scheme information the database records (scheme, status, source, CanSNPer2 version) and exit")


	convert_database = parser.add_argument_group("Convert from the format before 2.3.0")
	convert_database.add_argument('--convert',   action='store_true',                           help="bring a SNP source file (--annotation) or a database (--database) of an older format to the current one, written to -o: positions along the joined reference sequences become positions on each reference sequence, and a SNP file becomes a canSNP file")
	convert_database.add_argument('--refdir',    metavar='',                                    help="folder with the reference files <genome>.fna used by --convert; missing files are downloaded")

	modify_database = parser.add_argument_group("Database modifications")
	modify_database.add_argument('--mod_file',     metavar='',                                 	help="File with modifications/update to the tree")
	modify_database.add_argument('--parent',         metavar='',                                 	help="Node (or nodes matching tree file) from which to update/replace/remove")
	modify_database.add_argument('--remove',     action='store_true',                         	help="Remove the whole branch rooted at the --parent node, with its nodes, links and SNP annotations")
	modify_database.add_argument('--replace',     action='store_true',                          help="replace node")


	export_database = parser.add_argument_group("Export database")
	export_database.add_argument('--export',  action='store_true',					help="Export the database to text files")
	export_database.add_argument('--export_format', metavar='', default=None,  			choices=supported_output, help="source: the tree, SNP and reference files that --create reads; newick: the tree; tab is the old name of source (default source), only with --export")
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
		genome_files, table_genomes, database_genomes, references_from_file, references_from_database, check_target,
		check_table_format, check_database_format, canSNPs_target, database_stage, table_stage,
		upgrade_database, upgrade_table)
	if bool(args.annotation) == bool(args.database):
		logger.error("--convert converts either a SNP source file (--annotation) or a database (--database), give one of them!")
		return 1
	if not args.outdir:
		logger.error("--convert writes the converted file to a folder, give it with -o/--outdir!")
		return 1
	source = args.annotation or args.database
	sequences = set()
	'''what the input needs: its positions moved onto each reference sequence (the format before
		2.3.0, which needs the reference files), or only the rank, score and typed columns of the
		canSNP file (2.4.0)'''
	try:
		if args.annotation:
			stage = table_stage(source)
		else:
			try:
				check_database(source, current_format=False)
			except ConnectionError as e:
				raise ConversionError(str(e))
			stage = database_stage(source)
	except (ConversionError, ConnectionError, sqlite3.Error, OSError) as e:
		logger.error(str(e))
		return 1
	if stage == "joined" and not args.refdir:
		logger.error("--convert checks each SNP against its reference, give the folder with the reference files <genome>.fna with --refdir!")
		return 1
	target = canSNPs_target(args.outdir, source) if args.annotation else os.path.join(
		args.outdir, os.path.basename(source))
	'''a converted or upgraded database records that this version wrote it; the scheme and its
		status are not known from the old database'''
	scheme_info = build_info(status="unknown", source="converted from {name}".format(name=os.path.basename(source)))
	try:
		check_target(source, target)
		if args.annotation:
			check_table_format(source) ## before any reference is fetched
			if stage == "joined":
				references = references_from_file(args.references) if args.references else None
				files = genome_files(args.refdir, table_genomes(source), references)
				joined = os.path.join(args.outdir, "." + os.path.basename(source) + ".positions")
				count = convert_table(source, joined, files, sequences=sequences)
				try:
					upgrade_table(joined, target)
				finally:
					if os.path.exists(joined):
						os.remove(joined)
			else:
				count = upgrade_table(source, target)
		else:
			check_database_format(source) ## before any reference is fetched
			if stage == "joined":
				references = references_from_database(source)
				if args.references:
					references.update(references_from_file(args.references))
				files = genome_files(args.refdir, database_genomes(source), references)
				count = convert_database(source, target, files, sequences=sequences, scheme_info=scheme_info)
			else:
				count = upgrade_database(source, target, scheme_info=scheme_info)
	except (ConversionError, ConnectionError, sqlite3.Error, OSError) as e:
		'''every failure of the conversion is one error message and exit code 1, not a traceback'''
		logger.error(str(e))
		return 1
	snps = "canSNP" if count == 1 else "canSNPs"
	if stage == "joined":
		print("converted {n} {snps} of {name} to positions on {m} reference {noun}, written to {target}".format(
			n=count, snps=snps, name=os.path.basename(source), m=len(sequences),
			noun="sequence" if len(sequences) == 1 else "sequences", target=target))
	else:
		print("upgraded {n} {snps} of {name} to the canSNP format, written to {target}".format(
			n=count, snps=snps, name=os.path.basename(source), target=target))
	return 0


def check_mode(args, logger):
	'''The modes are --create, --import, --export, --info, --convert and a modification
		(--mod_file, --remove, --replace); the source files (--tree, --annotation, --references)
		and the options of a mode are not modes. More than one mode used to be resolved by
		precedence, one of them silently dropped. Returns False after an error message.'''
	modes = []
	for flag, given in (("--create", args.create), ("--import", args.do_import), ("--export", args.export),
			("--info", args.info), ("--convert", args.convert)):
		if given:
			modes.append(flag)
	modification = [flag for flag, given in (("--mod_file", args.mod_file), ("--remove", args.remove),
		("--replace", args.replace)) if given]
	if modification:
		modes.append(" ".join(modification))
	if len(modes) > 1:
		logger.error("{given} cannot be combined, give one of --create, --import, --export, --info, --convert "
			"or a modification (--mod_file, --remove) at a time".format(given=" and ".join(modes)))
		return False
	if args.export_format is not None and not args.export:
		logger.error("--export_format {format} only applies with --export, add --export or remove --export_format".format(
			format=args.export_format))
		return False
	if args.replace and not args.mod_file:
		logger.error("--replace replaces a node with the tree in --mod_file, give --mod_file (and --parent) with --replace!")
		return False
	return True


def check_folder(database, logger):
	'''The folder of a database that is to be written must exist. Returns False after an error message.'''
	folder = os.path.dirname(database) or "."
	if not os.path.isdir(folder):
		logger.error("The folder {folder} for the database {db} does not exist, make the folder or choose "
			"another path".format(folder=folder, db=database))
		return False
	return True


def modify(args, logger):
	'''--mod_file, --annotation, --references (each alone or together) and --remove change an existing database.
		The loaders used to skip a SNP on an unknown node or genome and a short row with a warning,
		and to ignore a failed insert, and the command ended with exit status 0 and a database
		that had changed in part. Now the changes are made on a copy beside the database: the
		tree is changed first, the SNP and references files are then checked against the tree
		and the reference genomes of the copy (SourceCheck.check_modification), and the copy
		replaces the database only when every step succeeded and it passes check_database.
		A run that fails leaves the database byte-identical. A run that succeeds marks the
		scheme_info of the database as modified (SchemeInfo.modified_info), so a changed
		database no longer claims to be the published scheme. A database given as a symbolic link
		is resolved first, so the file it points to is replaced and the link is kept. Returns
		the exit code.'''
	database = args.database
	target = os.path.realpath(database)
	work = None
	modifier = None
	try:
		fd, work = tempfile.mkstemp(prefix=".{name}.".format(name=os.path.basename(target)), suffix=".part",
			dir=os.path.dirname(target))
		os.close(fd)
		shutil.copyfile(target, work)
		shutil.copymode(target, work) ## the database keeps its mode, not the 0600 of mkstemp
		options = dict(vars(args))
		options["database"] = work
		modifier = ModifyCanSNPer2Database(**options)
		if args.remove:
			modifier.remove_branch(args.parent)
		else:
			modifier.update_database()
			problems = check_modification(work, annotation=args.annotation or None,
				references=args.references or None, name=database)
			if problems:
				raise SourceFileError("The files cannot be added to {db}, {n} problem(s):\n  {list}".format(
					db=database, n=len(problems), list="\n  ".join(problems)))
			if args.references:
				'''the references first, so the new SNPs can be placed on new genomes'''
				logger.info("Load new genomes!")
				modifier.load_genome_reference_file(args.references)
			if args.annotation:
				logger.info("Load new SNP annotations!")
				modifier.load_annotation_file(args.annotation)
		modifier.close()
		modifier = None
		check_database(work)
		'''the modified tree must still be a tree: one root, one parent per node and every node
			reachable from the root. Links that gave a node a second parent or made a cycle were
			written, and the database then could not be exported'''
		conn = sqlite3.connect(work)
		try:
			tree_paths(conn)
		except ExportError as e:
			raise TreeError("the modification would leave {db} with a tree that is not valid: {e}".format(db=database, e=e))
		finally:
			conn.close()
		'''the changed database no longer matches the scheme files it was built from, so its
			scheme_info is marked: status unknown, commit +modified, the date and the version'''
		conn = sqlite3.connect(work)
		try:
			write_scheme_info(conn, modified_info(read_scheme_info(work)))
		finally:
			conn.close()
		os.replace(work, target)
	except (InputError, TreeError, DatabaseTreeError, SourceFileError, DatabaseWriteError, ConnectionError,
			sqlite3.Error, OSError) as e:
		'''messages name the database, not the copy that was worked on'''
		message = getattr(e, "message", None) or str(e)
		if work is not None:
			message = message.replace(work, database)
		logger.error(message)
		logger.error("{db} was not changed".format(db=database))
		return 1
	finally:
		if modifier is not None:
			modifier.close()
		if work is not None and os.path.exists(work):
			os.remove(work)
	return 0


def command():
	'''Modify the CanSNPer2 database, if it doesn´t exist (create is added create database from files)'''
	args = get_arguments()
	logger = get_logger(args)
	logger.debug(args)
	if not check_mode(args, logger):
		exit(1)
	args.export_format = args.export_format or "source"
	scheme_options = [option for option, value in (("--scheme_name", args.scheme_name),
		("--scheme_status", args.scheme_status), ("--scheme_source", args.scheme_source)) if value is not None]
	if scheme_options and not (args.create or args.do_import):
		logger.error("{options} {verb} recorded by --create and --import only, they cannot be "
			"used with an existing database".format(options=", ".join(scheme_options),
			verb="is" if len(scheme_options) == 1 else "are"))
		exit(1)
	for option, value in (("--scheme_name", args.scheme_name), ("--scheme_source", args.scheme_source)):
		if value is not None and any(c in value for c in "\t\n\r"):
			logger.error("{option} contains a control character (tab, newline or carriage return), "
				"remove it; the value is recorded in the database and shown by --info".format(option=option))
			exit(1)
	if args.scheme_status is not None and args.scheme_status not in STATUSES:
		logger.error("--scheme_status {status} is not a known status, give verified, draft or unknown".format(
			status=args.scheme_status))
		exit(1)
	scheme_info = {"scheme": args.scheme_name, "status": args.scheme_status, "source": args.scheme_source}
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
		if not check_folder(args.database, logger):
			exit(1)
		try:
			counts = import_from_sources(args.tree, args.annotation, args.references,
				args.database, tmpdir=args.tmpdir, verbose=args.verbose, scheme_info=scheme_info)
		except DatabaseImportError as e:
			logger.error(str(e))
			exit(1)
		print("Imported {tree}, {snp}, {references} into {dst} (current format)".format(
			tree=args.tree, snp=args.annotation, references=args.references, dst=args.database))
		print("{nodes} nodes, {tree} links, {snp_annotation} canSNPs ({typed} typed), "
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
		if args.info:
			info = read_scheme_info(args.database)
			if not info:
				print(MISSING_TEXT)
			for key in ordered_keys(info):
				## an empty value, such as known_issues of a scheme without any, is printed as none
				print("{key}\t{value}".format(key=key, value=info[key] if info[key] not in (None, "") else "none"))
			counts = database_counts(args.database)
			if counts:
				print("")
				for key, value in counts:
					print("{key}\t{value}".format(key=key, value=value))
			exit()
	elif args.info:
		logger.error("--info prints what an existing database records, it cannot be combined with --create!")
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
		exit(modify(args, logger))
	if args.export: ## Dump database to file
		outdir = args.outdir or "." ## without -o the export goes to the current directory (issue #132)
		log_scheme_notice(args.database, logger) ## the scheme_info table itself is not exported
		if args.export_format == "newick":
			logger.info("Export Newick tree!")
			NT = NewickTree(args.database, outdir=outdir)
			if args.outdir:
				'''-o is documented as the export folder, so write the tree to a file instead of only
					printing it to stdout when it is given (issue #131)'''
				newick_file = os.path.join(args.outdir,"{db}_tree.nwk".format(db=os.path.basename(args.database).rsplit(".",1)[0]))
				try:
					write_outputs({newick_file: "{tree}\n".format(tree=NT.newickTree)}) ## never overwrites, no partial file
				except ExportError as e:
					logger.error("The database cannot be exported: {error}".format(error=e))
					exit(1)
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
			print("Import back with: CanSNPer2-database --import --database <new.db> --tree {tree} "
				"--annotation {snp} --references {references}".format(**files))
		exit()
	if not args.create:
		if not (args.mod_file or args.annotation or args.references):
			logger.error("Nothing to do with {db}: give a mode, --export, --info or --remove --parent, or the changes to "
				"make, --mod_file, --annotation or --references!".format(db=args.database))
			exit(1)
		'''--annotation and --references need no --mod_file. A modification without any source is a usage error, not a crash (issue #150); a source
			file without the expected columns is refused, not a crash in the loader (issue #154)'''
		exit(modify(args, logger))
	elif args.annotation or args.tree or args.references:
		'''a database can also be created from only a tree or only a reference list (issue #129)'''
		if not args.database:
			logger.error("--create needs the new database to write, give it with --database!")
			exit(1)
		if os.path.exists(args.database):
			'''--create into an existing file added the tree, the references and the SNPs a second
				time; into a database from before 2.3.0 every SNP insert also failed and the
				command still ended with exit status 0. An existing file is refused, as --import does,
				and a database from before 2.3.0 is named as such with the conversion hint.'''
			try:
				check_database(args.database)
			except ConnectionError as e:
				logger.error(str(e))
				exit(1)
			logger.error("The database {db} already exists, --create builds a new database, "
				"choose a new name or remove it first".format(db=args.database))
			exit(1)
		if not check_folder(args.database, logger):
			exit(1)
		try:
			CanSNPDB = CanSNPDatabase(scheme_info=scheme_info, **vars(args))
		except (InputError,TreeFileError,SourceFileError,DatabaseWriteError,sqlite3.Error,OSError) as e:
			'''a tree, node, link, reference or SNP that could not be written ends the build'''
			logger.error(str(e))
			if os.path.exists(args.database):
				'''the file did not exist before this run, so a database that was started is incomplete'''
				os.remove(args.database)
			exit(1)
		except Exception:
			'''an unexpected error is reported by the top-level handler; the started database
				is removed first, as for the known errors above'''
			if os.path.exists(args.database):
				os.remove(args.database)
			raise
	else:
		logger.error("No datafile supplied, nothing to process!")
		exit(1)


def main():
	'''the entry point of the command: command() with the handler of unexpected errors'''
	return run_main(command)


if __name__ == "__main__":
	'''without this guard "python -m CanSNPer2.SNPDatabase" parsed its arguments and then did nothing (issue #128)'''
	main()
