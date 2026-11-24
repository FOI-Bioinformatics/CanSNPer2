#!/usr/bin/env python3
'''
cansnper2 - one command with a verb for every job (2.6.0).

This is the interface described in the CLI audit: cansnper2 type, selftest,
db create|import|export|info|convert|modify|remove, scheme list|get, ref get.
It calls the same run functions as the four old commands (CanSNPer2,
CanSNPer2-database, CanSNPer2-download, CanSNPer2-test), which keep working
with their old arguments; the new names and hyphenated spellings live here
only. The option names after every dest= are the internal names of the
matching old parser, so the shared run sees exactly what it always saw.
'''

import argparse
import logging

from CanSNPer2 import __version__
from CanSNPer2.modules.TopLevel import run_main
from CanSNPer2.modules.Aligners import ALIGNERS, DEFAULT_ALIGNER

import CanSNPer2.CanSNPerTree as typing_command
import CanSNPer2.SNPDatabase as database_command
import CanSNPer2.DownloadGenomes as download_command
import CanSNPer2.selftest as selftest_command


def add_logging_options(parser, logs_dest, logs_default):
	'''-v/--debug/--quiet and the log folder, the same three levels the old commands
		have (--supress becomes --quiet); logs_dest and logs_default are the internal
		name and default of the log folder of this command'''
	group = parser.add_argument_group("Logging options")
	group.add_argument('--log-dir', dest=logs_dest, metavar='DIR',
		default=logs_default, help="folder for the log files")
	group.add_argument('-v', '--verbose', action='store_const', const=logging.INFO,
		help="print process information")
	group.add_argument('--debug', action='store_const', const=logging.DEBUG,
		help="print debug information, and the traceback of an unexpected error")
	group.add_argument('--quiet', dest='supress', action='store_const',
		const=logging.ERROR, default=logging.WARNING, help="do not print warnings")


def add_scheme_source_options(parser):
	'''where schemes are read from: CanSNPer2-data on GitHub, another repository,
		or a local folder laid out the same way'''
	group = parser.add_argument_group("Scheme source options")
	group.add_argument('--repo', metavar='NAME',
		help="GitHub repository laid out as CanSNPer2-data, owner/name "
			"(default FOI-Bioinformatics/CanSNPer2-data)")
	group.add_argument('--branch', metavar='REF',
		help="branch, tag or commit of the repository; in CanSNPer2-data master "
			"holds verified schemes, dev also draft schemes")
	group.add_argument('--data-dir', dest='data_dir', metavar='DIR',
		help="local directory laid out as CanSNPer2-data to read schemes from "
			"instead of a repository")


def add_type_subparser(subparsers):
	parser = subparsers.add_parser('type',
		help="type query assemblies against a CanSNP database",
		description='type query assemblies against a CanSNP database',
		formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	required = parser.add_argument_group("Required arguments")
	required.add_argument('query', nargs=argparse.ZERO_OR_MORE, metavar='query',
		help="file(s) to align (fasta), or one list file with a file path per line")
	required.add_argument('-d', '--database', metavar='DB',
		help="CanSNP database")

	output = parser.add_argument_group("Output options")
	output.add_argument('-o', '--outdir', metavar='DIR',
		help="output directory; the results overview summary.html is written "
			"there, the per-sample files in its samples folder")
	output.add_argument('--save-tree', dest='save_tree', action='store_true',
		help="also write a standalone tree page per sample "
			"(samples/<sample>_tree.html, the overview already shows every tree)")
	output.add_argument('--no-snp-files', dest='no_snpfiles',
		action='store_false',
		help="do not write the per-sample snp files, the call is still printed")
	output.add_argument('--summary', action='store_true',
		help="output a summary file and tree with all called SNPs")
	output.add_argument('--node-labels', dest='annotations', metavar='FILE',
		help="TSV file annotating tree nodes: <snp name>\\t<node|clade>\\t<rank>"
			"\\t<label>[\\t<color>]")

	run = parser.add_argument_group("Run options")
	run.add_argument('--ref-dir', dest='refdir', metavar='DIR',
		help="folder with the reference genomes of the database")
	run.add_argument('--allow-missing-references', dest='allow_missing_references',
		action='store_true',
		help="type even when reference genomes have no file in --ref-dir")
	run.add_argument('--workdir', metavar='DIR',
		help="working directory the other folder options are relative to")
	run.add_argument('--min-required-hits', dest='min_required_hits',
		type=typing_command.positive_int,
		help="minimum number of derived SNPs on the path to the final SNP "
			"for a call")
	run.add_argument('--strictness', type=typing_command.fraction,
		help="fraction of the SNPs on the path that must be called (0-1)")
	run.add_argument('--mask', type=typing_command.non_negative_int,
		help="ignore SNPs closer than this many bases to an alignment edge")
	run.add_argument('--keep-going', dest='keep_going', action='store_true',
		help="continue with the other samples when one fails")
	run.add_argument('--rerun', action='store_true',
		help="type samples again although their result file exists")
	run.add_argument('--aligner', choices=sorted(ALIGNERS),
		help="program that aligns the query against each reference")
	run.add_argument('--alignment-timeout', dest='alignment_timeout',
		type=typing_command.minutes, metavar='MINUTES',
		help="stop the alignments of a query that take longer than this "
			"(0 for no limit)")
	run.add_argument('--aligner-path', dest='aligner_path', metavar='PATH',
		help="the aligner program, or the folder that holds it "
			"(default: found in PATH)")
	run.add_argument('--skip-alignment', dest='skip_alignment',
		action='store_true',
		help="do not align, read the alignments already in --tmp-dir")
	run.add_argument('--keep-temp', dest='keep_temp', action='store_true',
		help="keep the temporary alignment files")
	run.add_argument('--tmp-dir', dest='tmpdir', metavar='DIR', default=None,
		help="folder for temporary alignments (default: CanSNPer2 in the "
			"system temporary folder, see TMPDIR)")

	add_logging_options(parser, logs_dest='logdir', logs_default='logs/')
	'''the internal options the new interface does not offer keep their old
		defaults: --read_input and the hidden --no_save_tree were dropped,
		--version is a global option of cansnper2'''
	parser.set_defaults(read_input=False, no_save_tree=False, version=False,
		query=[], database=None, outdir="results", save_tree=False,
		no_snpfiles=True, summary=False, annotations=None,
		refdir="references/", allow_missing_references=False, workdir="./",
		min_required_hits=3, strictness=0.7, mask=0, keep_going=False,
		rerun=False, aligner=DEFAULT_ALIGNER, alignment_timeout=60,
		aligner_path="", skip_alignment=False, keep_temp=False)
	return parser


def database_defaults(parser, **mode):
	'''every internal option of the database command with its default, and the
		one mode of this verb set to True'''
	settings = dict(database=None, tree=False, annotation=False,
		references=False, source_type="CanSNPer", create=False,
		do_import=False, scheme_name=None, scheme_status=None,
		scheme_source=None, info=False, convert=False, refdir=None,
		mod_file=None, parent=None, remove=False, replace=False,
		export=False, export_format=None, outdir=None, tmpdir="/tmp",
		version=False)
	settings.update(mode)
	parser.set_defaults(**settings)


def add_scheme_info_options(parser):
	'''the scheme identity recorded in a new database'''
	group = parser.add_argument_group("Scheme information (recorded in the database)")
	group.add_argument('--scheme-name', dest='scheme_name', metavar='NAME')
	group.add_argument('--scheme-status', dest='scheme_status', metavar='STATUS',
		help="verified, draft or unknown")
	group.add_argument('--scheme-source', dest='scheme_source', metavar='SRC',
		help="where the scheme files came from")


def add_source_file_options(parser):
	'''the three source files of db create and db import; --annotation becomes
		--cansnps because the 2.4.0 file it reads lists canSNPs'''
	parser.add_argument('--tree', metavar='FILE',
		help="tree source file (written by db export)")
	parser.add_argument('--cansnps', dest='annotation', metavar='FILE',
		help="canSNP source file (<name>_canSNPs.txt)")
	parser.add_argument('--references', metavar='FILE',
		help="file listing all reference genomes")


def add_db_subparsers(subparsers):
	db = subparsers.add_parser('db',
		description='create, import, export, inspect and modify CanSNP databases',
		help="create, import, export, info, convert, modify, remove")
	db_sub = db.add_subparsers(dest='db_action', metavar='<action>', required=True)
	parsers = {}

	create = db_sub.add_parser('create', help="create a database from source files",
		formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	create.add_argument('-d', '--database', metavar='DB',
		help="the new database to write")
	add_source_file_options(create)
	add_scheme_info_options(create)
	add_logging_options(create, logs_dest='logs', logs_default='logs')
	parsers['db create'] = create

	imp = db_sub.add_parser('import',
		help="build a database from the source files written by db export",
		formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	imp.add_argument('-d', '--database', metavar='DB',
		help="the new database to write")
	add_source_file_options(imp)
	add_scheme_info_options(imp)
	add_logging_options(imp, logs_dest='logs', logs_default='logs')
	parsers['db import'] = imp

	info = db_sub.add_parser('info',
		help="print the scheme information a database records",
		formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	info.add_argument('-d', '--database', metavar='DB', help="the database")
	add_logging_options(info, logs_dest='logs', logs_default='logs')
	parsers['db info'] = info

	export = db_sub.add_parser('export',
		help="export a database to text files",
		formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	export.add_argument('-d', '--database', metavar='DB', help="the database")
	export.add_argument('--format', dest='export_format', metavar='FORMAT',
		choices=database_command.supported_output,
		help="source: the files db import reads; newick: the tree only; "
			"tab is the old name of source")
	export.add_argument('-o', '--outdir', metavar='DIR',
		help="folder for the exported files (default: the current folder)")
	add_logging_options(export, logs_dest='logs', logs_default='logs')
	parsers['db export'] = export

	convert = db_sub.add_parser('convert',
		help="bring a database or canSNP file of a format before 2.3.0 to the "
			"current one",
		formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	convert.add_argument('-d', '--database', metavar='DB',
		help="the database to convert (or give --cansnps instead)")
	convert.add_argument('--cansnps', dest='annotation', metavar='FILE',
		help="the canSNP source file to convert (or give --database instead)")
	convert.add_argument('--ref-dir', dest='refdir', metavar='DIR',
		help="folder with the reference files <genome>.fna; missing files "
			"are downloaded")
	convert.add_argument('-o', '--outdir', metavar='DIR',
		help="folder for the converted file")
	add_logging_options(convert, logs_dest='logs', logs_default='logs')
	parsers['db convert'] = convert

	modify = db_sub.add_parser('modify',
		help="change the tree or the sources of an existing database",
		formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	modify.add_argument('-d', '--database', metavar='DB', help="the database")
	modify.add_argument('--mod-file', dest='mod_file', metavar='FILE',
		help="file with the tree modifications")
	modify.add_argument('--parent', metavar='NODE',
		help="node from which to update or replace")
	modify.add_argument('--replace', action='store_true',
		help="replace the node given with --parent")
	modify.add_argument('--cansnps', dest='annotation', metavar='FILE',
		help="canSNP source file with changes to the canSNP annotations")
	modify.add_argument('--references', metavar='FILE',
		help="file with reference genomes to add or update")
	add_logging_options(modify, logs_dest='logs', logs_default='logs')
	parsers['db modify'] = modify

	remove = db_sub.add_parser('remove',
		help="remove the whole branch rooted at a node",
		formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	remove.add_argument('-d', '--database', metavar='DB', help="the database")
	remove.add_argument('--parent', metavar='NODE',
		help="the branch root to remove")
	add_logging_options(remove, logs_dest='logs', logs_default='logs')
	parsers['db remove'] = remove

	'''the internal options no verb shows get their old defaults here; the options
		the verbs do show already carry their default on the option itself'''
	for name, sub in parsers.items():
		database_defaults(sub, **DB_MODES[name])
	return parsers


DB_MODES = {'db create': dict(create=True), 'db import': dict(do_import=True),
	'db info': dict(info=True), 'db export': dict(export=True),
	'db convert': dict(convert=True), 'db modify': {},
	'db remove': dict(remove=True)}


def add_scheme_subparsers(subparsers):
	scheme = subparsers.add_parser('scheme',
		help="list the schemes of CanSNPer2-data; get builds one's database")
	scheme_sub = scheme.add_subparsers(dest='scheme_action', metavar='<action>',
		required=True)

	listing = scheme_sub.add_parser('list',
		help="list the schemes and their status, and exit")
	add_scheme_source_options(listing)
	download_defaults(listing, list_schemes=True)
	parsers = {'scheme list': listing}

	get = scheme_sub.add_parser('get',
		help="fetch a scheme, build its database and download its references")
	get.add_argument('scheme', metavar='NAME', help="the scheme to build")
	get.add_argument('-d', '--database', metavar='DB',
		help="the database to build (default <scheme>.db)")
	add_scheme_source_options(get)
	get.add_argument('--overwrite', action='store_true',
		help="replace an existing database")
	get.add_argument('--no-references', dest='no_references',
		action='store_true', help="build the database without downloading "
			"its references")
	get.add_argument('--allow-unverified', dest='allow_unverified',
		action='store_true',
		help="build a verified scheme although its files cannot be checked "
			"against a manifest")
	add_logging_options(get, logs_dest='logs', logs_default='logs')
	download_defaults(get)
	parsers['scheme get'] = get
	return parsers


def add_ref_subparsers(subparsers):
	ref = subparsers.add_parser('ref',
		help="get the reference genomes a database needs from NCBI")
	ref_sub = ref.add_subparsers(dest='ref_action', metavar='<action>',
		required=True)
	get = ref_sub.add_parser('get',
		help="download the reference genomes of a database")
	get.add_argument('-d', '--database', metavar='DB',
		help="the database whose references are fetched")
	get.add_argument('-s', '--source', choices=["genbank", "refseq"],
		help="source for the downloads")
	get.add_argument('-o', '--outdir', metavar='DIR',
		help="folder for the reference genomes")
	add_logging_options(get, logs_dest='logs', logs_default='logs')
	download_defaults(get)
	return {'ref get': get}


def download_defaults(parser, **mode):
	'''every internal option of the download command with its default, plus the
		one verb of it that this subcommand runs'''
	settings = dict(database=None, scheme=None, branch=None, repo=None,
		data_dir=None, list_schemes=False, overwrite=False,
		no_references=False, allow_unverified=False, source="genbank",
		outdir="references", verbose=None, debug=None)
	settings.update(mode)
	parser.set_defaults(**settings)


def build_parsers():
	'''the top-level parser and a dict of every leaf subparser, keyed by its verb
		names, so a verb gets its own usage and error messages. The main commands
		are the top of the help and every verb names the subcommands or options it
		uses, which are shown in detail only in the help of that verb'''
	parser = argparse.ArgumentParser(
		prog='cansnper2',
		description='CanSNPer2 - one command with a verb for every job:\n'
			'type assemblies against a CanSNP database, build and maintain the\n'
			'databases, and fetch the schemes and reference genomes they need.',
		epilog='The commands CanSNPer2, CanSNPer2-database, CanSNPer2-download and\n'
			'CanSNPer2-test do the same jobs with their old arguments and keep\n'
			'working; cansnper2 runs the same code with one spelling for the shared\n'
			'options. Run cansnper2 <command> --help for the arguments of a verb;\n'
			'examples in README.md.',
		formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument('--version', action='version',
		version="CanSNPer2 - version {version}".format(version=__version__))
	subparsers = parser.add_subparsers(dest='verb', metavar='<command>',
		required=True)

	parsers = {'type': add_type_subparser(subparsers)}

	selftest_parser = subparsers.add_parser('selftest',
		help="check the installation: python imports and aligner programs "
			"in PATH")
	selftest_parser.add_argument('--debug', action='store_true',
		help="print the traceback of an unexpected error")

	parsers.update(add_db_subparsers(subparsers))
	parsers.update(add_scheme_subparsers(subparsers))
	parsers.update(add_ref_subparsers(subparsers))
	return parser, parsers


def run_command(args, parsers):
	'''call the same run function the matching old command calls'''
	if args.verb == 'type':
		typing_command.run_command(args, parsers['type'])
	elif args.verb == 'selftest':
		selftest_command.run_command()
	elif args.verb == 'db':
		database_command.run_command(args)
	elif args.verb == 'scheme':
		key = 'scheme ' + args.scheme_action
		download_command.run_command(args, parsers[key])
	elif args.verb == 'ref':
		download_command.run_command(args, parsers['ref get'])


def command():
	'''the entry point of cansnper2: one parser tree, one verb dispatch'''
	parser, parsers = build_parsers()
	args = parser.parse_args()
	run_command(args, parsers)


def main():
	'''the entry point of the command: command() with the handler of unexpected errors'''
	return run_main(command)


if __name__ == '__main__':
	main()
