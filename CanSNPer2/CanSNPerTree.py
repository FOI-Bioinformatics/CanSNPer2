#!/usr/bin/env python3
'''
CanSNPer2: A toolkit for SNP-typing NGS data.
Copyright (C) 2019 David Sundell @ FOI bioinformatics group

VERSION 2.2.2 release of CanSNPer2

The second release of CanSNPer (CanSNPer2) is exclusively written for python3
CanSNPer2 is simplified from CanSNPer1 stripped to only perform
required tasks. It is also written in a modular form with separate
classes for each task which allows future extentions such as a sequence read
input option planned during 2020.

Requires an aligner, by default progressiveMauve

This program is free software: you can redistribute and/or modify
the code under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''

oname = __name__
__name__ 		= "CanSNPer2"
from CanSNPer2 import __version__, COMMANDS_HELP
__author__ 		= "David Sundell"
__credits__ 	= ["David Sundell"]
__license__ 	= "GPLv3"
__maintainer__ 	= "FOI bioinformatics group"
__email__ 		= ["bioinformatics@foi.se", "david.sundell@foi.se"]
__date__ 		= "2020-04-23"
__status__ 		= "Production"

from CanSNPer2.modules.CanSNPer2 import CanSNPer2, CanSNPer2Error, RUN_MARKER
from CanSNPer2.modules.Aligners import ALIGNERS, DEFAULT_ALIGNER
from CanSNPer2.modules.DatabaseConnection import check_database, ConnectionError
from CanSNPer2.modules.LogSetup import setup_logging
from CanSNPer2.modules.SchemeInfo import log_scheme_notice
import logging
import logging.config
import os,sys
## Basic config file

'''CanSNPer2 settings'''
import argparse


def resolve_tmpdir(path):
	'''The working folder is a CanSNPer2 subfolder of --tmpdir. It used to be added only
		when --tmpdir was exactly /tmp, so --skip_alignment looked in <tmpdir> for a custom tmpdir
		but in <tmpdir>/CanSNPer2 for the default (issue #177). A run directory that CanSNPer2
		created and kept (--keep_temp) carries a marker file and is used as it is, so its
		alignments can be read again with --skip_alignment.'''
	tmpdir = path.rstrip("/")
	if os.path.basename(tmpdir) == "CanSNPer2" or os.path.isfile(os.path.join(tmpdir, RUN_MARKER)):
		return tmpdir
	return tmpdir + "/CanSNPer2"


def fraction(value):
	'''argparse type for --strictness: a number from 0 to 1'''
	try:
		number = float(value)
	except ValueError:
		raise argparse.ArgumentTypeError("{value} is not a number".format(value=value))
	if not 0 <= number <= 1:
		raise argparse.ArgumentTypeError("{value} is not between 0 and 1".format(value=value))
	return number


def main():
	'''Initiate CanSNPer2 object'''
	parser = argparse.ArgumentParser(
		description='CanSNPer2 - type query assemblies (or reads) against a CanSNP database',
		epilog=COMMANDS_HELP,
		formatter_class=argparse.RawDescriptionHelpFormatter)
	required_arguments = parser.add_argument_group("Required arguments")
	required_arguments.add_argument('query', nargs=argparse.ZERO_OR_MORE, 	metavar='query', 		help="File(s) to align (fasta)")
	required_arguments.add_argument('-db',  '--database', 	metavar='', 							help='CanSNP database')

	output_options = parser.add_argument_group("Output options")
	output_options.add_argument('-o', '--outdir', 			metavar='DIR', default="results",		help="Output directory")
	output_options.add_argument('--save_tree',				action='store_true',					help='Save tree as PDF using ETE3 (default False)')
	output_options.add_argument('--no_snpfiles', 			action='store_false',					help="Don´t save output files.")
	output_options.add_argument('--summary',				action='store_true',					help="Output a summary file and tree with all called SNPs\nnot affected by no_snpfiles")
	output_options.add_argument('--annotations', metavar='FILE', default=None, help="TSV file annotating tree nodes: <snp name>\t<node|clade>\t<rank>\t<label>[\t<color>]. 'node' draws a badge under the SNP dot, 'clade' a dashed box around its whole subtree; rank is a level carried by the annotation (type, subspecies, marker...), never a tree depth")

	run_options = parser.add_argument_group("Run options")
	run_options.add_argument('--refdir', 			metavar='', default="references/",	help="Specify reference directory")
	run_options.add_argument('--workdir',			metavar='',	default="./",			help="Change workdir default (./)")
	run_options.add_argument('--read_input', 		action='store_true', 				help="Select if input is reads not fasta")
	run_options.add_argument('--min_required_hits', type=int, default=3, 				help="Minimum sequential hits to call a SNP!")
	run_options.add_argument('--strictness', 		type=fraction, default=0.7,			help="Fraction of the SNPs on the path to the final SNP that must be called; SNPs without alignment data count against it, an ancestral SNP stops the path (0-1, default 0.7)")
	run_options.add_argument('--mask',				type=int, default=0,					help="Ignore SNPs closer than this many bases to the edge of an alignment (default 0, off)")
	run_options.add_argument('--keep_going', 		action='store_true', 				help="If Error occurs, continue with the rest of samples")
	run_options.add_argument('--rerun', 			action='store_true', 				help="Rerun already processed files (else skip if result file exists)")
	run_options.add_argument('--aligner',			metavar='', default=DEFAULT_ALIGNER, choices=sorted(ALIGNERS),	help="Program that aligns the query against each reference: {names} (default {default})".format(names=", ".join(sorted(ALIGNERS)),default=DEFAULT_ALIGNER))
	run_options.add_argument('--alignment_timeout',	type=float, default=60,	metavar='MINUTES',	help="Stop the alignments of a query that take longer than this, in minutes (0 for no limit, default 60)")
	run_options.add_argument('--aligner_path',		metavar='', default="",				help="The aligner program, or the folder that holds it (default: found in PATH)")

	'''Remove the two below when script is complete, possibly keep as hidden for debug'''
	run_options.add_argument('--skip_alignment', '--skip_mauve', dest='skip_alignment',	action='store_true', 	help="Do not align, read the alignments already in --tmpdir (--skip_mauve is the old name)")
	run_options.add_argument('--keep_temp',			action='store_true', 				help="keep temporary files")


	debugopts = parser.add_argument_group("Logging and debug options")
	debugopts.add_argument('--tmpdir', 			metavar='', default="/tmp/CanSNPer2",						help="Specify tmp directory default (/tmp/CanSNPer2)")
	debugopts.add_argument('--logdir', 			metavar='', default="logs/", 								help="Specify log directory")
	debugopts.add_argument('--verbose',			action='store_const', const=logging.INFO,					help="Verbose output")
	debugopts.add_argument('--debug',			action='store_const', const=logging.DEBUG,					help="Debug output")
	debugopts.add_argument('--supress',	action='store_const', const=logging.ERROR,	default=logging.WARNING,	help="Supress warnings")
	parser.add_argument("--version", action='store_true', help=argparse.SUPPRESS)

	args = parser.parse_args()
	if len(sys.argv)==1:
		parser.print_help()
		parser.exit()
	'''
		Setup logging and debug options
	'''
	if args.version:
		print("CanSNPer2 - version {version}".format(version=__version__))
		exit()
	if args.read_input:
		'''a flag that is not implemented must not leave with status 0: a script that passed
		--read_input looked like it had typed something (issue #227)'''
		print("Read input is not implemented in this version, fasta input is required.")
		exit(1)
	if not args.database:
		print("A database is required: CanSNPer2 --database database.db query.fasta")
		exit(1)
	'''fail with a clear message on a bad database path instead of creating an empty database'''
	try:
		check_database(args.database)
	except ConnectionError as e:
		print(e, file=sys.stderr) ## an error message, so a script that reads stdout does not take it as output
		exit(1)

	logval = args.supress
	if args.debug:
		logval = args.debug
	elif args.verbose:
		logval = args.verbose

	if args.skip_alignment and not args.keep_temp: ## Make sure keep temp is default if --skip_alignment is used
		args.keep_temp = True

	'''Make sure log directory path exists'''
	logdircreated = False
	if not os.path.exists(args.logdir):
		logdircreated = True
		os.makedirs(args.logdir)
	args.logdir = args.logdir.rstrip("/")+"/" ## a logdir without a trailing slash would glue straight onto the log file name (issue #71)


	logpath = setup_logging(args.logdir, "CanSNPer2", level=logval) ## shared log naming for all commands (issue #217)
	logger = logging.getLogger(__name__)
	if logdircreated: logger.info("Creating logs directory {logdir}".format(logdir=args.logdir))
	log_scheme_notice(args.database, logger) ## which scheme is typed against, before any work starts

	'''Fix temp path, see resolve_tmpdir (issue #177)'''
	if not args.tmpdir.rstrip("/"):
		parser.error("--tmpdir cannot be empty")
	args.tmpdir = resolve_tmpdir(args.tmpdir)
	if not os.path.exists(args.tmpdir):
		logger.debug("Tmp directory {path} did not exist tmp, add subfolder(s) for CanSNPer2 ".format(path=args.tmpdir))
		os.makedirs(args.tmpdir,exist_ok=True)

	## Add log level and log file path
	logger.debug(args)

	CanSNPer2_obj = CanSNPer2(args.query,
									refdir=args.refdir,
									verbose=args.verbose,
									outdir=args.outdir,
									tmpdir=args.tmpdir,
									logdir=args.logdir.rstrip("/")+"/",
									skip_alignment=args.skip_alignment,
									aligner=args.aligner,
									aligner_path=args.aligner_path,
									alignment_timeout=args.alignment_timeout,
									save_tree=args.save_tree,
									keep_temp= args.keep_temp,
									workdir=args.workdir,
									export=args.no_snpfiles,
									database=args.database,
									min_required_hits=args.min_required_hits,
									keep_going=args.keep_going,
									mask=args.mask,
									rerun=args.rerun,
									summary=args.summary,
									annotations=args.annotations,
									strictness=args.strictness
	)
	print("Alignments and temporary files are kept in {path}".format(path=CanSNPer2_obj.tmpdir))

	'''Run CanSNPer2'''
	try:
		CanSNPer2_obj.run(database=args.database)
	except KeyboardInterrupt:
		'''Ctrl-C or SIGTERM: the alignment processes are already stopped in run_alignments, report it
			plainly instead of dumping a traceback (issue #151)'''
		logger.warning("CanSNPer2 was interrupted, the {program} processes were stopped".format(program=args.aligner))
		exit(130)
	except CanSNPer2Error as e:
		'''a refused run (for example two samples with the same name) is reported as a message,
			not as a traceback, and leaves no run directory behind'''
		logger.error(e.value)
		if not args.keep_temp:
			CanSNPer2_obj.cleanup()
		exit(1)
	'''a run that typed nothing (all NA, all skipped) is a failure and must not exit 0,
		the reason was already logged by run() (issue #191); a sample whose earlier result
		was used again, made from the same input, counts as typed'''
	if not CanSNPer2_obj.typed and not CanSNPer2_obj.reused and (CanSNPer2_obj.not_typed or CanSNPer2_obj.skipped):
		exit(1)

if oname=="__main__":
	main()
