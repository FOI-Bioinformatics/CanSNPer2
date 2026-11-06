#!/usr/bin/env python3
'''
CanSNPer2: A toolkit for SNP-typing NGS data.
Copyright (C) 2019 David Sundell @ FOI bioinformatics group

VERSION 2.1.0 release of CanSNPer2

The second release of CanSNPer (CanSNPer2) is exclusively written for python3
CanSNPer2 is simplified from CanSNPer1 stripped to only perform
required tasks. It is also written in a modular form with separate
classes for each task which allows future extentions such as a sequence read
input option planned during 2020.

Requires progressiveMauve and ete3

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
from CanSNPer2 import __version__
__author__ 		= "David Sundell"
__credits__ 	= ["David Sundell"]
__license__ 	= "GPLv3"
__maintainer__ 	= "FOI bioinformatics group"
__email__ 		= ["bioinformatics@foi.se", "david.sundell@foi.se"]
__date__ 		= "2020-04-23"
__status__ 		= "Production"

from CanSNPer2.modules.CanSNPer2 import CanSNPer2
from CanSNPer2.modules.DatabaseConnection import check_database, ConnectionError
from CanSNPer2.modules.LogSetup import setup_logging
import logging
import logging.config
import os,sys
## Basic config file

'''CanSNPer2 settings'''
import argparse


def main():
	'''Initiate CanSNPer2 object'''
	parser = argparse.ArgumentParser(description='CanSNPer2 ')
	required_arguments = parser.add_argument_group("Required arguments")
	required_arguments.add_argument('query', nargs=argparse.ZERO_OR_MORE, 	metavar='query', 		help="File(s) to align (fasta)")
	required_arguments.add_argument('-db',  '--database', 	metavar='', 							help='CanSNP database')

	output_options = parser.add_argument_group("Output options")
	output_options.add_argument('-o', '--outdir', 			metavar='DIR', default="results",		help="Output directory")
	output_options.add_argument('--save_tree',				action='store_true',					help='Save tree as PDF using ETE3 (default False)')
	output_options.add_argument('--no_snpfiles', 			action='store_false',					help="Don´t save output files.")
	output_options.add_argument('--summary',				action='store_true',					help="Output a summary file and tree with all called SNPs\nnot affected by no_snpfiles")

	run_options = parser.add_argument_group("Run options")
	run_options.add_argument('--refdir', 			metavar='', default="references/",	help="Specify reference directory")
	run_options.add_argument('--workdir',			metavar='',	default="./",			help="Change workdir default (./)")
	run_options.add_argument('--read_input', 		action='store_true', 				help="Select if input is reads not fasta")
	run_options.add_argument('--min_required_hits', type=int, default=3, 				help="Minimum sequential hits to call a SNP!")
	run_options.add_argument('--strictness', 		type=float, default=0.7,			help="Percent of snps in path reqired for calling SNP (default 0.7)")
	run_options.add_argument('--mask',				type=int, default=0,					help="Ignore SNPs closer than this many bases to the edge of an alignment (default 0, off)")
	run_options.add_argument('--keep_going', 		action='store_true', 				help="If Error occurs, continue with the rest of samples")
	run_options.add_argument('--rerun', 			action='store_true', 				help="Rerun already processed files (else skip if result file exists)")

	'''Remove the two below when script is complete, possibly keep as hidden for debug'''
	run_options.add_argument('--skip_mauve' ,		action='store_true', 				help="If xmfa files already exists skip step")
	run_options.add_argument('--keep_temp',			action='store_true', 				help="keep temporary files")


	debugopts = parser.add_argument_group("Logging and debug options")
	debugopts.add_argument('--tmpdir', 			metavar='', default="/tmp/CanSNPer2",						help="Specify reference directory")
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
		print(e)
		exit(1)

	logval = args.supress
	if args.debug:
		logval = args.debug
	elif args.verbose:
		logval = args.verbose

	if args.skip_mauve and not args.keep_temp: ## Make sure keep temp is default if --skip_mauve is used
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

	'''Fix temp path: the working folder is always a CanSNPer2 subfolder of --tmpdir.
		It used to be added only when --tmpdir was exactly /tmp, so --skip_mauve looked in
		<tmpdir> for a custom tmpdir but in <tmpdir>/CanSNPer2 for the default, and the
		cleanup at the end of a run removed everything in the folder the user had pointed
		--tmpdir at (issue #177)'''
	tmpdir = args.tmpdir.rstrip("/")
	if not tmpdir:
		parser.error("--tmpdir cannot be empty")
	if os.path.basename(tmpdir) != "CanSNPer2":
		args.tmpdir = tmpdir + "/CanSNPer2"
	else:
		args.tmpdir = tmpdir
	print("Alignments and temporary files are kept in {path}".format(path=args.tmpdir))
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
									skip_mauve=args.skip_mauve,
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
									strictness=args.strictness
	)

	'''Run CanSNPer2'''
	try:
		CanSNPer2_obj.run(database=args.database)
	except KeyboardInterrupt:
		'''Ctrl-C or SIGTERM: the alignment processes are already stopped in run_mauve, report it
			plainly instead of dumping a traceback (issue #151)'''
		logger.warning("CanSNPer2 was interrupted, the progressiveMauve processes were stopped")
		exit(130)
	'''a run that typed nothing (all NA, all skipped) is a failure and must not exit 0,
		the reason was already logged by run() (issue #191)'''
	if not CanSNPer2_obj.typed and (CanSNPer2_obj.not_typed or CanSNPer2_obj.skipped):
		exit(1)

if oname=="__main__":
	main()
