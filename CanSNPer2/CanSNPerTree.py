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
from CanSNPer2.modules.TopLevel import run_main
from CanSNPer2 import __version__, COMMANDS_HELP
__author__ 		= "David Sundell"
__credits__ 	= ["David Sundell"]
__license__ 	= "GPLv3"
__maintainer__ 	= "FOI bioinformatics group"
__email__ 		= ["bioinformatics@foi.se", "david.sundell@foi.se"]
__date__ 		= "2020-04-23"
__status__ 		= "Production"

from CanSNPer2.modules.CanSNPer2 import CanSNPer2, CanSNPer2Error, RUN_MARKER, ensure_directory
from CanSNPer2.modules.Aligners import ALIGNERS, DEFAULT_ALIGNER
from CanSNPer2.modules.ParsePileup import PLATFORMS, DEFAULT_PLATFORM
from CanSNPer2.modules.DatabaseConnection import check_database, ConnectionError
from CanSNPer2.modules.LogSetup import setup_logging
from CanSNPer2.modules.SchemeInfo import log_scheme_notice
import logging
import logging.config
import os,sys,tempfile
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


def default_tmpdir_candidates(outdir):
	'''Folders tried in turn when --tmpdir is not given: the system temporary folder (which
		follows TMPDIR), then a folder in the output directory. /tmp can be read-only, full, or hold a
		CanSNPer2 folder owned by another user, which used to stop the run.'''
	return [os.path.join(tempfile.gettempdir(), "CanSNPer2"), os.path.join(outdir, ".cansnper2_tmp")]


def fraction(value):
	'''argparse type for --strictness: a number from 0 to 1'''
	try:
		number = float(value)
	except ValueError:
		raise argparse.ArgumentTypeError("{value} is not a number".format(value=value))
	if not 0 <= number <= 1:
		raise argparse.ArgumentTypeError("{value} is not between 0 and 1".format(value=value))
	return number


def positive_int(value):
	'''argparse type for --min_required_hits: a whole number of at least 1'''
	try:
		number = int(value)
	except ValueError:
		raise argparse.ArgumentTypeError("{value} is not a whole number".format(value=value))
	if number < 1:
		raise argparse.ArgumentTypeError("{value} is below 1".format(value=value))
	return number


def non_negative_int(value):
	'''argparse type for --mask: a whole number of 0 or more'''
	try:
		number = int(value)
	except ValueError:
		raise argparse.ArgumentTypeError("{value} is not a whole number".format(value=value))
	if number < 0:
		raise argparse.ArgumentTypeError("{value} is below 0".format(value=value))
	return number


def minutes(value):
	'''argparse type for --alignment_timeout: a finite number of minutes, 0 or more. A negative
		value stopped every alignment at once, and nan meant no limit.'''
	try:
		number = float(value)
	except ValueError:
		raise argparse.ArgumentTypeError("{value} is not a number".format(value=value))
	if number != number or number in (float("inf"), float("-inf")) or number < 0:
		raise argparse.ArgumentTypeError("{value} is not a number of minutes of 0 or more".format(value=value))
	return number


def command():
	'''Initiate CanSNPer2 object'''
	parser = argparse.ArgumentParser(
		description='CanSNPer2 - type query assemblies (or reads) against a CanSNP database',
		epilog=COMMANDS_HELP,
		formatter_class=argparse.RawDescriptionHelpFormatter)
	required_arguments = parser.add_argument_group("Required arguments")
	required_arguments.add_argument('query', nargs=argparse.ZERO_OR_MORE, 	metavar='query', 		help="Assemblies (fasta) or sequencing reads (fastq, may be gzip-compressed) to type, or one list file with a file path per line (a list must be the only query argument). Paired-end read files are found by their names (_R1/_R2, _1/_2)")
	required_arguments.add_argument('-db',  '--database', 	metavar='', 							help='CanSNP database')

	output_options = parser.add_argument_group("Output options")
	output_options.add_argument('-o', '--outdir', 			metavar='DIR', default="results",		help="Output directory")
	'''the results overview (summary.html) carries the tree of every sample, so the separate
		per-sample tree pages became opt-in again: --save_tree writes them on request.
		--no_save_tree keeps working and does nothing, so command lines that give it (issue #78)
		keep working.'''
	output_options.add_argument('--save_tree',				action='store_true',					help="Also write a standalone tree page per sample as samples/<sample>_tree.html (the results overview already shows every sample's tree)")
	output_options.add_argument('--no_save_tree',			action='store_true', help=argparse.SUPPRESS)
	output_options.add_argument('--no_snpfiles', 			action='store_false',					help="Do not write the per-sample result files (samples/<sample>_snps.txt, samples/<sample>_not_called.txt); the call is still printed")
	output_options.add_argument('--summary',				action='store_true',					help="Output a summary file and tree with all called SNPs\nnot affected by no_snpfiles")
	output_options.add_argument('--annotations', metavar='FILE', default=None, help="TSV file annotating tree nodes: <snp name>\t<node|clade>\t<rank>\t<label>[\t<color>]. 'node' draws a badge under the SNP dot, 'clade' a dashed box around its whole subtree; rank is a level carried by the annotation (type, subspecies, marker...), never a tree depth")

	run_options = parser.add_argument_group("Run options")
	run_options.add_argument('--refdir', 			metavar='', default="references/",	help="Specify reference directory")
	run_options.add_argument('--allow_missing_references', action='store_true',		help="Type even when reference genomes of the database have no <genome>.fna in --refdir; one warning names them, and their SNPs are not compared (default: stop with an error)")
	run_options.add_argument('--workdir',			metavar='',	default="./",			help="Change workdir default (./)")
	'''reads are recognised by their content; --read_input is kept so that command lines that give it work'''
	run_options.add_argument('--read_input', 		action='store_true', 				help=argparse.SUPPRESS)
	run_options.add_argument('--min_required_hits', type=positive_int, default=3, 		help="Minimum number of derived SNPs on the path from the root to the final SNP for a call (at least 1, default 3)")
	run_options.add_argument('--strictness', 		type=fraction, default=0.7,			help="Fraction of the SNPs on the path to the final SNP that must be called; SNPs without alignment data count against it, an ancestral SNP stops the path (0-1, default 0.7)")
	run_options.add_argument('--mask',				type=non_negative_int, default=0,		help="Ignore SNPs closer than this many bases to the edge of an alignment (default 0, off)")
	run_options.add_argument('--keep_going', 		action='store_true', 				help="If Error occurs, continue with the rest of samples")
	run_options.add_argument('--rerun', 			action='store_true', 				help="Rerun already processed files (else skip if result file exists)")
	run_options.add_argument('--aligner',			metavar='', default=DEFAULT_ALIGNER, choices=sorted(ALIGNERS),	help="Program that aligns the query against each reference: {names} (default {default})".format(names=", ".join(sorted(ALIGNERS)),default=DEFAULT_ALIGNER))
	run_options.add_argument('--alignment_timeout',	type=minutes, default=60,	metavar='MINUTES',	help="Stop the alignments of a query that take longer than this, in minutes (0 for no limit, default 60)")
	run_options.add_argument('--aligner_path',		metavar='', default="",				help="The aligner program, or the folder that holds it (default: found in PATH)")

	read_options = parser.add_argument_group("Read typing options", "Used when the queries are sequencing reads: minimap2 maps them to each reference and samtools mpileup counts the bases at the canSNP positions. Thresholds not given take the defaults of --platform")
	read_options.add_argument('--platform',		choices=sorted(PLATFORMS), default=DEFAULT_PLATFORM,	help="Sequencing platform of the reads: illumina (short reads, paired or single), ont or hifi (default illumina)")
	read_options.add_argument('--samplesheet',	metavar='FILE', default=None,			help="Tab separated file with the header sample, reads1 and optionally reads2 and platform, one sample per row; instead of query arguments")
	read_options.add_argument('--min_depth',		type=positive_int, default=None,		help="Reads needed at a canSNP for a call (default: illumina 5, ont 5, hifi 5)")
	read_options.add_argument('--min_allele_fraction', type=fraction, default=None,		help="Fraction of the reads that must carry a base for the canSNP to be called as that base; below it the site is mixed and counts as no data (default: illumina 0.8, ont 0.7, hifi 0.8)")
	read_options.add_argument('--min_mapq',		type=non_negative_int, default=None,	help="Minimum mapping quality of a read (default 20)")
	read_options.add_argument('--min_baseq',		type=non_negative_int, default=None,	help="Minimum base quality (default: illumina 20, ont 10, hifi 20)")
	read_options.add_argument('--mixed_sites',	type=positive_int, default=3,			help="Mixed canSNP sites at which a sample is reported as possibly mixed (default 3)")
	read_options.add_argument('--threads',		type=positive_int, default=4,			help="Threads of minimap2 and samtools sort for each reference (default 4)")

	'''Remove the two below when script is complete, possibly keep as hidden for debug'''
	run_options.add_argument('--skip_alignment', '--skip_mauve', dest='skip_alignment',	action='store_true', 	help="Do not align, read the alignments already in --tmpdir (--skip_mauve is the old name)")
	run_options.add_argument('--keep_temp',			action='store_true', 				help="keep temporary files")


	debugopts = parser.add_argument_group("Logging and debug options")
	debugopts.add_argument('--tmpdir', 			metavar='', default=None,						help="Specify tmp directory (default: CanSNPer2 in the system temporary folder, see TMPDIR; a folder in --outdir is used when that cannot be written to)")
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
	if not args.database:
		print("A database is required: CanSNPer2 --database database.db query.fasta")
		exit(1)
	if not args.query and not args.samplesheet:
		'''a run without queries used to warn, exit 0 and leave an empty run directory'''
		print("No query files were given: CanSNPer2 --database database.db query.fasta", file=sys.stderr)
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

	'''--logdir, --tmpdir and --outdir are relative to --workdir. The log and temporary folders used to
		be created relative to the current folder before CanSNPer2 resolved them against the workdir,
		so the main log and an empty tmp folder ended up outside it.'''
	in_workdir = lambda path: os.path.abspath(os.path.join(args.workdir, path))
	args.logdir = in_workdir(args.logdir)
	if args.tmpdir is not None and args.tmpdir.strip():
		args.tmpdir = in_workdir(args.tmpdir)

	'''Make sure log directory path exists'''
	logdircreated = False
	if not os.path.exists(args.logdir):
		logdircreated = True
	try:
		ensure_directory(args.logdir, "--logdir")
	except CanSNPer2Error as e:
		print(e.value, file=sys.stderr)
		exit(1)
	args.logdir = args.logdir.rstrip("/")+"/" ## a logdir without a trailing slash would glue straight onto the log file name (issue #71)


	logpath = setup_logging(args.logdir, "CanSNPer2", level=logval) ## shared log naming for all commands (issue #217)
	logger = logging.getLogger(__name__)
	if logdircreated: logger.info("Creating logs directory {logdir}".format(logdir=args.logdir))
	log_scheme_notice(args.database, logger) ## which scheme is typed against, before any work starts

	'''Fix temp path, see resolve_tmpdir (issue #177)'''
	if args.tmpdir is not None and not args.tmpdir.strip():
		parser.error("--tmpdir cannot be empty")
	if args.tmpdir is None:
		'''no --tmpdir: use the first folder that can be written to'''
		problems = []
		for candidate in default_tmpdir_candidates(in_workdir(args.outdir)):
			candidate = resolve_tmpdir(candidate)
			try:
				ensure_directory(candidate, "--tmpdir")
				args.tmpdir = candidate
				break
			except CanSNPer2Error as e:
				problems.append(e.value)
				logger.warning("{0}; trying another location".format(e.value))
		else:
			logger.error("No folder for temporary files could be used ({0}); give one with --tmpdir".format("; ".join(problems)))
			exit(1)
	else:
		args.tmpdir = resolve_tmpdir(args.tmpdir)
		if args.skip_alignment and not os.path.isdir(args.tmpdir):
			'''alignments are only read with --skip_alignment: a missing folder is not created'''
			logger.error("--skip_alignment reads the alignments in --tmpdir, but {path} does not exist".format(path=args.tmpdir))
			exit(1)
		if not os.path.exists(args.tmpdir):
			logger.debug("Tmp directory {path} did not exist tmp, add subfolder(s) for CanSNPer2 ".format(path=args.tmpdir))
		try:
			ensure_directory(args.tmpdir, "--tmpdir")
		except CanSNPer2Error as e:
			logger.error(e.value)
			exit(1)

	## Add log level and log file path
	logger.debug(args)

	try:
		CanSNPer2_obj = CanSNPer2(args.query,
									refdir=args.refdir,
									allow_missing_references=args.allow_missing_references,
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
									strictness=args.strictness,
									platform=args.platform,
									samplesheet=args.samplesheet,
									min_depth=args.min_depth,
									min_allele_fraction=args.min_allele_fraction,
									min_mapq=args.min_mapq,
									min_baseq=args.min_baseq,
									mixed_sites=args.mixed_sites,
									threads=args.threads
		)
	except CanSNPer2Error as e:
		logger.error(e.value) ## an unusable output, log or tmp folder is reported before anything is aligned
		exit(1)
	if args.keep_temp:
		print("Alignments and temporary files are kept in {path}".format(path=CanSNPer2_obj.tmpdir))

	'''Run CanSNPer2'''
	try:
		CanSNPer2_obj.run(database=args.database)
	except KeyboardInterrupt:
		'''Ctrl-C or SIGTERM: the alignment processes are already stopped in run_alignments, report it
			plainly instead of dumping a traceback (issue #151)'''
		logger.warning("CanSNPer2 was interrupted, the {program} processes were stopped".format(program=args.aligner))
		if not args.keep_temp:
			CanSNPer2_obj.cleanup()
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


def main():
	'''the entry point of the command: command() with the handler of unexpected errors'''
	return run_main(command)


if oname=="__main__":
	main()
