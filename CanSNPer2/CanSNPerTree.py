#!/usr/bin/env python3 -c
'''
CanSNPer2: A toolkit for SNP-typing NGS data.
Copyright (C) 2019 David Sundell @ FOI bioinformatics group

VERSION 2.0.1 First release of CanSNPer2

The second release of CanSNPer (CanSNPer2) is exclusively written for python3
CanSNPer2 is simplified from CanSNPer1 stripped to only perform
required tasks. It is also written in a modular form with separate
classes for each task which allows future extentions such as a sequence read
input option planned during 2020.

Requires progressiveMauve

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
__version__ 	= "2.0.1"
__author__ 		= "David Sundell"
__credits__ 	= ["David Sundell"]
__license__ 	= "GPLv3"
__maintainer__ 	= "FOI bioinformatics group"
__email__ 		= ["bioinformatics@foi.se", "david.sundell@foi.se"]
__date__ 		= "2019-07-09"
__status__ 		= "Production"

from CanSNPer2.modules.CanSNPer2 import CanSNPer2
import logging
import logging.config
import os,sys
## Basic config file

'''CanSNPer2 settings'''
import argparse
parser = argparse.ArgumentParser(description='CanSNPer2 ')
required_arguments = parser.add_argument_group("Required arguments")
required_arguments.add_argument('query', nargs=argparse.ZERO_OR_MORE, 	metavar='query', 						help="File(s) to align (fasta)")
required_arguments.add_argument('-db',  '--database', 	metavar='', 						help='CanSNP database')
#parser.add_argument('--initiate', 			action='store_true', 				help="Use to add a CanSNP info table to FlexTaxDatabase")
#parser.add_argument('--load_snp_annotation',metavar='', default=False, 			help="On initiation supply annotation for the snps")

output_options = parser.add_argument_group("Output options")
output_options.add_argument('-o', 	 '--outdir', 	metavar='', default="results",		help="Output directory")
output_options.add_argument('--snpfile', 			metavar='', default="snpfile.txt",	help="specify name of export")
output_options.add_argument('--save_tree',			type=bool, 	default=True, 			help='Save tree as PDF using ETE3 (default True)')
output_options.add_argument('--export', 			action='store_true',				help="export a txt file with SNP results")

run_options = parser.add_argument_group("Run options")
run_options.add_argument('--refdir', 			metavar='', default="references/",	help="Specify reference directory")
run_options.add_argument('--workdir',			metavar='',	default="./",			help="Change workdir default (./)")
run_options.add_argument('--read_input', 		action='store_true', 				help="Select if input is reads not fasta")

'''Remove the two below when script is complete, possibly keep as hidden for debug'''
run_options.add_argument('--skip_mauve' ,		action='store_true', 				help="If xmfa files already exists skip step")
run_options.add_argument('--keep_going', 		action='store_true', 				help="If Error occurs, continue with the rest of samples")
run_options.add_argument('--keep_temp',			action='store_true', 				help="keep temporary files")

debugopts = parser.add_argument_group("Logging and debug options")
debugopts.add_argument('--tmpdir', 			metavar='', default="/tmp/CanSNPer2",						help="Specify reference directory")
debugopts.add_argument('--logs', 			metavar='', default="logs/", 								help="Specify log directory")
debugopts.add_argument('--verbose',			action='store_const', const=logging.INFO,					help="Verbose output")
debugopts.add_argument('--debug',			action='store_const', const=logging.DEBUG,					help="Debug output")
debugopts.add_argument('--supress',	action='store_const', const=logging.ERROR,	default=logging.WARNING,	help="Supress warnings")
parser.add_argument('--organism', 			metavar='', default="Francisella", 	help="Specify organism")  ## Will probably be depricated and decided by database selection

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

logval = args.supress
if args.debug:
	logval = args.debug
elif args.verbose:
	logval = args.verbose

from datetime import date,time
t = time()
today = date.today()
logpath = args.logs+"CanSNPer2-"+today.strftime("%b-%d-%Y")+"{}.log"
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

'''Fix temp path''' ## If tmp is in standard root folder create sub folder for CanSNPer2
if args.tmpdir.rstrip("/") == "/tmp": ## Only root tmp name was chosen
	logger.debug("Tmp directory was selected as root tmp, add subfolder for CanSNPer2 to path!")
	args.tmpdir = args.tmpdir.rstrip()+"CanSNPer2"
if not os.path.exists(args.tmpdir):
	logger.debug("Tmp directory {path} did not exist tmp, add subfolder(s) for CanSNPer2 ".format(path=args.tmpdir))
	os.makedirs(args.tmpdir,exist_ok=True)

## Add log level and log file path
logger.debug(args)

def main():
	'''Initiate CanSNPer2 object'''
	CanSNPer2_obj = CanSNPer2(args.query,
									refdir=args.refdir,
									verbose=args.verbose,
									outdir=args.outdir,
									tmpdir=args.tmpdir,
									logdir=args.logs.rstrip("/")+"/",
									skip_mauve=args.skip_mauve,
									save_tree=args.save_tree,
									keep_temp= args.keep_temp,
									workdir=args.workdir,
									export=args.export,
									snpfile=args.snpfile,
									database=args.database,
									#initiate=args.initiate,
									#annotation=args.load_snp_annotation,
					)

	'''Run CanSNPer2'''
	CanSNPer2_obj.run(database=args.database,organism=args.organism)

if oname=="__main__":
	main()
