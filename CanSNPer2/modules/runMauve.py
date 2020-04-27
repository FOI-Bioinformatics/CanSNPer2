import sys, os
import sqlite3

'''
DatabaseConnections and DatabaseFunctions are classes written to simplify database work using sqlite3
This particular script is specially adapted to suit the use for CanSNPer web database!
'''

__version__ = "0.0.1"
__author__ = "David Sundell"
__credits__ = ["David Sundell"]
__license__ = "GPLv3"
__maintainer__ = "FOI bioinformatics group"
__email__ = ["bioinformatics@foi.se", "david.sundell@foi.se"]
__date__ = "2019-04-19"
__status__ = "Production"
__partof__ = "CanSNPer2"

from subprocess import Popen,PIPE,STDOUT
from ParseXMFA import ParseXMFA

class MauveError(Exception):
	"""docstring for MauveE"""
	def __init__(self, value):
		self.value = value
	def __str__(self):
		return repr(self.value)


class RunMauve(object):
	"""docstring for RunMauve"""
	def __init__(self, query, mauve_path="",tmpdir=".tmp", refdir="references",verbose=False,**kwargs):
		super(RunMauve, self).__init__()
		self.query = query
		self.verbose = verbose
		self.refdir = refdir
		self.mauve_path = mauve_path
		self.xmfa_files = []

		'''Create log and tmpdir if they do not exist'''
		self.tmpdir = tmpdir
		if not os.path.exists(self.tmpdir):
			if self.verbose: print("Creating tmp directory {tmpdir}".format(tmpdir=self.tmpdir))
			os.makedirs(self.tmpdir)

		self.logdir = kwargs["logdir"]
		if not os.path.exists(self.tmpdir):
			if self.verbose: print("Creating logs directory {tmpdir}".format(tmpdir=self.tmpdir))
			os.makedirs(self.tmpdir)

		self.outdir = kwargs["outdir"]
		if not os.path.exists(self.outdir):
			if self.verbose: print("Creating output directory {outdir}".format(outdir=self.outdir))
			os.makedirs(self.outdir)

		'''Fetch other key word arguments'''
		self.skip_mauve = kwargs["skip_mauve"]
		self.save_tree = kwargs["save_tree"]
		self.keep_temp = kwargs["keep_temp"]

	def get_xmfa(self):
		'''Return references to aligned files'''
		return self.xmfa_files

	def get_references(self,selected=False):
		'''references must be available in the refdir'''
		return [ref for ref in os.listdir(self.refdir) if ref.endswith(".fa")]

	def get_tempfiles(self):
		'''List all files in the tmp directory'''
		return [ref for ref in os.listdir(self.tmpdir)]

	def run_mauve(self, commands,logs):
		'''Execute mauve commands'''
		processes = [] 	#process container
		log_f = {}		#log container
		error = False	#Variable for errors
		for i in range(len(commands)): #Loop through commands,
			command = commands[i]
			if self.verbose: print(command)			## In verbose mode print the actual mauve command
			p = Popen(command.split(" "),  stdout=PIPE, stderr=STDOUT)  ##Split command to avoid shell=True pipe stdout and stderr to stdout
			log_f[p.stdout.fileno()] = logs[i]      ## Store the reference to the correct log file
			processes.append(p)
		while processes:
			for p in processes: ## Loop through processes
				exitcode = p.poll()  #get exitcode for each process
				if exitcode is not None: ## if there is no exitcode the program is still running

					## When a process is finished open the log file and write stdout/stderr to log file
					with open(log_f[p.stdout.fileno()], "w") as log:
						print(p.stdout.read().decode('utf-8'),file=log)

					### IF the exitcode is not 0 print a warning and ask user to read
					if exitcode != 0:
						error = True
						print("WARNING progressiveMauve finished with a non zero exitcode script will terminate when all processes are finished read {log} for more info".format(log=log_f[p.stdout.fileno()]))
					## Remove process from container
					processes.remove(p)
		if error:  ## Error handling regarding mauve subprocesses, stop script if any of them fails
			raise MauveError("At least of the progressiveMauve processes did not run correclty check the log files for more information")
		return

	def align(self, query, references=[]):
		'''Run mauve as a subprocess'''
		commands =[]  	# store execute command
		logs = []		# store log filepath
		if len(references) == 0: 	## If specific references are not given fetch references from the reference folder
			references = self.get_references()
		for ref in references:  	## For each reference in the reference folder align to query
			ref_name = ref.rsplit(".",1)[0] ## remove file ending
			self.query_name = os.path.basename(query).rsplit(".",1)[0]	## get name of file and remove ending
			xmfa_output = "{tmpdir}/{ref}_{target}.xmfa".format(tmpdir=self.tmpdir.rstrip("/"),ref=ref_name,target=self.query_name)
			ref_file = "{refdir}/{ref}".format(refdir=self.refdir, ref=ref)
			log_file = "{logdir}/{ref}_{target}.mauve.log".format(logdir=self.logdir,ref=ref_name,target=self.query_name)

			'''Create run command for mauve'''
			command = "{mauve_path}progressiveMauve --output {xmfa} {ref_fasta} {target_fasta}".format(
							mauve_path 		= self.mauve_path,
							xmfa 			= xmfa_output,
							ref_fasta 		= ref_file,
							target_fasta 	= query
			)
			commands.append(command)			## Mauve command
			logs.append(log_file)				## Store log files for each alignment
			self.xmfa_files.append(xmfa_output)	## Store the path to xmfa files as they will be used later
		if not self.skip_mauve: ### If mauve command was already run before skip step
			self.run_mauve(commands,logs)
		if self.verbose: print("Alignments for {query} complete!".format(query=query))
		return

	def read_query_input(self,query_file):
		'''If query input is a text file, parse file'''
		with open(query_file, "r") as f:
			query = [query.strip() for query in f]
		return query

	def run(self,database,organism):
		'''Run RunMauve'''
		print("Run alignments to references using progressiveMauve")

		'''Read query input file if a txt file is supplied insead of fasta files'''
		if self.query[0].endswith(".txt"):
			self.query=self.read_query_input(self.query)

		'''Main function of RunMauve
				1. Align sequences with progressiveMauve
				2. Parse XMFA files for SNPs
				3. Create a tree visualising which SNPs were found
				3b. If requested create a list with SNPs
				4. Clean up tmp directory
				'''
		xmfa = ParseXMFA(database=args.database,verbose=self.verbose)  ## Create XMFA object and connect to database
		for q in self.query:  ## For each target file
			if self.verbose: print("Running alignments on query: {query}".format(query=q))
			'''For each query fasta align to all CanSNPer reference for selected organism'''
			self.align(q)

			# '''Find all SNPs found in each reference'''
			# self.find_snps(xmfa=xmfa,organism=args.organism,save_tree=self.save_tree)

			# '''Clean references to aligned xmfa files between queries'''
			# self.xmfa_files = []
		'''Clean up temporary folder when all alignments and trees has been printed!'''
		if not self.keep_temp:
			self.cleanup()
		return self.xmfa_files

if __name__=="__main__":  #main should be moved to the file __main__ for setup.py
	import argparse
	parser = argparse.ArgumentParser(description='Parse xmfa files and extract non overlapping sequences')
	parser.add_argument('query', nargs='+', metavar='', help="Takes multiple fasta files or a txt file with one fasta file per line")
	parser.add_argument('-db','--database', 	metavar='', help='CanSNP database')
	parser.add_argument('--save_tree',	type=bool, default=True , help='Save tree as PDF using ETE3 (default True)')
	parser.add_argument('--organism', 	metavar='', default="Francisella", 	help="Specify organism")
	parser.add_argument('--refdir', 	metavar='', default="references",		help="Specify reference directory")
	parser.add_argument('--outdir', 	metavar='', default="results",	help="Specify output directory")
	parser.add_argument('--tmpdir', 	metavar='', default=".tmp",		help="Specify reference directory")
	parser.add_argument('--logs', 		metavar='', default="logs", 	help="Specify log directory")
	parser.add_argument('--verbose',	action='store_true',help="print process info, default no output")
	parser.add_argument('--skip_mauve' ,action='store_true', help="If xmfa files already exists skip step")
	parser.add_argument('--keep_temp' ,action='store_true', help="keep temporary files")
	args = parser.parse_args()
	if args.verbose: print(args)

	CanSNPer = RunMauve(args.query,
								refdir=args.refdir,
								verbose=args.verbose,
								outdir=args.outdir,
								tmpdir=args.tmpdir,
								logdir=args.logs,
								skip_mauve=args.skip_mauve,
								save_tree=args.save_tree,
								keep_temp= args.keep_temp,
								)
	CanSNPer.run(args.database,args.organism)
