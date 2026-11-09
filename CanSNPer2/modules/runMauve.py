import os
import shutil

'''
RunMauve aligns query genomes against every reference in a reference directory with
progressiveMauve and writes the xmfa alignments to the output directory. It is the
alignment step of CanSNPer2 on its own: no database is read and no SNPs are called
(issue #172).
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

from subprocess import Popen,STDOUT

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

		'''Create log and tmpdir if they do not exist

			the working folder is always a CanSNPer2 subfolder of --tmpdir, which is what
			CanSNPer2 --skip_mauve reads (issue #177) and it keeps cleanup from removing a
			folder that holds the user data'''
		tmpdir = tmpdir.rstrip("/")
		if os.path.basename(tmpdir) != "CanSNPer2":
			tmpdir = tmpdir + "/CanSNPer2"
		self.tmpdir = tmpdir
		print("Alignments and temporary files are kept in {tmpdir}".format(tmpdir=self.tmpdir))
		if not os.path.exists(self.tmpdir):
			if self.verbose: print("Creating tmp directory {tmpdir}".format(tmpdir=self.tmpdir))
			os.makedirs(self.tmpdir)

		self.logdir = kwargs.get("logdir","logs/") ## missing kwargs get defaults instead of a bare KeyError (issue #115)
		if not os.path.exists(self.logdir):
			if self.verbose: print("Creating logs directory {logdir}".format(logdir=self.logdir))
			os.makedirs(self.logdir)

		self.outdir = kwargs.get("outdir","results")
		if not os.path.exists(self.outdir):
			if self.verbose: print("Creating output directory {outdir}".format(outdir=self.outdir))
			os.makedirs(self.outdir)

		'''Fetch other key word arguments (save_tree is accepted for compatibility with the
			CanSNPer2 call pattern but this module only aligns, it never draws trees)'''
		self.skip_mauve = kwargs.get("skip_mauve",False)
		self.keep_temp = kwargs.get("keep_temp",False)

	def get_xmfa(self):
		'''Return references to aligned files'''
		return self.xmfa_files

	def get_references(self,selected=False):
		'''references must be available in the refdir, both .fna (what CanSNPer2-download writes)
			and the older .fa endings are aligned, an empty folder used to yield an empty command
			list so the run finished silently without producing a single alignment (issue #172)'''
		try:
			refs = sorted(ref for ref in os.listdir(self.refdir) if ref.endswith((".fna",".fa",".fasta")))
		except OSError as e:
			raise MauveError("The reference directory {refdir} could not be read: {err}".format(refdir=self.refdir,err=e))
		if not refs:
			raise MauveError("No reference genomes (.fna/.fa) were found in {refdir}, run CanSNPer2-download --outdir {refdir} first".format(refdir=self.refdir))
		return refs

	def get_tempfiles(self):
		'''List all files in the tmp directory'''
		return [ref for ref in os.listdir(self.tmpdir)]

	def run_mauve(self, commands,logs):
		'''Execute mauve commands'''
		processes = [] 	#process container
		logfiles = []	#open log file container, parallel to processes
		error = False	#Variable for errors
		for i in range(len(commands)): #Loop through commands,
			command = commands[i]
			if self.verbose: print(command)			## In verbose mode print the actual mauve command
			## stdout/stderr go straight to the log file: a PIPE fills up and blocks the
			## process while this loop only polls exitcodes, which hangs the run (issue #85)
			log = open(logs[i], "w")
			p = Popen(command,  stdout=log, stderr=STDOUT)  ##command is an argument list, no shell and no fragile splitting (issue #85)
			logfiles.append(log)
			processes.append(p)
		while processes:
			for idx in range(len(processes)-1, -1, -1): ## Loop through processes, walking backwards so removals are safe
				p = processes[idx]
				exitcode = p.poll()  #get exitcode for each process
				if exitcode is not None: ## if there is no exitcode the program is still running

					## When a process is finished close its log file
					logfiles[idx].close()

					### IF the exitcode is not 0 print a warning and ask user to read
					if exitcode != 0:
						error = True
						print("WARNING progressiveMauve finished with a non zero exitcode script will terminate when all processes are finished read {log} for more info".format(log=logs[idx]))
					## Remove process from container
					processes.pop(idx)
					logfiles.pop(idx)
		if error:  ## Error handling regarding mauve subprocesses, stop script if any of them fails
			raise MauveError("At least of the progressiveMauve processes did not run correclty check the log files for more information")
		return

	def stage_for_mauve(self,path,prefix):
		'''progressiveMauve writes a multi megabyte "<input>.sslist" progress file next to every input
			genome, which pollutes the query and reference folders and makes them required to be
			writable; same fix as CanSNPer2: symlink into the tmpdir and let cleanup remove it (issue #148)'''
		source = os.path.abspath(path)
		staged = os.path.join(self.tmpdir, "{prefix}_{name}".format(prefix=prefix,name=os.path.basename(source)))
		if not os.path.lexists(staged):
			try:
				os.symlink(source,staged)
			except (OSError,NotImplementedError,AttributeError):
				'''a tmpdir on a filesystem without symlinks, fall back to a copy'''
				if self.verbose: print("Could not symlink {source} into {tmpdir}, copying the genome instead".format(source=source,tmpdir=self.tmpdir))
				shutil.copyfile(source,staged)
		return staged

	def align(self, query, references=[]):
		'''Run mauve as a subprocess'''
		commands =[]  	# store execute command
		logs = []		# store log filepath
		if len(references) == 0: 	## If specific references are not given fetch references from the reference folder
			references = self.get_references()
		query_file = self.stage_for_mauve(query,"query")
		for ref in references:  	## For each reference in the reference folder align to query
			ref_name = ref.rsplit(".",1)[0] ## remove file ending
			self.query_name = os.path.basename(query).rsplit(".",1)[0]	## get name of file and remove ending
			xmfa_output = "{tmpdir}/{ref}_{target}.xmfa".format(tmpdir=self.tmpdir.rstrip("/"),ref=ref_name,target=self.query_name)
			ref_file = self.stage_for_mauve("{refdir}/{ref}".format(refdir=self.refdir, ref=ref),"ref")
			log_file = "{logdir}/{ref}_{target}.mauve.log".format(logdir=self.logdir,ref=ref_name,target=self.query_name)

			'''Create run command for mauve as an argument list, so paths with spaces stay single arguments (issue #85)'''
			command = ["{mauve_path}progressiveMauve".format(mauve_path=self.mauve_path),
							"--output", xmfa_output,
							ref_file,
							query_file]
			commands.append(command)			## Mauve command
			logs.append(log_file)				## Store log files for each alignment
			self.xmfa_files.append(xmfa_output)	## Store the path to xmfa files as they will be used later
		if not self.skip_mauve: ### If mauve command was already run before skip step
			self.run_mauve(commands,logs)
		if self.verbose: print("Alignments for {query} complete!".format(query=query))
		return

	def read_query_input(self,query_file):
		'''If query input is a text file, parse file

			query_file is the query list from the command line, so the text file is its
			first entry; passing the list itself to open() raised TypeError and the
			documented "one fasta per line" input never worked (issue #181). Blank lines
			are skipped, they would otherwise fail the file existence check.'''
		if isinstance(query_file,list):
			query_file = query_file[0]
		with open(query_file, "r") as f:
			query = [query.strip() for query in f if query.strip()]
		if not query:
			raise MauveError("The query input file {query_file} contained no filepaths".format(query_file=query_file))
		return query

	def cleanup(self):
		'''Remove files in temporary folder'''
		if self.verbose: print("Clean up temporary files... ")
		files = self.get_tempfiles()
		for f in files:
			path = os.path.join(self.tmpdir,f)
			if os.path.isdir(path): ## sub directories would break os.remove and leave tmpdir non-empty for os.rmdir (issue #86/#87)
				shutil.rmtree(path)
			else:
				os.remove(path)
		shutil.rmtree(self.tmpdir, ignore_errors=True)
		if self.verbose: print("Done!")
		return

	def run(self,database=False):
		'''Run RunMauve

			Alignment only: every query is aligned against every reference in the refdir and the
			xmfa files are written to outdir. SNP calling is not done here (CanSNPer2 does that),
			so the database and the organism flag are kept for compatibility with the CanSNPer2
			call pattern and are not used, and the alignments are copied out of the tmpdir before
			cleanup, which used to remove them and leave the run with nothing at all (issue #172)'''
		print("Run alignments to references using progressiveMauve")

		'''Read query input file if a txt file is supplied insead of fasta files'''
		if not self.query:
			raise MauveError("No query input was given, nothing to align")
		if self.query[0].endswith(".txt"):
			self.query=self.read_query_input(self.query)

		for q in self.query:  ## For each target file
			if not os.path.exists(q):
				'''mauve would only fail in a subprocess and the log file is easy to miss (issue #172)'''
				raise MauveError("Input file: {q} was not found!".format(q=q))
			if self.verbose: print("Running alignments on query: {query}".format(query=q))
			'''For each query fasta align to all CanSNPer reference for selected organism'''
			self.align(q)

		aligned = []
		for xmfa in self.xmfa_files:
			if not os.path.exists(xmfa):
				print("WARNING no alignment written for {xmfa}, check the mauve logs in {logdir}".format(xmfa=xmfa,logdir=self.logdir))
				continue
			target = os.path.join(self.outdir,os.path.basename(xmfa))
			if os.path.abspath(target) != os.path.abspath(xmfa):
				shutil.copyfile(xmfa,target)
			aligned.append(target)
		if not aligned:
			raise MauveError("No alignments were produced, check the mauve logs in {logdir}".format(logdir=self.logdir))
		print("{n} alignments written to {outdir}".format(n=len(aligned),outdir=self.outdir))
		'''Clean up temporary folder when all alignments has been printed!'''
		if not self.keep_temp:
			self.cleanup()
		return aligned

if __name__=="__main__":  #main should be moved to the file __main__ for setup.py
	import argparse
	parser = argparse.ArgumentParser(description='Align query genomes against all CanSNPer2 references with progressiveMauve and write the xmfa files to outdir (alignment only, no SNP calling)')
	parser.add_argument('query', nargs='+', metavar='', help="Takes multiple fasta files or a txt file with one fasta file per line")
	parser.add_argument('-db','--database', 	metavar='', help='CanSNP database (kept for compatibility, alignment only so not used)')
	parser.add_argument('--save_tree',	type=bool, default=True , help='Save tree as PDF using ETE3 (not used, CanSNPer2 creates the trees)')
	parser.add_argument('--organism', 	metavar='', default="Francisella", 	help="Specify organism (not used, the refdir selects the references)")
	parser.add_argument('--refdir', 	metavar='', default="references",		help="Specify reference directory")
	parser.add_argument('--outdir', 	metavar='', default="results",	help="Directory the xmfa alignment files are written to")
	parser.add_argument('--tmpdir', 	metavar='', default=".tmp",		help="Specify temporary directory")
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
	try:
		CanSNPer.run(args.database)
	except MauveError as e:
		'''a missing reference folder or a bad query path is a usage error, not a traceback (issue #172)'''
		print("ERROR: {err}".format(err=e))
		exit(1)
