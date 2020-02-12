#!/usr/bin/env python3 -c
'''
CanSNPer2 module: A toolkit for SNP-typing using NGS data.
Copyright (C) 2019 David Sundell @ FOI bioinformatics group
'''
__name__        = "CanSNPer2"
__version__     = "2.0.1"
__author__      = "David Sundell"
__credits__     = ["David Sundell"]
__license__     = "GPLv3"
__maintainer__  = "FOI bioinformatics group"
__email__       = ["bioinformatics@foi.se", "david.sundell@foi.se"]
__date__        = "2019-07-09"
__status__      = "Production"

import os
import logging
logger = logging.getLogger(__name__)

## import CanSNPer2 specific modules
from CanSNPer2.modules.ParseXMFA import ParseXMFA
from CanSNPer2.modules.NewickTree import NewickTree


## import standard python libraries for subprocess and multiprocess
from subprocess import Popen,PIPE,STDOUT
from multiprocessing import Process, Queue
from time import sleep

class Error(Exception):
	"""docstring for Error"""
	def __init__(self, value):
		self.value = value
	def __str__(self):
		return repr(self.value)

class MauveError(Error):
	"""docstring for MauveE"""
	pass





class CanSNPer2(object):
	"""docstring for CanSNPer2"""
	def __init__(self, query, mauve_path="",tmpdir=".tmp", refdir="references",database="CanSNPer.fdb", verbose=False,keep_going=False,**kwargs):
		super(CanSNPer2, self).__init__()
		self.query = query
		self.verbose = verbose
		self.refdir = refdir
		self.mauve_path = mauve_path
		self.xmfa_files = []
		self.export = kwargs["export"]
		self.snpfile = kwargs["snpfile"]
		self.database = database
		#self.initiate = kwargs["initiate"]
		#self.annotation = kwargs["annotation"]

		'''Create log and tmpdir if they do not exist'''
		self.workdir = kwargs["workdir"]
		if not os.path.exists(self.workdir):
			logger.info("Creating working directory {workdir}".format(workdir=self.workdir))
			os.makedirs(self.workdir)

		try:
			os.chdir(self.workdir)
		except FileNotFoundError:
			exit("The directory {workdir} could not be found!".format(workdir=self.workdir))

		self.tmpdir = tmpdir
		if not os.path.exists(self.tmpdir):
			logger.info("Creating tmp directory {tmpdir}".format(tmpdir=self.tmpdir))
			os.makedirs(self.tmpdir)

		self.logdir = kwargs["logdir"]
		if not os.path.exists(self.logdir):
			logger.info("Creating logs directory {logdir}".format(logdir=self.logdir))
			os.makedirs(self.logdir)

		self.outdir = kwargs["outdir"]
		if not os.path.exists(self.outdir):
			logger.info("Creating output directory {outdir}".format(outdir=self.outdir))
			os.makedirs(self.outdir)

		'''Fetch other key word arguments'''
		self.skip_mauve = kwargs["skip_mauve"]
		self.save_tree = kwargs["save_tree"]
		self.keep_temp = kwargs["keep_temp"]
		self.keep_going = keep_going

	'''CanSNPer2 get functions'''

	def get_xmfa(self):
		'''Return references to aligned files'''
		return self.xmfa_files

	def get_references(self,selected=False):
		'''references must be available in the refdir'''
		return [ref for ref in os.listdir(self.refdir) if ref.endswith(".fna")]

	def get_tempfiles(self):
		'''List all files in the tmp directory'''
		return [ref for ref in os.listdir(self.tmpdir)]

	'''ProgressiveMauve alignment functions'''

	def run_mauve(self, commands,logs):
		'''Run mauve
			ProgressiveMauve is an alignment program for small genomes
			this function will execute mauve commands as paralell subprocesses

			To run subprocesses securely this functino uses the Popen command includint standard pipes
			stderr will be passed to the stdout pipe to reduce the number of log files required

			This function currently does not support any limitation of number of processes spawned,
			modern OS will however mostly distribute subprocesses efficiently, it will also be limited by the
			number of references supplied. A warning message will be given if the number of input references
			exceeds the number of references in the database
		'''
		processes = []  #process container
		log_f = {}      #log container
		error = 0   #Variable for errors
		logger.info("Starting progressiveMauve on {n} references".format(n=len(commands)))
		for i in range(len(commands)): #Loop through commands,
			command = commands[i]
			logger.info(command)            ## In verbose mode print the actual mauve command
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

					### IF the exitcode is not 0 print a warning and ask user to read potential error messages
					if exitcode == 11:
						logger.warning("WARNING progressiveMauve finished with a exitcode: {exitcode}\nThis is typically showing up for bad genome names, script will continue...")
					elif exitcode != 0:
						if not self.keep_going:
							logger.error("Error: exitcode-{exitcode}".format(exitcode=exitcode),file=log)
						error += 1
						logger.warning("WARNING progressiveMauve finished with a non zero exitcode: {exitcode}\nThe script will terminate when all processes are finished read {log} for more info".format(log=log_f[p.stdout.fileno()],exitcode=exitcode))
					## Remove process from container
					processes.remove(p)
		if not self.keep_going:
			if error:  ## Error handling regarding mauve subprocesses, stop script if any of them fails
				logger.error("Error: {errors} progressiveMauve processes did not run correctly check log files for more information".format(errors=error))
				raise MauveError("Error: {errors} progressiveMauve processes did not run correctly check log files for more information".format(errors=error))
		else:
			return 1
		return 0

	def align(self, query, references=[]):
		'''Align sequences and run mauve as subprocess'''
		commands =[]    # store execute command
		logs = []       # store log filepath
		if len(references) == 0:    ## If specific references are not given fetch references from the reference folder
			references = self.get_references()
		for ref in references:      ## For each reference in the reference folder align to query
			ref_name = ref.rsplit(".",1)[0] ## remove file ending
			#self.query_name = os.path.basename(query).rsplit(".",1)[0] ## get name of file and remove ending
			xmfa_output = "{tmpdir}/{ref}_{target}.xmfa".format(tmpdir=self.tmpdir.rstrip("/"),ref=ref_name,target=self.query_name)
			ref_file = "{refdir}/{ref}".format(refdir=self.refdir, ref=ref)
			log_file = "{logdir}/{ref}_{target}.mauve.log".format(logdir=self.logdir,ref=ref_name,target=self.query_name)

			'''Create run command for mauve'''
			command = "{mauve_path}progressiveMauve --output {xmfa} {ref_fasta} {target_fasta}".format(
							mauve_path      = self.mauve_path,
							xmfa            = xmfa_output,
							ref_fasta       = ref_file,
							target_fasta    = query
			)
			commands.append(command)            ## Mauve command
			logs.append(log_file)               ## Store log files for each alignment
			self.xmfa_files.append(xmfa_output) ## Store the path to xmfa files as they will be used later
		if not self.skip_mauve: ### If mauve command was already run before don´t run mauve return xmfa paths
			ret = self.run_mauve(commands,logs)
			if ret != 0:
				return []
			logger.info("Alignments for {query} complete!".format(query=query))
		return self.xmfa_files

	'''Functions'''

	def create_tree(self,SNPS,name):
		'''This function uses ETE3 to color the SNP tree in the database with SNPS found in the reference database
			and outputs a pdf file
		'''
		newickTree = NewickTree(self.database,name,self.outdir)
		newickTree.draw_ete3_tree(SNPS)
		logger.info("{outdir}/{name}_tree.pdf".format(outdir =self.outdir, name=name))
		return

	def read_query_textfile_input(self,query_file):
		'''If query input is a text file, parse file'''
		with open(query_file, "r") as f:
			query = [query.strip() for query in f]
		return query

	def cleanup(self):
		'''Remove files in temporary folder'''
		logger.info("Clean up temporary files... ")
		files = self.get_tempfiles()
		for f in files:
			os.remove(os.path.join(self.tmpdir,f))
		os.rmdir(self.tmpdir)
		logger.info("Done!")
		return

	def parse_xmfa(XMFA_obj, xmfa_file, organism,results=[]):
		'''Process xmfa file using ParseXMFA object'''
		XMFA_obj.run(xmfa_file, organism)
		snps = XMFA_obj.get_snps()
		results.put(snps)
		if self.export:
			allSNP = XMFA_obj.get_snp_info()
			results.put(allSNP)
		return results

	def find_snps(self,XMFA_obj,xmfa_file,organism,results=[],export_results=[]):
		'''Align sequences to references and return SNPs'''
		XMFA_obj.run(xmfa_file, organism)
		snps = XMFA_obj.get_snps()
		results.put(snps)
		if self.export:
			allSNPs = XMFA_obj.get_snp_info()
			export_results.put(allSNPs)
		return results

	def find_snps_multiproc(self,xmfa_obj,xmfa_files,organism,export=False):
		'''function to run genomes in paralell'''
		jobs = []
		SNPS = {}
		SNP_info = []
		result_queue = Queue()
		export_queue = Queue()
		for xmfa_file in xmfa_files:
			p = Process(target=self.find_snps, args=(xmfa_obj,xmfa_file,organism ,result_queue,export_queue))
			p.start()
			jobs.append(p)
			sleep(0.05) ## A short sleep to make sure all threads do not initiate access to the database file simultanously
		for job in jobs:
			job.join()
		for j in jobs:
			## Merge SNP dictionaries
			SNPS = dict(**SNPS, **result_queue.get())
			if self.export:
				### join SNP info files
				SNP_info+= export_queue.get()
		return SNPS,SNP_info

	def run(self,database,organism):
		'''Run CanSNPer2'''
		#print("Running CanSNPer2 version-{version}".format(version=__version__))
		logger.info("Running CanSNPer2 version-{version}".format(version=__version__))

		'''Read query input file if a txt file is supplied insead of fasta files'''
		if self.query[0].endswith(".txt"):
			logger.info("Textfile input was found, parsing filepaths in {q} file".format(q=self.query[0]))
			self.query=self.read_query_textfile_input(self.query)


		'''Main function of CanSNPer2
				1. Align sequences with progressiveMauve
				2. Parse XMFA files and find SNPs
				3. Create a tree visualising SNPs found in sequence
				3b. If requested create a list with SNPs and their status
				4. Clean up tmp directory
				'''

		''' Create ParseXMFA object'''
		parse_xmfa_obj = ParseXMFA(
					database=database,
					export=self.export,
					snpfile=self.snpfile,
					verbose=self.verbose)  ## Create XMFA object and connect to database
		'''Walk through the list of queries supplied'''
		if not self.skip_mauve: print("Run {n} alignments to references using progressiveMauve".format(n=len(self.query)))
		for q in self.query:            ## For each query file_path
			qfile = q.rsplit("/")[-1]   ## Remove path from query name
			self.query_name = os.path.basename(q).rsplit(".",1)[0]  ## get name of file and remove ending

			logger.info("Running CanSNPer2 on {query}".format(query=qfile))

			if not self.skip_mauve: ### If mauve command was already run before skip step
				logger.info("Run mauve alignments")

			'''For each query fasta align to all CanSNP references the reference folder
				if skip_mauve parameter is True this the align function will only format xmfa file paths
			'''
			xmfa_files = self.align(q)

			if len(xmfa_files) == 0: ## if keep going is set and mauve exits with an error continue to next sequence
				logger.debug("Mauve exited with a non zero exit status, continue with next sample!")
				logger.warning("Mauve error skip {sample}".format(q))
				continue
			'''Parse Mauve XMFA output and find SNPs; returns SNPS (for the visual tree) and SNP_info (text file output)'''
			logger.info("Finding SNPs")
			SNPS,SNP_info = self.find_snps_multiproc(xmfa_obj=parse_xmfa_obj,xmfa_files=xmfa_files,organism=organism,export=self.export)

			'''If file export is requested print the result for each SNP location to file'''
			if self.export:
				outputfile = "{outdir}/{xmfa}_{snpfile}".format(outdir=self.outdir,xmfa=qfile.strip(".fa"),snpfile=self.snpfile)

				logger.info("Printing SNP info to {file}".format(file=outputfile))

				'''Print SNPs to tab separated file'''
				with open(outputfile,"w") as snplist_out:
					print("\t".join(["Name","Reference","Pos","Ancestral base","Derived base", "Target base"]),file=snplist_out)
					for snp in SNP_info:
						print("\t".join(snp),file=snplist_out)

			'''If save tree is requested print tree using ETE3 prints a pdf tree output'''
			if self.save_tree:
				self.create_tree(SNPS,self.query_name)

			'''Clean references to aligned xmfa files between queries if several was supplied'''
			self.xmfa_files = []

		'''Finally clean up temporary folder when all alignments and trees has been printed!'''
		if not self.keep_temp: ## if keep temp is turned on do not remove away alignments
			self.cleanup()

		logger.info("CanSNPer2 finished successfully, files can be found in {outdir}".format(outdir=self.outdir+"/"))
