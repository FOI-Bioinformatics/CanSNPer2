'''
CanSNPer2 module: A toolkit for SNP-typing using NGS data.
Copyright (C) 2019 David Sundell @ FOI bioinformatics group
'''

import os
import logging
logger = logging.getLogger(__name__)

## import CanSNPer2 specific modules
from CanSNPer2.modules.ParseXMFA import ParseXMFA
from CanSNPer2.modules.NewickTree import NewickTree
from CanSNPer2.CanSNPerTree import __version__


## import standard python libraries for subprocess and multiprocess
from subprocess import Popen,PIPE,STDOUT
from multiprocessing import Process, Queue
from time import sleep,time

class Error(Exception):
	"""docstring for Error"""
	def __init__(self, value):
		self.value = value
	def __str__(self):
		return repr(self.value)

class MauveError(Error):
	"""docstring for MauveE"""
	pass

class CanSNPer2Error(Error):
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
		self.database = database
		self.min_required_hits = kwargs["min_required_hits"]
		self.strictness = kwargs["strictness"]
		self.rerun = kwargs["rerun"]

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

		if kwargs["summary"]:
			self.summary_set = set()
			self.called_genome = {}
			self.summary = True
		else:
			self.summary = False

		self.no_export = False

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
		retvalue=0
		processes = []  #process container
		log_f = {}	  #log container
		error = 0   #Variable for errors
		logger.info("Starting progressiveMauve on {n} references".format(n=len(commands)))
		for i in range(len(commands)): #Loop through commands,
			command = commands[i]
			logger.debug(command)			## In verbose mode print the actual mauve command
			p = Popen(command.split(" "),  stdout=PIPE, stderr=STDOUT)  ##Split command to avoid shell=True pipe stdout and stderr to stdout
			log_f[p.stdout.fileno()] = logs[i]	  ## Store the reference to the correct log file
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
						logger.warning("WARNING progressiveMauve finished with a exitcode: {exitcode}".format(exitcode=exitcode))
						logger.debug("This progressiveMauve error is showing up for bad genomes containing short repetative contigs or sequence contains dashes".format(exitcode=exitcode))
						retvalue=11
					elif exitcode == -6:
						logger.warning("Input sequence is not free of gaps, replace gaps with N and retry!!")
						retvalue=6
					elif exitcode != 0:
						if not self.keep_going:
							logger.error("Error: exitcode-{exitcode}".format(exitcode=exitcode))
						error += 1
						logger.warning("WARNING progressiveMauve finished with a non zero exitcode: {exitcode}\nThe script will terminate when all processes are finished read {log} for more info".format(log=log_f[p.stdout.fileno()],exitcode=exitcode))
					## Remove process from container
					processes.remove(p)
		if retvalue == 6 or retvalue == 11: ## This is done if sequence is not free of gaps
			return retvalue
		elif not self.keep_going:
			if error:  ## Error handling regarding mauve subprocesses, stop script if any of them fails
				logger.error("Error: {errors} progressiveMauve processes did not run correctly check log files for more information".format(errors=error))
				raise MauveError("Error: {errors} progressiveMauve processes did not run correctly check log files for more information".format(errors=error))
		else:
			retvalue=1
		return retvalue

	def create_mauve_command(self,query,references=[]):
		'''Mauve commands'''
		commands =[]	# store execute command
		logs = []	   # store log filepath
		if len(references) == 0:	## If specific references are not given fetch references from the reference folder
			references = self.get_references()
		for ref in references:	  ## For each reference in the reference folder align to query
			ref_name = ref.rsplit(".",1)[0] ## remove file ending
			#self.query_name = os.path.basename(query).rsplit(".",1)[0] ## get name of file and remove ending
			xmfa_output = "{tmpdir}/{ref}_{target}.xmfa".format(tmpdir=self.tmpdir.rstrip("/"),ref=ref_name,target=self.query_name)
			ref_file = "{refdir}/{ref}".format(refdir=self.refdir, ref=ref)
			log_file = "{logdir}/{ref}_{target}.mauve.log".format(logdir=self.logdir,ref=ref_name,target=self.query_name)

			'''Create run command for mauve'''
			command = "{mauve_path}progressiveMauve --output {xmfa} {ref_fasta} {target_fasta}".format(
							mauve_path	  = self.mauve_path,
							xmfa			= xmfa_output,
							ref_fasta	   = ref_file,
							target_fasta	= query
			)
			commands.append(command)			## Mauve command
			logs.append(log_file)			   ## Store log files for each alignment
			self.xmfa_files.append(xmfa_output) ## Store the path to xmfa files as they will be used later
		return commands,logs

	def align(self, query, references=[]):
		'''Align sequences and run mauve as subprocess'''
		commands,logs = self.create_mauve_command(query,references)
		if not self.skip_mauve: ### If mauve command was already run before don´t run mauve return xmfa paths
			ret = self.run_mauve(commands,logs)
			if int(ret) == 11:
				'''This error might be caused by mauve not handling dash(-) in the sequence, change the incoming sequence and replace - with N and retry once'''
				logger.warning("Mauve ran into an error with sequence {seq} may contain a dash, replace with N characters and retry mauve".format(seq=query))
				query_base = os.path.basename(query)
				logger.debug("sed 's/-/N/g' {query} > {tmpdir}/{query_base}.tmp".format(query=query,query_base=query_base,tmpdir=self.tmpdir))
				self.xmfa_files = []
				os.system("sed 's/-/N/g' {query} > {tmpdir}/{query_base}.tmp".format(query=query,query_base=query_base,tmpdir=self.tmpdir))
				new_query = "{tmpdir}/{query}.tmp".format(query=query_base,tmpdir=self.tmpdir)
				logger.debug("New query: {nq}".format(nq=new_query))
				commands,logs = self.create_mauve_command(new_query,references)
				ret = self.run_mauve(commands,logs)
			if ret != 0:
				return []
			logger.info("Alignments for {query} complete!".format(query=query))
		return self.xmfa_files

	'''Functions'''

	def create_tree(self,SNPS,name,called_snps,save_tree,min_required_hits,strictness=0.7, summary=False):
		'''This function uses ETE3 to color the SNP tree in the database with SNPS found in the reference database
			and outputs a pdf file
		'''
		newickTree = NewickTree(self.database,name,self.outdir,min_required_hits=min_required_hits, strictness=strictness)
		final_snp = newickTree.draw_ete3_tree(SNPS,called_snps,save_tree,summary=summary)
		logger.info("{outdir}/{name}_tree.pdf".format(outdir =self.outdir, name=name))
		return final_snp

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

	def parse_xmfa(XMFA_obj, xmfa_file,results=[]):
		'''Process xmfa file using ParseXMFA object'''
		XMFA_obj.run(xmfa_file)
		results.put(XMFA_obj.get_snps())
		export_results.put(XMFA_obj.get_snp_info())
		called_snps.put(XMFA_obj.get_called_snps())
		return results
	
	def find_snps(self,XMFA_obj,xmfa_file,results=[],export_results=[],called_snps=[]):
		'''Align sequences to references and return SNPs'''
		# execute snp calling
		XMFA_obj.run(xmfa_file)
		#/
		# stream output into Queues (if the complete job output is put in then the queue gets filled when having "large" outputs. It has to be streamed to the queue so that the queue can simultanously be cleared)
		for result in XMFA_obj.get_snps().items(): # returns dictionary
			results.put(result)
		for result in XMFA_obj.get_snp_info(): # returns array
			export_results.put(result)
		for result in XMFA_obj.get_called_snps(): # returns array
			called_snps.put(result)
		#/
		# when finished, put a "stop-signal"
		results.put('XMFA_obj_finish_worker_'+XMFA_obj.reference)
		export_results.put('XMFA_obj_finish_worker_'+XMFA_obj.reference)
		called_snps.put('XMFA_obj_finish_worker_'+XMFA_obj.reference)
		#/

	def find_snps_multiproc(self,xmfa_obj,xmfa_files,export=False):
		'''function to run genomes in paralell'''
		jobs = []
		SNPS = {}
		SNP_info = []
		called_snps = []
		result_queue = Queue()
		export_queue = Queue()
		called_queue = Queue()
		for enum,xmfa_file in enumerate(xmfa_files):
			p = Process(target=self.find_snps, args=(xmfa_obj,xmfa_file ,result_queue,export_queue,called_queue))
			p.start()
			jobs.append(p)
			sleep(0.05) ## A short sleep to make sure all threads do not initiate access to the database file simultanously
		
		## Parse output queues from jobs continuously
		finish_signals = {'SNPS':set(),'SNP_info':set(),'called_snps':set()} # will keep track when all workers are finished
		tic = time() # "tic-toc" clock-sound
		while True:
			# check SNP queue
			if not result_queue.empty():
				tmp_data = result_queue.get()
				if type(tmp_data) == str and tmp_data.find('XMFA_obj_finish_worker_') != -1:
					finish_signals['SNPS'].add(tmp_data.split('_')[-1])
					logger.info("Output queue for 'SNPS' finished for reference {reference}".format(reference=tmp_data.split('_')[-1]))
				else:
					SNPS[tmp_data[0]] = tmp_data[1]
			#/
			# check SNP info queue
			if not export_queue.empty():
				tmp_data = export_queue.get()
				if type(tmp_data) == str and tmp_data.find('XMFA_obj_finish_worker_') != -1:
					finish_signals['SNP_info'].add(tmp_data.split('_')[-1])
					logger.info("Output queue for 'SNP_info' finished for reference {reference}".format(reference=tmp_data.split('_')[-1]))
				else:
					SNP_info.append(tmp_data)
			#/
			# check called snps queue
			if not called_queue.empty():
				tmp_data = called_queue.get()
				if type(tmp_data) == str and tmp_data.find('XMFA_obj_finish_worker_') != -1:
					finish_signals['called_snps'].add(tmp_data.split('_')[-1])
					logger.info("Output queue for 'called_snps' finished for reference {reference}".format(reference=tmp_data.split('_')[-1]))
				else:
					called_snps.append(tmp_data)
			#/
			# check if finished (all xmfa files has returned the "finish worker" in all output queues)
			finish_signals_called = []
			for queue_name,signals in finish_signals.items():
				if len(signals) == len(xmfa_files):
					finish_signals_called.append(queue_name)
			if len(finish_signals) == len(finish_signals_called):
				logger.info("All output queues finished!")
				break
			#/
			# check if we have computed for X amount of time, then break and assume something is wrong
			max_time_seconds = 5*60 # 5*60 => 5 minutes
			toc = time()
			time_spent = toc - tic
			if time_spent > max_time_seconds:
				logger.warning("Maximum time spent reached to capture queue output. This is not expected to happen and might mean that your output data is corrupted.")
				break
			#/
		##/
		## Wait until processes terminated (their output should already have been processed in the queue's)
		for job in jobs:
			job.join()
		##/
		return SNPS,SNP_info,called_snps

	def read_result_dir(self):
		'''Fetch all final SNPs from result file after run to produce summary'''
		snplist = set()
		for root, dirs, files in os.walk(self.outdir, topdown=False):
			for f in files:
				if f.endswith("_snps.txt"):
					with open(os.path.join(root,f)) as snpfile:
						for row in snpfile:
							if row.startswith("SNP path:"):
								snp = row.split(";")[-1].strip()
								snplist.add(snp)
								self.called_genome[snp] = f.rsplit(".fasta",1)[0]
		return snplist

	def print_summary(self):
		'''Create summary file and tree'''
		summarypath = "{outdir}/snps_summary.txt".format(outdir=self.outdir)
		SNPS = self.read_result_dir()
		with open(summarypath,"w") as summaryout:
			for snp in SNPS:
				print("{query}: {SNP}".format(query=self.called_genome[snp], SNP=snp),file=summaryout)

		'''Print summary tree showing all unique SNPs in the final tree'''
		self.create_tree([],"summary",SNPS,True,min_required_hits=self.min_required_hits,summary=True)

	def run(self,database):
		'''Run CanSNPer2'''
		logger.info("Running CanSNPer2 version-{version}".format(version=__version__))
		if len(self.query) > 0:
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
						verbose=self.verbose)  ## Create XMFA object and connect to database
			'''Walk through the list of queries supplied'''
			if not self.skip_mauve: print("Run {n} alignments to references using progressiveMauve".format(n=len(self.query)))
			for q in self.query:			## For each query file_path
				try:
					self.query_name = os.path.basename(q).rsplit(".",1)[0]  ## get name of file and remove ending

					qfile = q.rsplit("/")[-1]   ## Remove path from query name
					if not os.path.exists(q):
						raise FileNotFoundError("Input file: {qfile} was not found!".format(qfile=q))

					outputfile = "{outdir}/{xmfa}_snps.txt".format(outdir=self.outdir,xmfa=self.query_name)
					if os.path.exists(outputfile) and not self.rerun:
						logger.debug("{outputfile} already exits, skip!".format(outputfile=outputfile))
						continue
					logger.info("Running CanSNPer2 on {query}".format(query=qfile))
					if not self.skip_mauve: ### If mauve command was already run before skip step
						logger.info("Run mauve alignments")

					'''For each query fasta align to all CanSNP references the reference folder
						if skip_mauve parameter is True this the align function will only format xmfa file paths
					'''
					xmfa_files = self.align(q)
					logger.debug(xmfa_files)
					if len(xmfa_files) == 0: ## if keep going is set and mauve exits with an error continue to next sequence
						logger.debug("Mauve exited with a non zero exit status, continue with next sample!")
						logger.warning("Mauve error skip {sample}".format(q))
						self.xmfa_files = []
						continue
					'''Parse Mauve XMFA output and find SNPs; returns SNPS (for the visual tree) and SNP_info (text file output)'''
					logger.info("Find SNPs")
					try:
						SNPS,SNP_info,called_snps = self.find_snps_multiproc(xmfa_obj=parse_xmfa_obj,xmfa_files=xmfa_files,export=True)
					except FileNotFoundError:
						logger.warning("One or several xmfa files were not found for {qfile} continue with next file".format(qfile=qfile))
						self.xmfa_files = []
						continue
					'''If file export is requested print the result for each SNP location to file'''
					if self.export:
						outputfile = "{outdir}/{xmfa}_not_called.txt".format(outdir=self.outdir,xmfa=self.query_name)
						outputfile2 = "{outdir}/{xmfa}_snps.txt".format(outdir=self.outdir,xmfa=self.query_name)

						logger.info("Printing SNP info of non called SNPs to {file}".format(file=outputfile))
						self.csnpdict = {}
						'''Print SNPs to tab separated file'''
						with open(outputfile,"w") as snplist_out:
							print("\t".join(["Name","Reference","Pos","Ancestral base","Derived base", "Target base"]),file=snplist_out)
							for snp in SNP_info:
								if snp[0] in called_snps:
									self.csnpdict[snp[0]] = snp
								else:
									print("\t".join(snp),file=snplist_out)

					'''If save tree is requested print tree using ETE3 prints a pdf tree output'''
					SNP = "NA" ## Default message if SNP cannot be confirmed
					final_snp,message,called = self.create_tree(SNPS,self.query_name,called_snps,self.save_tree,min_required_hits=self.min_required_hits,strictness=self.strictness)
					if final_snp:
						SNP = final_snp[1]
						if not final_snp[2][0]:  ## if snp was never confirmed print NA
							SNP = "NA"
						if self.export:
							with open(outputfile2, "w") as called_out:
								if True:
									print("\t".join(["Name","Reference","Pos","Ancestral base","Derived base", "Target base"]),file=called_out)
									for snp in called:
										print("\t".join(self.csnpdict[snp[1]]),file=called_out)
								print("SNP path: {path}".format(path=";".join([snp[1] for snp in called])),file=called_out)
								print("Final SNP: {snp} found/depth: {found}/{depth}".format(snp=SNP,depth=int(final_snp[0]),found=final_snp[2][1]),file=called_out)
						logger.info("Final SNP: {snp} found/depth: {found}/{depth}".format(snp=SNP,depth=int(final_snp[0]),found=final_snp[2][1]))
					else:
						if self.export:
							with open(outputfile2, "a") as called_out:
								print("Final SNP: {snp}".format(snp=SNP), file=called_out)
						logger.info(message)
					if self.summary and SNP != "NA":
						self.summary_set |= set([SNP])
						self.called_genome[SNP] = self.query_name
					if self.export:
						print("{query}: {SNP}".format(query=self.query_name, SNP=SNP))
					'''Clean references to aligned xmfa files between queries if several was supplied'''
					self.xmfa_files = []
				except:
					if not self.keep_going:
						raise CanSNPer2Error("A file did not run correctly exit CanSNPer2 (use --keep_going to continue with next file!)")
					logger.debug("An error occured during processing of {file}".format(file=self.query_name))

		if self.summary:
			self.print_summary()
		'''Finally clean up temporary folder when all alignments and trees has been printed!'''
		if not self.keep_temp and len(self.query) > 0: ## if keep temp is turned on do not remove away alignments also if no input files were given
			self.cleanup()

		logger.info("CanSNPer2 finished successfully, files can be found in {outdir}".format(outdir=self.outdir+"/"))
