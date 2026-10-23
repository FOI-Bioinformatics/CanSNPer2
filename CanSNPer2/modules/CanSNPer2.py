'''
CanSNPer2 module: A toolkit for SNP-typing using NGS data.
Copyright (C) 2019 David Sundell @ FOI bioinformatics group
'''

import os
import shutil
import logging
logger = logging.getLogger(__name__)

## import CanSNPer2 specific modules
from CanSNPer2.modules.ParseXMFA import ParseXMFA
from CanSNPer2.modules.NewickTree import NewickTree
from CanSNPer2.modules.DatabaseConnection import CanSNPdbFunctions
from CanSNPer2 import __version__


## import standard python libraries for subprocess and multiprocess
from subprocess import Popen,PIPE,STDOUT
from multiprocessing import Process, Queue
from time import sleep,time
import signal

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
		self.export = kwargs.get("export", False)  ## export is optional, defaults to False (issue #51)
		self.database = database
		self.min_required_hits = kwargs.get("min_required_hits",3)
		self.strictness = kwargs.get("strictness",0.7)
		self.rerun = kwargs.get("rerun",False)

		'''Create log and tmpdir if they do not exist'''
		self.workdir = kwargs.get("workdir","./")
		if not os.path.exists(self.workdir):
			logger.info("Creating working directory {workdir}".format(workdir=self.workdir))
			os.makedirs(self.workdir)

		'''Resolve the working paths against workdir instead of changing the
			process working directory, which is global state (issue #2)'''
		self.tmpdir = os.path.join(self.workdir, tmpdir)
		if not os.path.exists(self.tmpdir):
			logger.info("Creating tmp directory {tmpdir}".format(tmpdir=self.tmpdir))
			os.makedirs(self.tmpdir)

		self.logdir = os.path.join(self.workdir, kwargs.get("logdir","logs"))
		if not os.path.exists(self.logdir):
			logger.info("Creating logs directory {logdir}".format(logdir=self.logdir))
			os.makedirs(self.logdir)

		self.outdir = os.path.join(self.workdir, kwargs.get("outdir","results"))
		if not os.path.exists(self.outdir):
			logger.info("Creating output directory {outdir}".format(outdir=self.outdir))
			os.makedirs(self.outdir)

		self.refdir = os.path.join(self.workdir, self.refdir)

		'''Fetch other key word arguments'''
		self.skip_mauve = kwargs.get("skip_mauve",False)
		self.save_tree = kwargs.get("save_tree",False)
		self.keep_temp = kwargs.get("keep_temp",False)
		self.keep_going = keep_going

		if kwargs.get("summary",False):
			self.summary_set = set()
			self.called_genome = {}
			self.summary = True
		else:
			self.summary = False

	'''CanSNPer2 get functions'''

	def get_xmfa(self):
		'''Return references to aligned files'''
		return self.xmfa_files

	def get_references(self,selected=False):
		'''references must be available in the refdir, and only the references listed in the
			database are aligned against: a shared reference folder otherwise makes mauve align a
			genome against unrelated species, which wastes hours and fails with a non zero exitcode (issue #124)
		'''
		refs = sorted(ref for ref in os.listdir(self.refdir) if ref.endswith(".fna"))
		if not refs:
			'''without references there is nothing to align against, say so instead of failing later
				as "Mauve error skip" for every sample (issue #124)'''
			raise CanSNPer2Error("No reference .fna files were found in {refdir}, run CanSNPer2-download --database {database} --outdir {refdir} first".format(refdir=self.refdir,database=self.database))
		try:
			db = CanSNPdbFunctions(self.database,verbose=self.verbose)
			db_refs = set()
			for genome,strain in db.query('''SELECT genome,strain FROM snp_references''').fetchall():
				db_refs.update([genome,strain])
		except Exception as e:
			logger.warning("Could not read references from the database, using all fasta files in {refdir}: {err}".format(refdir=self.refdir,err=e))
			return refs
		matched = [ref for ref in refs if ref.rsplit(".",1)[0] in db_refs]
		if not matched:
			logger.warning("None of the {n} reference files in {refdir} matched the references in the database, using all of them!".format(n=len(refs),refdir=self.refdir))
			return refs
		skipped = [ref.rsplit(".",1)[0] for ref in refs if ref not in matched]
		if skipped:
			logger.warning("Skipping {n} reference file(s) in {refdir} that are not in the database: {names}".format(n=len(skipped),refdir=self.refdir,names=", ".join(skipped)))
		return matched

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
		logfiles = []   #open log file container, parallel to processes
		error = 0   #Variable for errors
		logger.info("Starting progressiveMauve on {n} references".format(n=len(commands)))
		for i in range(len(commands)): #Loop through commands,
			command = commands[i]
			logger.debug(command)			## In verbose mode print the actual mauve command
			## stdout/stderr go straight to the log file: a pipe buffer would fill up and
			## block mauve while this loop only polls for exitcodes, hanging the run (issue #68)
			log = open(logs[i], "w")
			p = Popen(command,  stdout=log, stderr=STDOUT)  ##command is an argument list, no shell and no fragile splitting (issue #3)
			logfiles.append(log)
			processes.append(p)
		## SIGTERM (timeout, cluster schedulers, kill) ends the process at once and would leave every
		## progressiveMauve child running unattended, so route it through the cleanup below (issue #151)
		def stop_on_sigterm(signum,frame):
			raise KeyboardInterrupt()
		old_sigterm = None
		try:
			old_sigterm = signal.signal(signal.SIGTERM, stop_on_sigterm)
		except ValueError:
			logger.debug("SIGTERM handling is only possible in the main thread, no child cleanup installed")
		try:
			while processes:
				finished = False
				for idx in range(len(processes)-1, -1, -1): ## Loop through processes, walking backwards so removals are safe
					p = processes[idx]
					exitcode = p.poll()  #get exitcode for each process
					if exitcode is not None: ## if there is no exitcode the program is still running
						finished = True
						## When a process is finished close its log file
						logfiles[idx].close()
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
								logger.error("progressiveMauve failed with exitcode {exitcode} for {log}".format(exitcode=exitcode,log=logs[idx]))
							error += 1
							logger.warning("WARNING progressiveMauve finished with a non zero exitcode: {exitcode}\nThe script will terminate when all processes are finished read {log} for more info".format(log=logs[idx],exitcode=exitcode))
						## Remove process from container
						processes.pop(idx)
						logfiles.pop(idx)
				if not finished: ## the poll loop would otherwise run flat out while mauve works (issue #152)
					sleep(0.5)
		except KeyboardInterrupt:
			'''a cancelled CanSNPer2 left every progressiveMauve child running unattended, one of them
				was still aligning sixteen hours after its parent had been stopped (issue #151)'''
			logger.warning("CanSNPer2 was interrupted, stopping {n} running progressiveMauve processes".format(n=len(processes)))
			for p in processes:
				if p.poll() is None:
					p.terminate()
			for p in processes:
				try:
					p.wait(timeout=10)
				except Exception:
					p.kill()
			for log in logfiles:
				log.close()
			raise
		finally:
			'''the SIGTERM handler belongs to this alignment round only, put the previous one back'''
			if old_sigterm is not None:
				signal.signal(signal.SIGTERM, old_sigterm)
		if retvalue == 6 or retvalue == 11: ## This is done if sequence is not free of gaps
			return retvalue
		elif not self.keep_going:
			if error:  ## Error handling regarding mauve subprocesses, stop script if any of them fails
				logger.error("Error: {errors} progressiveMauve processes did not run correctly check log files for more information".format(errors=error))
				raise MauveError("Error: {errors} progressiveMauve processes did not run correctly check log files for more information".format(errors=error))
		elif error:
			'''keep_going only means that a failing sample is skipped instead of stopping the run,
				a run where every process succeeded must still report success (issue #126)'''
			retvalue=1
		return retvalue

	def stage_for_mauve(self,path,prefix):
		'''progressiveMauve writes a multi megabyte "<input>.sslist" progress file next to every input
			genome, which polluted the reference and query folders and made them required to be writable.
			There is no mauve option that turns this off, so the sequences are symlinked into the tmpdir
			instead and the progress files land there and are removed by cleanup() (issue #148)
		'''
		source = os.path.abspath(path)
		staged = os.path.join(self.tmpdir, "{prefix}_{name}".format(prefix=prefix,name=os.path.basename(source)))
		if not os.path.lexists(staged):
			try:
				os.symlink(source,staged)
			except (OSError,NotImplementedError,AttributeError):
				'''a tmpdir on a filesystem without symlinks, fall back to a copy'''
				logger.debug("Could not symlink {source} into {tmpdir}, copying the genome instead".format(source=source,tmpdir=self.tmpdir))
				shutil.copyfile(source,staged)
		return staged

	def create_mauve_command(self,query,references=[]):
		'''Mauve commands'''
		commands =[]	# store execute command
		logs = []	   # store log filepath
		if len(references) == 0:	## If specific references are not given fetch references from the reference folder
			references = self.get_references()
		query_file = self.stage_for_mauve(query,"query")
		for ref in references:	  ## For each reference in the reference folder align to query
			ref_name = ref.rsplit(".",1)[0] ## remove file ending
			xmfa_output = "{tmpdir}/{ref}_{target}.xmfa".format(tmpdir=self.tmpdir.rstrip("/"),ref=ref_name,target=self.query_name)
			ref_file = self.stage_for_mauve("{refdir}/{ref}".format(refdir=self.refdir, ref=ref),"ref")
			log_file = "{logdir}/{ref}_{target}.mauve.log".format(logdir=self.logdir,ref=ref_name,target=self.query_name)

			'''Create run command for mauve as an argument list, so paths with spaces stay single arguments (issue #3)'''
			command = ["{mauve_path}progressiveMauve".format(mauve_path=self.mauve_path),
							"--output", xmfa_output,
							ref_file,
							query_file]
			commands.append(command)			## Mauve command
			logs.append(log_file)			   ## Store log files for each alignment
			self.xmfa_files.append(xmfa_output) ## Store the path to xmfa files as they will be used later
		return commands,logs

	def align(self, query, references=[]):
		'''Align sequences and run mauve as subprocess'''
		'''the name is set here, from the query actually being aligned, so that create_mauve_command
			never depends on run() having set it before and the gap retry below keeps the original
			sample name for its xmfa files (issue #160)'''
		self.query_name = os.path.basename(query).rsplit(".",1)[0]
		commands,logs = self.create_mauve_command(query,references)
		if not self.skip_mauve: ### If mauve command was already run before don´t run mauve return xmfa paths
			ret = self.run_mauve(commands,logs)
			if int(ret) in (6,11):
				'''exitcode 11 is often mauve not handling dashes(-) in the sequence, and run_mauve already
					reports exitcode 6 as "not free of gaps, replace gaps with N and retry", so retry the same
					way for both instead of only for 11 (issue #127)'''
				logger.warning("Mauve ran into an error with sequence {seq} may contain a dash, replace with N characters and retry mauve".format(seq=query))
				query_base = os.path.basename(query)
				new_query = os.path.join(self.tmpdir, query_base + ".tmp")
				logger.debug("Replacing all '-' with 'N' in {query}, writing copy to {tmp}".format(query=query,tmp=new_query))
				self.xmfa_files = []
				with open(query, "r") as seq_in, open(new_query, "w") as seq_out:
					for line in seq_in:
						seq_out.write(line.replace("-", "N"))
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

	def parse_xmfa(self, XMFA_obj, xmfa_file, results, export_results, called_snps):
		'''Process xmfa file using ParseXMFA object, streaming the three result types into their queues'''
		XMFA_obj.run(xmfa_file)
		results.put(XMFA_obj.get_snps())
		export_results.put(XMFA_obj.get_snp_info())
		called_snps.put(XMFA_obj.get_called_snps())
		return results
	
	def find_snps(self,XMFA_obj,xmfa_file,results=None,export_results=None,called_snps=None):
		'''Align sequences to references and return SNPs'''
		# default to None instead of mutable list defaults (issue #6); fresh queues when not supplied
		if results is None: results = Queue()
		if export_results is None: export_results = Queue()
		if called_snps is None: called_snps = Queue()
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
								'''the result files are named {query}_snps.txt, strip that suffix to get the query name (issue #123)'''
								self.called_genome[snp] = f[:-len("_snps.txt")]
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
						logger.warning("Mauve error skip {sample}".format(sample=q))
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
					tree_result = self.create_tree(SNPS,self.query_name,called_snps,self.save_tree,min_required_hits=self.min_required_hits,strictness=self.strictness)
					if not isinstance(tree_result, tuple):  ## create_tree returns None when it only draws a summary tree, no SNP call to consume
						tree_result = (False,"No SNP call was produced.",[])
					final_snp,message,called = tree_result
					if final_snp:
						SNP = final_snp[1]
						if not final_snp[2][0]:  ## if snp was never confirmed print NA
							SNP = "NA"
						if self.export:
							with open(outputfile2, "w") as called_out:
								if called:
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
				except Exception as e:
					'''log the real failure before wrapping it, otherwise the cause is invisible (issue #121)'''
					logger.debug("Traceback for {file}".format(file=self.query_name),exc_info=True)
					logger.error("Error while processing {file}: {err}".format(file=self.query_name,err=e))
					if not self.keep_going:
						raise CanSNPer2Error("A file did not run correctly exit CanSNPer2 (use --keep_going to continue with next file!)") from e
					logger.warning("An error occured during processing of {file}, continuing!".format(file=self.query_name))

		if self.summary:
			self.print_summary()
		'''Finally clean up temporary folder when all alignments and trees has been printed!'''
		if not self.keep_temp and len(self.query) > 0: ## if keep temp is turned on do not remove away alignments also if no input files were given
			self.cleanup()

		logger.info("CanSNPer2 finished successfully, files can be found in {outdir}".format(outdir=self.outdir.rstrip("/")+"/")) ## a trailing slash was already there would print dir// (issue #76)
