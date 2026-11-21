'''
CanSNPer2 module: A toolkit for SNP-typing using NGS data.
Copyright (C) 2019 David Sundell @ FOI bioinformatics group
'''

import os
import shutil
import hashlib
import tempfile
import logging
logger = logging.getLogger(__name__)

## import CanSNPer2 specific modules
from CanSNPer2.modules.Aligners import Aligner, get_aligner
from CanSNPer2.modules.NewickTree import NewickTree
from CanSNPer2.modules.DatabaseConnection import CanSNPdbFunctions
from CanSNPer2.modules.ReferenceCheck import check_reference, duplicate_names, reference_snps
from CanSNPer2 import __version__


## import standard python libraries for subprocess and multiprocess
from subprocess import Popen,PIPE,STDOUT
import multiprocessing
from queue import Empty
from time import sleep,time
import signal

# the workers are forked and inherit the parent's memory, which is what the fork start method
# gives. Python 3.14 made forkserver the default on linux, where the worker target has to be
# pickled and "cannot pickle sqlite3.Connection" killed every run, so the context is requested
# explicitly (issue #186). Queues must come from the same context as the processes. The database
# connection is opened inside the worker now, from the path ParseXMFA carries, so no cursor is
# shared between processes (issue #49).
_MP_CTX = multiprocessing.get_context("fork")

'''column header of the per sample SNP tables, both <sample>_snps.txt and
	<sample>_not_called.txt. Sequence is the reference sequence the SNP lies on and Pos the
	1-based position on that sequence, as the database stores them; the last column carries
	the "edge" flag for SNPs that sit near an alignment edge (issue #48)'''
SNP_FILE_COLUMNS = ["Name","Reference","Sequence","Pos","QuerySeq","QueryPosition","Ancestral base","Derived base","Target base","Edge flag"]
'''index of the edge flag in a result row'''
EDGE_COLUMN = SNP_FILE_COLUMNS.index("Edge flag")

'''file that marks a private run directory created by CanSNPer2 below --tmpdir; cleanup()
	only removes a directory it created itself, and the command line uses a marked
	directory as it is so a kept run can be reused with --skip_alignment'''
RUN_MARKER = ".cansnper2_run"

class Error(Exception):
	"""docstring for Error"""
	def __init__(self, value):
		self.value = value
	def __str__(self):
		return repr(self.value)

class AlignmentError(Error):
	"""an alignment process failed"""
	pass

'''the name from the time progressiveMauve was the only aligner'''
MauveError = AlignmentError

class CanSNPer2Error(Error):
	"""docstring for CanSNPer2Error"""
	pass

class CanSNPer2(object):
	"""docstring for CanSNPer2"""
	def __init__(self, query, mauve_path="",tmpdir=".tmp", refdir="references",database="CanSNPer.fdb", verbose=False,keep_going=False,**kwargs):
		super(CanSNPer2, self).__init__()
		self.query = query
		self.verbose = verbose
		self.refdir = refdir
		'''the aligner is a name from Aligners.ALIGNERS or an Aligner object; aligner_path is
			the program or the folder that holds it (mauve_path is its old name)'''
		aligner = kwargs.get("aligner")
		if isinstance(aligner, Aligner):
			self.aligner = aligner
		else:
			self.aligner = get_aligner(aligner, kwargs.get("aligner_path") or mauve_path)
		self.mauve_path = mauve_path
		self.alignment_files = []      ## alignment of the current query against each reference
		self.alignment_references = [] ## reference genome name of each file in alignment_files
		self.export = kwargs.get("export", False)  ## export is optional, defaults to False (issue #51)
		self.database = database
		self.min_required_hits = kwargs.get("min_required_hits",3)
		self.strictness = kwargs.get("strictness",0.7)
		self.rerun = kwargs.get("rerun",False)

		'''what each sample ended up as, so the command can tell a typed run from an
			all NA run instead of always claiming success (issue #191)'''
		self.typed = {}
		self.not_typed = {}
		self.skipped = {}
		'''per sample outcome for the summary.html dashboard: severity, reasons and the SNP
			states the embedded tree paints when the sample is selected'''
		self.sample_reports = {}
		'''the SNP path per sample, so the summary can print it for every sample (issue #24)'''
		self.result_paths = {}
		'''samples whose result file from an earlier run was used again; only these are read
			back for the summary: a sample skipped because its alignment failed used to show
			the call of whatever result file an earlier run had left under its name'''
		self.reused = set()
		self.query_sha256 = None
		'''how long the queue loop may sit without a single item before it calls the run
			stuck; a worker that keeps delivering is never cut off (issue #190)'''
		self.queue_idle_timeout = kwargs.get("queue_idle_timeout", 5*60)
		'''how long a worker may stay alive after the collection loop ended before it is
			stopped: joining forever is what turned a stalled worker into a run that never
			returned (issue #41)'''
		self.worker_join_timeout = kwargs.get("worker_join_timeout", 2*60)
		'''SNPs closer to an alignment edge than this are not evaluated (issue #175)'''
		self.mask = int(kwargs.get("mask",0) or 0)
		'''how long one round of alignments (the query against all references, in parallel) may
			take before the aligners are stopped, in minutes; 0 means no limit. progressiveMauve
			never returns on some inputs, which used to hang the run'''
		self.alignment_timeout = float(kwargs.get("alignment_timeout",60) or 0)*60
		self.checked_references = {} ## reference file -> problem or None, checked once per run

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
		'''An aligning run works in its own run_* directory below the shared tmpdir. All runs
			used to write their staged genomes and alignments straight into the shared folder
			(by default /tmp/CanSNPer2), and cleanup() removed every file in it, so two runs on
			the same machine deleted each other's files. With --skip_alignment the given folder is
			where existing alignments are read from, it is used as it is and never removed.'''
		self.own_tmpdir = False
		'''skip_mauve is the old name of skip_alignment'''
		self.skip_alignment = bool(kwargs.get("skip_alignment",False) or kwargs.get("skip_mauve",False))
		if not self.skip_alignment:
			self.tmpdir = tempfile.mkdtemp(prefix="run_", dir=self.tmpdir)
			with open(os.path.join(self.tmpdir, RUN_MARKER), "w"):
				pass
			self.own_tmpdir = True
			logger.info("Temporary files of this run are written to {tmpdir}".format(tmpdir=self.tmpdir))

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
		self.save_tree = kwargs.get("save_tree",False)
		'''tree annotations for the html views: a file path (load_annotations parses it) or
			an already built list of {node, kind, rank, label} dicts'''
		annotations = kwargs.get("annotations",None)
		if isinstance(annotations,str):
			from CanSNPer2.modules.NewickTree import load_annotations
			annotations = load_annotations(annotations)
		self.annotations = annotations or []
		self.keep_temp = kwargs.get("keep_temp",False)
		self.keep_going = keep_going

		if kwargs.get("summary",False):
			self.summary = True
		else:
			self.summary = False

	'''old attribute names, from the time progressiveMauve was the only aligner'''
	@property
	def skip_mauve(self):
		return self.skip_alignment

	@property
	def xmfa_files(self):
		return self.alignment_files

	@xmfa_files.setter
	def xmfa_files(self, files):
		self.alignment_files = files

	'''CanSNPer2 get functions'''

	def get_alignments(self):
		'''Return the paths of the alignment files of the current query'''
		return self.alignment_files

	get_xmfa = get_alignments

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

	'''Alignment functions'''

	def run_alignments(self, commands, logs, names=None):
		'''Run the alignment commands as parallel subprocesses

			To run subprocesses securely this function uses Popen with an argument list; stderr
			is passed to the stdout of the process, which goes to one log file per alignment.

			The number of processes is the number of references in the database; there is no
			further limit, the OS distributes the processes.

			names holds (reference, query) for each command, used in the messages. Returns 0
			when every process succeeded, a code of self.aligner.retry_codes when a process
			failed in a way the aligner retries, and 1 when a process failed with --keep_going;
			without --keep_going a failure raises AlignmentError.
		'''
		program = self.aligner.name
		if names is None:
			names = [self._names_from_command(command) for command in commands]
		retvalue=0
		processes = []  #process container
		logfiles = []   #open log file container, parallel to processes
		error = 0   #Variable for errors
		logger.info("Starting {program} on {n} references".format(program=program,n=len(commands)))
		running_names = []      #reference/query name of each running alignment, parallel to processes
		logpaths = []   #log path of each running alignment, parallel to processes
		for i in range(len(commands)): #Loop through commands,
			command = commands[i]
			logger.debug(command)			## In verbose mode print the actual command
			## stdout/stderr go straight to the log file: a pipe buffer would fill up and
			## block the aligner while this loop only polls for exitcodes, hanging the run (issue #68)
			log = open(logs[i], "w")
			p = Popen(command,  stdout=log, stderr=STDOUT)  ##command is an argument list, no shell and no fragile splitting (issue #3)
			logfiles.append(log)
			processes.append(p)
			running_names.append(names[i])
			logpaths.append(logs[i])
		## SIGTERM (timeout, cluster schedulers, kill) ends the process at once and would leave every
		## aligner child running unattended, so route it through the cleanup below (issue #151)
		def stop_on_sigterm(signum,frame):
			raise KeyboardInterrupt()
		old_sigterm = None
		try:
			old_sigterm = signal.signal(signal.SIGTERM, stop_on_sigterm)
		except ValueError:
			logger.debug("SIGTERM handling is only possible in the main thread, no child cleanup installed")
		started = time()
		try:
			while processes:
				if self.alignment_timeout and time() - started > self.alignment_timeout:
					'''an aligner that never returns used to hang the whole run'''
					for idx in range(len(processes)):
						ref_name, query_name = running_names[idx]
						logger.error("{program} did not finish within {minutes:g} minutes for genome {query} against reference {ref}, it was stopped (--alignment_timeout)".format(
								program=program, minutes=self.alignment_timeout/60, query=query_name, ref=ref_name))
						processes[idx].terminate()
					for idx, p in enumerate(processes):
						try:
							p.wait(timeout=10)
						except Exception:
							p.kill()
						logfiles[idx].close()
					error += len(processes)
					processes, logfiles, running_names, logpaths = [], [], [], []
					break
				finished = False
				for idx in range(len(processes)-1, -1, -1): ## Loop through processes, walking backwards so removals are safe
					p = processes[idx]
					exitcode = p.poll()  #get exitcode for each process
					if exitcode is not None: ## if there is no exitcode the program is still running
						finished = True
						## When a process is finished close its log file
						logfiles[idx].close()
						## the containers are popped as processes finish, so the running position
						## in these containers is the index for the names and logs, not the one
						## in the original commands list (issue #12)
						ref_name, query_name = running_names[idx]
						status, code = self.aligner.exit_status(exitcode)
						if status == "retry":
							logger.warning("{message} for genome {query} against reference {ref}".format(message=self.aligner.retry_message(code),query=query_name,ref=ref_name))
							retvalue=code
						elif status == "error":
							if not self.keep_going:
								logger.error("{program} failed with exitcode {exitcode} for {log}".format(program=program,exitcode=exitcode,log=logpaths[idx]))
							error += 1
							logger.warning("WARNING {program} finished with a non zero exitcode: {exitcode} for genome {query} against reference {ref}\nThe script will terminate when all processes are finished read {log} for more info".format(program=program,log=logpaths[idx],exitcode=exitcode,query=query_name,ref=ref_name))
						## Remove process from container
						processes.pop(idx)
						logfiles.pop(idx)
						running_names.pop(idx)
						logpaths.pop(idx)
				if not finished: ## the poll loop would otherwise run flat out while the aligner works (issue #152)
					sleep(0.5)
		except KeyboardInterrupt:
			'''a cancelled CanSNPer2 left every aligner child running unattended, one of them
				was still aligning sixteen hours after its parent had been stopped (issue #151)'''
			logger.warning("CanSNPer2 was interrupted, stopping {n} running {program} processes".format(n=len(processes),program=program))
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
		if retvalue in self.aligner.retry_codes:
			return retvalue
		elif not self.keep_going:
			if error:  ## stop the run if any of the alignments failed
				logger.error("Error: {errors} {program} processes did not run correctly check log files for more information".format(errors=error,program=program))
				raise AlignmentError("Error: {errors} {program} processes did not run correctly check log files for more information".format(errors=error,program=program))
		elif error:
			'''keep_going only means that a failing sample is skipped instead of stopping the run,
				a run where every process succeeded must still report success (issue #126)'''
			retvalue=1
		return retvalue

	def _names_from_command(self, command):
		'''(reference, query) from a progressiveMauve command [program, --output, xmfa, ref_file,
			query_file]: the staged input names carry "ref_"/"query_" prefixes, and a pre-cleaned
			copy (issues #26, #28) adds a "cleaned_<kind>_" layer in between. Only used when the
			caller gives no names (issue #12)'''
		def staged_name(arg,kind):
			base = os.path.basename(arg).rsplit(".",1)[0]
			for prefix in ("{0}_".format(kind), "cleaned_", "{0}_".format(kind)):
				if base.startswith(prefix):
					base = base[len(prefix):]
			return base
		try:
			return staged_name(command[3],"ref"), staged_name(command[4],"query")
		except IndexError:
			return "unknown reference", "unknown genome"

	def run_mauve(self, commands, logs):
		'''the name from the time progressiveMauve was the only aligner'''
		return self.run_alignments(commands, logs)

	def stage_input(self,path,prefix):
		'''link an input genome into the run directory, for aligners that write files next to
			their inputs (progressiveMauve, issue #148); other aligners read the file where it is'''
		if not self.aligner.stage_inputs:
			return path
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

	stage_for_mauve = stage_input

	def clean_input_fasta(self,path,kind):
		'''Replace every sequence character that is not A/C/G/T (or their lowercase) with N before
			anything is aligned: gaps ('-', '.', '*') made progressiveMauve fail with exitcode 6/11
			(the gap retry used to replace dashes only after a failure, and could not see dashes it
			could not "see", issue #28), and the ambiguous IUPAC codes (R,Y,M,K,S,W,B,D,H,V) can never
			match the target base of a SNP, so replacing them with N is the honest call (issue #26).
			Every replacement keeps the exact position of the base, so database SNP positions still
			match, and an already clean file is returned unchanged.
			Whitespace is not sequence: blank lines are dropped (progressiveMauve never returns on a
			file that starts with one) and spaces or tabs inside a sequence line are removed. They
			used to become N, which added bases and shifted every position after them.
		'''
		keep = set("ACGTacgtn")
		dirty = False
		lines = []
		names = []
		with open(path) as seq_in:
			for line in seq_in:
				if line.startswith(">"):
					lines.append(line)
					words = line[1:].split()
					names.append(words[0] if words else "")
					continue
				body = line.rstrip("\n")
				bases = "".join(body.split())
				if bases != body or not bases: ## an empty line is unchanged by the split but is dropped
					dirty = True
				if not bases:
					continue
				cleaned = "".join(ch if ch in keep or ch in "Nn" else "N" for ch in bases)
				if cleaned != bases:
					dirty = True
				lines.append(cleaned + "\n")
		'''sequences are identified by the first word of their header (aligners, result files);
			a query with two sequences of the same name gave bases and positions of the wrong
			contig. References are checked once per run in check_reference_file.'''
		if kind == "query":
			problem = duplicate_names(names)
			if problem:
				raise CanSNPer2Error("query {path} has {problem}; rename them, for example with seqkit rename".format(path=path, problem=problem))
		if not dirty:
			return path
		cleaned_path = os.path.join(self.tmpdir, "cleaned_{kind}_{name}".format(kind=kind,name=os.path.basename(path)))
		with open(cleaned_path,"w") as seq_out:
			seq_out.write("".join(lines))
		logger.info("Replaced characters that are not ACGT with N and removed whitespace in {name} before alignment (gaps, ambiguous nucleotides and blank lines, issues #26, #28)".format(name=os.path.basename(path)))
		return cleaned_path

	def check_reference_file(self,path,name):
		'''The reference must have every sequence the database places SNPs on, each at least as
			long as its largest SNP position, and the ancestral or derived base at its SNPs. A
			truncated, renamed or wrong reference otherwise gave missing or wrong genotypes
			without an error. Checked once per run.'''
		if path not in self.checked_references:
			if not hasattr(self, "reference_snps"):
				try:
					'''sqlite creates an empty file for a path that does not exist'''
					if not self.database or not os.path.isfile(self.database):
						raise FileNotFoundError("{database} is not a file".format(database=self.database))
					self.reference_snps = reference_snps(CanSNPdbFunctions(self.database,verbose=self.verbose))
				except Exception as e:
					logger.warning("Could not read the SNP positions from the database, the reference files are not checked for sequences, length and bases: {e}".format(e=e))
					self.reference_snps = {}
			self.checked_references[path] = check_reference(path, name, self.reference_snps.get(name))
		if self.checked_references[path]:
			raise CanSNPer2Error(self.checked_references[path])

	def create_alignment_commands(self,query,references=[]):
		'''One alignment command per reference. Returns the commands, their log files and
			(reference, query) names; the alignment file and reference name of each command
			are stored in alignment_files and alignment_references.'''
		commands =[]	# store execute command
		logs = []	   # store log filepath
		names = []     # (reference, query) of each command
		if len(references) == 0:	## If specific references are not given fetch references from the reference folder
			references = self.get_references()
		query_file = self.stage_input(self.clean_input_fasta(query,"query"),"query")
		for ref in references:	  ## For each reference in the reference folder align to query
			ref_name = ref.rsplit(".",1)[0] ## remove file ending
			output = self.aligner.alignment_path(self.tmpdir.rstrip("/"), ref_name, self.query_name)
			ref_path = "{refdir}/{ref}".format(refdir=self.refdir, ref=ref)
			'''with --skip_alignment the reference is checked when it is at hand: a renamed
				reference gave NA calls and nothing in the log. Alignments can be read again
				without the reference folder, the check is then left out.'''
			if not self.skip_alignment or os.path.isfile(ref_path):
				self.check_reference_file(ref_path, ref_name)
			ref_file = self.stage_input(self.clean_input_fasta("{refdir}/{ref}".format(refdir=self.refdir, ref=ref),"ref"),"ref")
			commands.append(self.aligner.command(ref_file, query_file, output))
			if not self.skip_alignment:
				self.aligner.prepare(ref_file, query_file, output)
			logs.append(self.aligner.log_path(self.logdir, ref_name, self.query_name))
			names.append((ref_name, self.query_name))
			self.alignment_files.append(output)        ## the alignment files are read later
			self.alignment_references.append(ref_name)
		return commands,logs,names

	def create_mauve_command(self,query,references=[]):
		'''the name from the time progressiveMauve was the only aligner'''
		commands,logs,_names = self.create_alignment_commands(query,references)
		return commands,logs

	def align(self, query, references=[]):
		'''Align the query against the references with the selected aligner and return the
			alignment files, an empty list when the alignments failed with --keep_going'''
		'''the name is set here, from the query actually being aligned, so that the commands
			never depend on run() having set it before and the retry below keeps the original
			sample name for its alignment files (issue #160)'''
		self.query_name = os.path.basename(query).rsplit(".",1)[0]
		self.alignment_files = []
		self.alignment_references = []
		commands,logs,names = self.create_alignment_commands(query,references)
		if not self.skip_alignment: ### alignments made before (--skip_alignment) are only read
			ret = self.run_alignments(commands,logs,names)
			if ret in self.aligner.retry_codes:
				new_query = self.aligner.prepare_retry(query, self.tmpdir)
				if new_query:
					logger.warning("{program} ran into an error with sequence {seq}, retrying with a cleaned copy".format(program=self.aligner.name,seq=query))
					self.alignment_files = []
					self.alignment_references = []
					commands,logs,names = self.create_alignment_commands(new_query,references)
					ret = self.run_alignments(commands,logs,names)
			if ret != 0:
				return []
			logger.info("Alignments for {query} complete!".format(query=query))
		return self.alignment_files

	'''Functions'''

	def create_tree(self,SNPS,name,called_snps,save_tree,min_required_hits,strictness=0.7, summary=False):
		'''This function colors the SNP tree of the database with the SNPS found in the query and
			writes the self-contained html viewer of it (NewickTree.draw_tree)
		'''
		newickTree = NewickTree(self.database,name,self.outdir,min_required_hits=min_required_hits, strictness=strictness)
		final_snp = newickTree.draw_tree(SNPS,called_snps,save_tree,summary=summary,annotations=self.annotations)
		logger.info("{outdir}/{name}_tree.html".format(outdir =self.outdir, name=name))
		return final_snp

	def read_query_textfile_input(self,query_file):
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
			raise CanSNPer2Error("The query input file {query_file} contained no filepaths".format(query_file=query_file))
		return query

	def cleanup(self):
		'''Remove the run directory of this run, and nothing else: other runs and other
			files below the shared tmpdir are left untouched. A folder this run did not
			create (the alignments read with --skip_alignment) is never removed.'''
		if not self.own_tmpdir:
			logger.debug("Not removing {tmpdir}, it was not created by this run".format(tmpdir=self.tmpdir))
			return
		logger.info("Clean up temporary files in {tmpdir}".format(tmpdir=self.tmpdir))
		shutil.rmtree(self.tmpdir, ignore_errors=True)
		logger.info("Done!")
		return

	def find_snps(self,XMFA_obj,xmfa_file,results=None,export_results=None,called_snps=None,reference=None):
		'''Read the SNPs of one alignment with the reader XMFA_obj (an AlignmentReader) and
			stream them into the queues. reference is the reference genome name of the
			alignment; without it the reader takes it from the file name.'''
		# default to None instead of mutable list defaults (issue #6); fresh queues when not supplied
		if results is None: results = _MP_CTX.Queue()
		if export_results is None: export_results = _MP_CTX.Queue()
		if called_snps is None: called_snps = _MP_CTX.Queue()
		# execute snp calling
		'''a worker that dies here never sends its finish signal, so the queue loop waited
			for the five minute time out and the sample was typed from whatever had
			arrived; the failure is reported and signalled at once (issue #179).
			Nothing of the half-read alignment is passed on and the failure is sent as its
			own signal, because the sample must not be typed from a parse that stopped
			midway (issue #50)'''
		try:
			if reference:
				XMFA_obj.run(xmfa_file, reference=reference)
			else:
				XMFA_obj.run(xmfa_file)
		except Exception as e:
			reference = getattr(XMFA_obj,"reference",os.path.basename(xmfa_file).split("_")[0])
			logger.error("SNP calling failed for {alignment}: {error}".format(alignment=xmfa_file,error=e))
			for queue in (results,export_results,called_snps):
				queue.put('XMFA_obj_failed_worker_{reference}: {error}'.format(reference=reference,error=e))
			return
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

	def find_snps_multiproc(self,xmfa_obj,xmfa_files,export=False,references=None):
		'''Read the alignments in parallel, one process per alignment file. xmfa_obj is the
			AlignmentReader of the aligner, references the reference genome name of each file
			(parallel to xmfa_files); without it the names are taken from the file names.'''
		'''the alignments are checked before the workers start: ParseXMFA.read_xmfa swallows a
			FileNotFoundError and simply returns no SNPs, so --skip_alignment against a wrong tmpdir
			(or alignments that were moved away) typed every sample NA and exited 0 (issue #176)'''
		missing = [xmfa for xmfa in xmfa_files if not os.path.isfile(xmfa)]
		if missing:
			raise FileNotFoundError("{n} of {m} alignment files were not found in {tmpdir}, for example {first}".format(n=len(missing),m=len(xmfa_files),tmpdir=self.tmpdir,first=missing[0]))
		if references is not None and len(references) != len(xmfa_files):
			raise ValueError("{n} reference names were given for {m} alignment files".format(n=len(references),m=len(xmfa_files)))
		jobs = []
		SNPS = {}
		SNP_info = []
		called_snps = []
		result_queue = _MP_CTX.Queue()
		export_queue = _MP_CTX.Queue()
		called_queue = _MP_CTX.Queue()
		for enum,xmfa_file in enumerate(xmfa_files):
			reference = references[enum] if references is not None else None
			p = _MP_CTX.Process(target=self.find_snps, args=(xmfa_obj,xmfa_file ,result_queue,export_queue,called_queue,reference))
			p.start()
			jobs.append(p)
			sleep(0.05) ## A short sleep to make sure all threads do not initiate access to the database file simultanously
		
		## Parse output queues from jobs continuously
		finish_signals = {'SNPS':0,'SNP_info':0,'called_snps':0} # how many workers sent their finish signal per queue
		failed_reports = [] ## workers that reported a failed SNP calling, they never count as finished (issue #50)
		idle_timeout_hit = False ## the loop ended because the queues went silent, see the join below (issue #51)
		'''the timeout used to be total wall clock time: a run with more alignments than can
			be parsed in five minutes stopped collecting on its own and typed the sample from
			partial data. It is an idle timeout now, a working worker may take as long as it
			needs, and only a queue that stays silent while no process is alive ends the loop
			(the loop also used to spin at full core between the empty queue checks)'''
		max_idle_seconds = self.queue_idle_timeout ## 5 minutes without a single queue item, by default
		last_activity = time() # "tic-toc" clock-sound
		def drain(queue, batch=1000):
			'''take everything the queue holds right now instead of one item per round: a
				custom database with thousands of SNPs per reference sends tens of thousands of
				items per alignment through these queues, and picking a single item per check
				let the queue pipes fill up until the workers blocked in put() (issue #41)'''
			items = []
			while len(items) < batch:
				try:
					items.append(queue.get_nowait())
				except Empty:
					break
			return items
		while True:
			queues_had_data = False
			# check SNP queue
			for tmp_data in drain(result_queue):
				queues_had_data = True
				if type(tmp_data) == str and tmp_data.find('XMFA_obj_failed_worker_') != -1:
					failed_reports.append(tmp_data)
				elif type(tmp_data) == str and tmp_data.find('XMFA_obj_finish_worker_') != -1:
					finish_signals['SNPS'] += 1
					logger.info("Output queue for 'SNPS' finished for reference {reference}".format(reference=tmp_data.split('_')[-1]))
				else:
					SNPS[tmp_data[0]] = tmp_data[1]
			#/
			# check SNP info queue
			for tmp_data in drain(export_queue):
				queues_had_data = True
				if type(tmp_data) == str and tmp_data.find('XMFA_obj_failed_worker_') != -1:
					failed_reports.append(tmp_data)
				elif type(tmp_data) == str and tmp_data.find('XMFA_obj_finish_worker_') != -1:
					finish_signals['SNP_info'] += 1
					logger.info("Output queue for 'SNP_info' finished for reference {reference}".format(reference=tmp_data.split('_')[-1]))
				else:
					SNP_info.append(tmp_data)
			#/
			# check called snps queue
			for tmp_data in drain(called_queue):
				queues_had_data = True
				if type(tmp_data) == str and tmp_data.find('XMFA_obj_failed_worker_') != -1:
					failed_reports.append(tmp_data)
				elif type(tmp_data) == str and tmp_data.find('XMFA_obj_finish_worker_') != -1:
					finish_signals['called_snps'] += 1
					logger.info("Output queue for 'called_snps' finished for reference {reference}".format(reference=tmp_data.split('_')[-1]))
				else:
					called_snps.append(tmp_data)
			#/
			# check if finished (all xmfa files has returned the "finish worker" in all output queues)
			finish_signals_called = []
			for queue_name,signals in finish_signals.items():
				if signals >= len(xmfa_files):
					finish_signals_called.append(queue_name)
			if len(finish_signals) == len(finish_signals_called):
				logger.info("All output queues finished!")
				break
			#/
			if queues_had_data: ## still delivering, the clock restarts (issue #190)
				last_activity = time()
				continue
			#/
			'''nothing arrived this round: keep waiting while any worker is alive, and stop
				either when they are all gone or when the queues have been silent too long'''
			alive = [job for job in jobs if job.is_alive()]
			if not alive:
				break
			if time() - last_activity > max_idle_seconds:
				logger.error("No queue output for {seconds} seconds while {n} SNP calling process(es) were still running".format(seconds=max_idle_seconds,n=len(alive)))
				idle_timeout_hit = True
				break
			sleep(0.1)
			#/
		##/
		## Wait until processes terminated (their output should already have been processed in the queue's).
		## The wait is bounded: joining forever is what made the whole command hang when one SNP calling
		## worker stalled, which is what "software hang with a large custom SNP database" looked like
		## (issue #41). Stuck workers are stopped here, and the incomplete result check below turns the
		## run into a reported failure instead of a silent freeze.
		## Once the idle timeout had already called the run stuck, waiting the full join timeout on top
		## of the idle wait only made one stalled sample cost seven minutes before --keep_going moved on,
		## so the grace period after a silent queue is short (issue #51).
		join_grace = 10 if idle_timeout_hit else self.worker_join_timeout
		for job in jobs:
			job.join(timeout=join_grace)
		stalled = [job for job in jobs if job.is_alive()]
		if stalled:
			logger.error("Stopping {n} SNP calling process(es) that were still running {seconds} seconds after the queues went silent".format(
					n=len(stalled), seconds=join_grace))
			for job in stalled:
				job.terminate()
			for job in stalled:
				job.join(timeout=10)
				if job.is_alive():
					job.kill()
		##/
		'''a worker whose parse failed mid-alignment used to send normal finish signals, so the
			sample was typed from a half-read alignment; the failure signal is its own report now
			and the run refuses the sample instead (issue #50)'''
		if failed_reports:
			reports = sorted(set(failed_reports))
			raise CanSNPer2Error("{n} of {m} alignments failed during SNP calling ({reports}), refusing to type from a half-read alignment".format(
					n=len(reports), m=len(xmfa_files), reports="; ".join(reports)))
		##/
		'''a worker that died before sending its finish signals used to leave the sample typed
			from whatever had arrived; an incomplete alignment set is an error now (issue #190)'''
		answered = min(finish_signals.values())
		if answered < len(xmfa_files):
			raise CanSNPer2Error("Only {n} of {m} alignments reported their SNPs ({files}), refusing to type from partial data".format(
					n=answered, m=len(xmfa_files), files=", ".join([os.path.basename(xmfa) for xmfa in xmfa_files[:3]])))
		return SNPS,SNP_info,called_snps

	def sample_name(self, path):
		'''the sample name of a query file: its file name without the last extension'''
		return os.path.basename(path).rsplit(".",1)[0]

	def check_sample_names(self, queries):
		'''Every result file, alignment and summary row is named after the sample, so two
			different query files with the same file name (a/x.fasta and b/x.fasta) were typed as
			one sample: the second was skipped as already done and reported with the call of the
			first. Different files with the same name stop the run before anything is aligned;
			the same file given twice is typed once.'''
		unique, seen, clashes = [], {}, {}
		for q in queries:
			name = self.sample_name(q)
			real = os.path.realpath(q)
			if name in seen:
				if seen[name] == real:
					logger.warning("{q} is given more than once, it is typed once".format(q=q))
					continue
				clashes.setdefault(name, [seen[name]]).append(real)
				continue
			seen[name] = real
			unique.append(q)
		if clashes:
			raise CanSNPer2Error("Samples are named after their file name, which must be unique: {clashes}; rename one or type them into separate --outdir".format(
					clashes="; ".join("{name} = {files}".format(name=name, files=", ".join(files)) for name, files in clashes.items())))
		return unique

	def file_sha256(self, path):
		digest = hashlib.sha256()
		with open(path, "rb") as fin:
			for block in iter(lambda: fin.read(1 << 20), b""):
				digest.update(block)
		return digest.hexdigest()

	def read_result_source(self, resultfile):
		'''(path, sha256) of the query a result file was made from, (None, None) for result
			files from versions that did not record it'''
		with open(resultfile) as handle:
			for row in handle:
				if row.startswith("Query file:"):
					words = row.split(":",1)[1].split()
					sha = next((w.split(":",1)[1] for w in words if w.startswith("sha256:")), None)
					path = " ".join(w for w in words if not w.startswith("sha256:")) or None
					return path, sha
		return None, None

	def query_source_line(self, query):
		'''the trailer line that records which input a result file was made from'''
		return "Query file: {path} sha256:{sha}".format(path=os.path.abspath(query), sha=self.query_sha256)

	def read_result_file(self,resultfile):
		'''Read the final SNP and the SNP path from a sample's result file, used for samples
			skipped because their result file already exists. The final SNP is read from the
			"Final SNP:" line, not from the last SNP of the path: that last SNP is exactly the
			one that can have failed its check, which made the summary disagree with the
			per-sample file (issue #36). An empty or missing file reads as NA/NA (issue #35
			makes empty NA files possible)'''
		final_snp = "NA"
		snp_path = "NA"
		if os.path.isfile(resultfile):
			with open(resultfile) as handle:
				for row in handle:
					if row.startswith("SNP path:"):
						path = row.split(":",1)[1].strip()
						if path and path != "NA":
							snp_path = path
					elif row.startswith("Final SNP:"):
						final_snp = row.split(":",1)[1].strip().split()[0]
		return final_snp,snp_path

	def print_summary(self):
		'''Create summary file and tree. Every sample of the run is listed, NA included, as
			three tab separated columns "ID final_snp snp_path" matching what the run reported
			on stdout: the old summary only contained samples that had a final SNP, silently
			omitting the rest of the run (issue #24)'''
		summarypath = "{outdir}/snps_summary.txt".format(outdir=self.outdir)
		final_snps = set()
		with open(summarypath,"w") as summaryout:
			for q in self.query:
				name = self.sample_name(q)
				if name in self.reused:
					final_snp,snp_path = self.read_result_file("{outdir}/{name}_snps.txt".format(outdir=self.outdir,name=name))
				else:
					final_snp = self.typed.get(name,"NA")
					snp_path = self.result_paths.get(name,"NA")
				print("{name}\t{final_snp}\t{snp_path}".format(name=name,final_snp=final_snp,snp_path=snp_path),file=summaryout)
				if final_snp != "NA":
					final_snps.add(final_snp)

		'''Print summary tree showing all unique SNPs in the final tree'''
		self.create_tree([],"summary",final_snps,True,min_required_hits=self.min_required_hits,summary=True)

	def add_sample_report(self,name,severity,final_snp,snp_path,reasons,states=None,called=None,edge=None,no_data=None):
		'''one row of the results dashboard: what the sample became, why a manual look may be
			wanted, and the SNP states the embedded tree paints when it is selected'''
		self.sample_reports[name] = {"name":name,"severity":severity,"final":final_snp,"path":snp_path,
				"reasons":list(reasons),"states":dict(states or {}),"called":list(called or []),
				"edge":sorted(edge or []),"no_data":sorted(no_data or [])}

	def write_dashboard(self):
		'''summary.html: the results overview for runs with many samples. The samples sit in a
			table with the ones that need a manual look first (failed, not called, calls that
			lean on the strictness tolerance over positions without data, calls whose path runs
			over SNPs flagged at an alignment edge); the database tree is embedded once and
			clicking a sample paints the tree with that sample's states, because the topology
			is the same for every sample'''
		if not self.sample_reports:
			return
		try:
			newickTree = NewickTree(self.database,"summary",self.outdir,min_required_hits=self.min_required_hits,strictness=self.strictness)
		except ValueError as e:
			logger.warning("No results overview written: {error}".format(error=e))
			return
		dashboard = newickTree.write_dashboard([self.sample_reports[name] for name in self.sample_reports],
				"{outdir}/summary.html".format(outdir=self.outdir), self.annotations)
		logger.info("Results overview: open {dashboard} in a browser".format(dashboard=dashboard))

	def run(self,database):
		'''Run CanSNPer2'''
		logger.info("Running CanSNPer2 version-{version}".format(version=__version__))
		if len(self.query) > 0:
			'''Read query input file if a txt file is supplied insead of fasta files'''
			if self.query[0].endswith(".txt"):
				logger.info("Textfile input was found, parsing filepaths in {q} file".format(q=self.query[0]))
				self.query=self.read_query_textfile_input(self.query)
			self.query = self.check_sample_names(self.query)

			'''Main function of CanSNPer2
					1. Align sequences with the selected aligner (progressiveMauve by default)
					2. Read the alignments and find SNPs
					3. Create a tree visualising SNPs found in sequence
					3b. If requested create a list with SNPs and their status
					4. Clean up tmp directory
					'''

			''' Create the reader of the alignment files; it connects to the database in the worker'''
			alignment_reader = self.aligner.reader(
						database=database,
						mask=self.mask,
						verbose=self.verbose)
			'''Walk through the list of queries supplied'''
			if not self.skip_alignment: print("Run alignments for {queries} queries against the references in {refdir} using {program}".format(queries=len(self.query),refdir=self.refdir,program=self.aligner.name))
			for q in self.query:			## For each query file_path
				try:
					self.query_name = self.sample_name(q)  ## get name of file and remove ending

					qfile = q.rsplit("/")[-1]   ## Remove path from query name
					if not os.path.exists(q):
						raise FileNotFoundError("Input file: {qfile} was not found!".format(qfile=q))

					outputfile = "{outdir}/{xmfa}_snps.txt".format(outdir=self.outdir,xmfa=self.query_name)
					self.query_sha256 = self.file_sha256(q)
					if os.path.exists(outputfile) and not self.rerun:
						'''a result is only used again when it was made from the same input: the
							name alone used to let a different genome with the same file name be
							reported with the old call'''
						old_path, old_sha = self.read_result_source(outputfile)
						if old_sha is not None and old_sha != self.query_sha256:
							raise CanSNPer2Error("{outputfile} exists for a different input ({old}); use --rerun to replace it or choose another --outdir".format(
									outputfile=outputfile, old=old_path or "unknown file"))
						if old_sha is None:
							logger.info("{outputfile} does not record its input (older version), it is used again as it is".format(outputfile=outputfile))
						logger.debug("{outputfile} already exits, skip!".format(outputfile=outputfile))
						self.skipped[self.query_name] = "{outputfile} already exists, use --rerun".format(outputfile=outputfile)
						self.reused.add(self.query_name)
						'''the dashboard lists reused samples too, but without states: the tree cannot
							paint a result it did not compute in this run'''
						old_final,old_path = self.read_result_file(outputfile)
						self.add_sample_report(self.query_name,"reused",old_final,old_path,
								["result from an earlier run, use --rerun to recompute it and see it in the tree"])
						continue
					logger.info("Running CanSNPer2 on {query}".format(query=qfile))
					if not self.skip_alignment: ### alignments made before (--skip_alignment) are only read
						logger.info("Run {program} alignments".format(program=self.aligner.name))

					'''For each query fasta align to all CanSNP references the reference folder
						with --skip_alignment the align function only lists the alignment file paths
					'''
					alignment_files = self.align(q)
					logger.debug(alignment_files)
					if len(alignment_files) == 0: ## if keep going is set and the aligner exits with an error continue to next sequence
						logger.debug("{program} exited with a non zero exit status, continue with next sample!".format(program=self.aligner.name))
						logger.warning("{program} error skip {sample}".format(program=self.aligner.name,sample=q))
						self.skipped[self.query_name] = "{program} produced no alignments".format(program=self.aligner.name)
						self.add_sample_report(self.query_name,"failed","NA","NA",["{program} produced no alignments for this sample".format(program=self.aligner.name)])
						self.alignment_files = []
						continue
					'''Read the alignments and find SNPs; returns SNPS (for the visual tree) and SNP_info (text file output)'''
					logger.info("Find SNPs")
					'''the query file maps the query coordinates of the alignments to its sequences when
						the staged copy named in the alignment header is gone (issue #23)'''
					alignment_reader.query_file = q
					try:
						SNPS,SNP_info,called_snps = self.find_snps_multiproc(xmfa_obj=alignment_reader,xmfa_files=alignment_files,export=True,references=list(self.alignment_references))
					except FileNotFoundError as e:
						'''the message has to carry the reason: the default log level only shows
							warnings, and the folder that was searched is what the user needs (issue #176, #177)'''
						logger.warning("One or several xmfa files were not found for {qfile}, continue with next file: {error}".format(qfile=qfile,error=e))
						self.skipped[self.query_name] = str(e)
						self.add_sample_report(self.query_name,"failed","NA","NA",["alignment files missing: {error}".format(error=e)])
						self.alignment_files = []
						continue
					'''If file export is requested print the result for each SNP location to file'''
					if self.export:
						outputfile = "{outdir}/{xmfa}_not_called.txt".format(outdir=self.outdir,xmfa=self.query_name)
						outputfile2 = "{outdir}/{xmfa}_snps.txt".format(outdir=self.outdir,xmfa=self.query_name)

						logger.info("Printing SNP info of non called SNPs to {file}".format(file=outputfile))
						self.csnpdict = {}
						'''Print SNPs to tab separated file'''
						with open(outputfile,"w") as snplist_out:
							print("\t".join(SNP_FILE_COLUMNS),file=snplist_out)
							for snp in SNP_info:
								if snp[0] in called_snps:
									self.csnpdict[snp[0]] = snp
								else:
									print("\t".join(snp),file=snplist_out)

					'''If save tree is requested print tree using ETE3 prints a pdf tree output'''
					SNP = "NA" ## Default message if SNP cannot be confirmed
					tree_result = self.create_tree(SNPS,self.query_name,called_snps,self.save_tree,min_required_hits=self.min_required_hits,strictness=self.strictness)
					if not isinstance(tree_result, tuple):  ## create_tree returns None when it only draws a summary tree, no SNP call to consume
						tree_result = (False,"No SNP call was produced.",[],[])
					final_snp,message,called,no_data_snps = tree_result
					self.result_paths[self.query_name] = ";".join([snp[1] for snp in called]) if called else "NA"
					if final_snp:
						SNP = final_snp[1]
						if not final_snp[2][0]:  ## if snp was never confirmed print NA
							SNP = "NA"
						if self.export:
							with open(outputfile2, "w") as called_out:
								if called:
									print("\t".join(SNP_FILE_COLUMNS),file=called_out)
									for snp in called:
										print("\t".join(self.csnpdict[snp[1]]),file=called_out)
								print("SNP path: {path}".format(path=";".join([snp[1] for snp in called])),file=called_out)
								print("Final SNP: {snp} found/depth: {found}/{depth}".format(snp=SNP,depth=int(final_snp[0]),found=final_snp[2][1]),file=called_out)
								print(self.query_source_line(q),file=called_out)
						logger.info("Final SNP: {snp} found/depth: {found}/{depth}".format(snp=SNP,depth=int(final_snp[0]),found=final_snp[2][1]))
					else:
						if self.export:
							'''an NA sample used to write a snps file with nothing but the final
								line: the SNPs that WERE called on the way to the failing point are
								real information, and with them in this file the two result files
								together list every evaluated position of the tree (issue #35)'''
							with open(outputfile2, "w") as called_out:
								called_ids = set(called_snps)
								called_rows = [row for row in SNP_info if row[0] in called_ids]
								if called_rows:
									print("\t".join(SNP_FILE_COLUMNS),file=called_out)
									for row in called_rows:
										print("\t".join(row),file=called_out)
								print("SNP path: {path}".format(path=";".join([snp[1] for snp in called]) if called else "NA"),file=called_out)
								print("Final SNP: {snp}".format(snp=SNP), file=called_out)
								print(self.query_source_line(q),file=called_out)
							logger.info(message)
					'''SNP positions of the tree that no alignment carried a base for had no data at all,
						and were missing from both result files, so the only place they showed up was the
						log; they are listed in the not called file with an N/A target base (issue #27)'''
					if self.export and no_data_snps:
						logger.info("Listing {n} SNP positions without any alignment data in {file}".format(n=len(no_data_snps),file=outputfile))
						with open(outputfile,"a") as nodata_out:
							for snp in no_data_snps:
								print("\t".join(snp),file=nodata_out)
					'''remember the outcome per sample: the exit code is decided from this
						(issue #191)'''
					if SNP != "NA":
						self.typed[self.query_name] = SNP
					else:
						self.not_typed[self.query_name] = message if not final_snp else "{snp} was not confirmed with enough hits".format(snp=final_snp[1])
					'''the dashboard row for this sample: severity and the reasons a call needs a
						manual look (positions on the path with no alignment data, SNPs flagged at an
						alignment edge, a competing unconfirmed path); states/called let the embedded
						tree paint the sample when it is selected'''
					called_path = [snp[1] for snp in called] if called else []
					edge_flags = set(row[0] for row in SNP_info if len(row) > EDGE_COLUMN and row[EDGE_COLUMN] == "edge")
					edge_on_path = [s for s in edge_flags if s in set(called_path)]
					path_missing = final_snp[4] if final_snp and len(final_snp) > 4 else []
					severity = "clean"
					reasons = []
					if SNP == "NA":
						severity = "not called"
						reasons.append(message if not final_snp else "{snp} was not confirmed with enough hits".format(snp=final_snp[1]))
					if path_missing:
						if severity == "clean":
							severity = "review"
						reasons.append("{k} SNP positions on the called path had no alignment data (the call leans on --strictness): {ids}".format(k=len(path_missing),ids=", ".join(path_missing)))
					if edge_on_path:
						if severity == "clean":
							severity = "review"
						reasons.append("SNPs flagged at an alignment edge on the called path: {ids}".format(ids=", ".join(sorted(edge_on_path))))
					self.add_sample_report(self.query_name,severity,SNP,";".join(called_path) if called_path else "NA",reasons,
							states=SNPS,called=called_snps,edge=sorted(s for s in edge_flags if s in set(called_snps)),
							no_data=[row[0] for row in no_data_snps])
					if self.export:
						print("{query}: {SNP}".format(query=self.query_name, SNP=SNP))
					'''Clean references to aligned xmfa files between queries if several was supplied'''
					self.alignment_files = []
				except Exception as e:
					'''log the real failure before wrapping it, otherwise the cause is invisible (issue #121)'''
					logger.debug("Traceback for {file}".format(file=self.query_name),exc_info=True)
					logger.error("Error while processing {file}: {err}".format(file=self.query_name,err=e))
					self.not_typed[self.query_name] = "failed: {err}".format(err=e)
					if not self.keep_going:
						raise CanSNPer2Error("A file did not run correctly exit CanSNPer2 (use --keep_going to continue with next file!)") from e
					logger.warning("An error occured during processing of {file}, continuing!".format(file=self.query_name))
					self.add_sample_report(self.query_name,"failed","NA","NA",["processing failed: {err}".format(err=e)])

		'''the run has to say what it actually produced: an all NA run used to end with
			"finished successfully" and exit code 0 (issue #191)'''
		if len(self.query) > 0:
			logger.info("{typed} of {n} samples were typed, {na} were not called, {skipped} were skipped".format(
					typed=len(self.typed), n=len(self.query), na=len(self.not_typed), skipped=len(self.skipped)))
			if not self.typed and not self.reused:
				reasons = ["{sample}: {why}".format(sample=sample,why=why) for sample,why in list(self.not_typed.items())+list(self.skipped.items())]
				logger.error("No sample was typed - {reasons}".format(reasons="; ".join(reasons)))
		else:
			logger.warning("No input files were supplied, nothing was typed!")

		if self.summary:
			self.print_summary()
		'''the results dashboard is written for every run: with many samples the list of final
			SNPs on stdout is no longer an overview, and the samples that need a manual look
			have to be found without reading every result file'''
		if len(self.query) > 0:
			self.write_dashboard()
		'''Finally clean up temporary folder when all alignments and trees has been printed!'''
		if not self.keep_temp and len(self.query) > 0: ## if keep temp is turned on do not remove away alignments also if no input files were given
			self.cleanup()
		elif self.own_tmpdir:
			logger.info("Alignments were kept in {tmpdir}, reuse them with --skip_alignment --tmpdir {tmpdir}".format(tmpdir=self.tmpdir))

		if self.typed or self.reused or not self.query:
			logger.info("CanSNPer2 finished successfully, files can be found in {outdir}".format(outdir=self.outdir.rstrip("/")+"/")) ## a trailing slash was already there would print dir// (issue #76)
