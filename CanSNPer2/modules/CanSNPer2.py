'''
CanSNPer2 module: A toolkit for SNP-typing using NGS data.
Copyright (C) 2019 David Sundell @ FOI bioinformatics group
'''

import os
import re
import shutil
import hashlib
import datetime
import tempfile
import logging
logger = logging.getLogger(__name__)

## import CanSNPer2 specific modules
from CanSNPer2.modules.Aligners import Aligner, get_aligner
from CanSNPer2.modules.NewickTree import NewickTree, database_provenance, provenance_lines, path_details
from CanSNPer2.modules.DatabaseConnection import CanSNPdbFunctions
from CanSNPer2.modules.ReferenceCheck import check_reference, duplicate_names, reference_snps, listing
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

SEQUENCE_LINE = re.compile(r"^[ACGTUNRYSWKMBDHVacgtunryswkmbdhv*-]+$") ## a line of nucleotide sequence


def ensure_directory(path, option):
	'''Create the folder when it is missing and check that it can be written to. A path that is a file,
		or a folder without write permission, used to be found only after the alignments were made.
		makedirs with exist_ok is safe when another run creates the same folder at the same time.
		Raises CanSNPer2Error naming the option.'''
	if os.path.exists(path) and not os.path.isdir(path):
		raise CanSNPer2Error("{option} {path} is not a directory; give another path".format(option=option, path=path))
	try:
		os.makedirs(path, exist_ok=True)
	except FileExistsError:
		raise CanSNPer2Error("{option} {path} is not a directory; give another path".format(option=option, path=path))
	except OSError as e:
		raise CanSNPer2Error("{option} {path} cannot be created: {error}; give another path".format(option=option, path=path, error=e.strerror or e))
	if not os.access(path, os.W_OK | os.X_OK):
		raise CanSNPer2Error("{option} {path} is not writable; change its permissions or give another path".format(option=option, path=path))


class CanSNPer2(object):
	"""docstring for CanSNPer2"""
	def __init__(self, query, mauve_path="",tmpdir=".tmp", refdir="references",database="CanSNPer.fdb", verbose=False,keep_going=False,**kwargs):
		super(CanSNPer2, self).__init__()
		self.query = query
		self.verbose = verbose
		self.refdir = refdir
		'''a reference genome of the database without its file in refdir stops the run,
			unless this is set: the run then continues with one warning'''
		self.allow_missing_references = kwargs.get("allow_missing_references", False)
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
		'''the support of the call per sample, found/n SNPs on the SNP path, NA when not called'''
		self.sample_supports = {}
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
		ensure_directory(self.workdir, "--workdir")

		'''Resolve the working paths against workdir instead of changing the
			process working directory, which is global state (issue #2)'''
		self.tmpdir = os.path.join(self.workdir, tmpdir)
		if not os.path.exists(self.tmpdir):
			logger.info("Creating tmp directory {tmpdir}".format(tmpdir=self.tmpdir))
		ensure_directory(self.tmpdir, "--tmpdir")
		self.logdir = os.path.join(self.workdir, kwargs.get("logdir","logs"))
		if not os.path.exists(self.logdir):
			logger.info("Creating logs directory {logdir}".format(logdir=self.logdir))
		ensure_directory(self.logdir, "--logdir")

		self.outdir = os.path.join(self.workdir, kwargs.get("outdir","results"))
		if not os.path.exists(self.outdir):
			logger.info("Creating output directory {outdir}".format(outdir=self.outdir))
		ensure_directory(self.outdir, "--outdir")

		'''An aligning run works in its own run_* directory below the shared tmpdir. All runs
			used to write their staged genomes and alignments straight into the shared folder
			(by default /tmp/CanSNPer2), and cleanup() removed every file in it, so two runs on
			the same machine deleted each other's files. With --skip_alignment the given folder is
			where existing alignments are read from, it is used as it is and never removed.'''
		self.own_tmpdir = False
		'''skip_mauve is the old name of skip_alignment'''
		self.skip_alignment = bool(kwargs.get("skip_alignment",False) or kwargs.get("skip_mauve",False))
		self.references = None ## the reference files, listed once by check_start
		if not self.skip_alignment:
			self.tmpdir = tempfile.mkdtemp(prefix="run_", dir=self.tmpdir)
			with open(os.path.join(self.tmpdir, RUN_MARKER), "w"):
				pass
			self.own_tmpdir = True
			logger.info("Temporary files of this run are written to {tmpdir}".format(tmpdir=self.tmpdir))

		self.refdir = os.path.join(self.workdir, self.refdir)

		'''Fetch other key word arguments'''
		'''the per-sample tree page is off by default since the results overview (summary.html)
			carries the tree of every sample; --save_tree writes them again on request'''
		self.save_tree = kwargs.get("save_tree",False)
		'''per-sample files (result tables, tree pages) go into a samples subfolder so that
			summary.html is the one file that greets you in --outdir (the default --outdir is
			itself called results, so the folder under it is not); the folder is created
			lazily, a run whose samples all fail leaves no empty folder behind'''
		self.results_dir = "{outdir}/samples".format(outdir=self.outdir)
		'''tree annotations for the html views: a file path (load_annotations parses it) or
			an already built list of {node, kind, rank, label} dicts. A file that cannot be read
			must not leave the run directory made above behind.'''
		annotations = kwargs.get("annotations",None)
		if isinstance(annotations,str):
			from CanSNPer2.modules.NewickTree import load_annotations
			try:
				annotations = load_annotations(annotations)
			except BaseException:
				self.cleanup()
				raise
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

	def read_database(self):
		'''The checks of the database before any sample is typed. Each of these used to show up
			once per sample, or not at all: a database without SNPs failed with a ValueError for
			every sample, SNPs that could not be read skipped the reference check with a
			warning, and a SNP on a node missing from the tree was never reported. Raises
			CanSNPer2Error; stores the SNPs of each reference genome in self.reference_snps'''
		database = self.database
		self.check_database_file()
		try:
			db = CanSNPdbFunctions(database,verbose=self.verbose)
			try:
				n_references = db.query("SELECT count(*) FROM snp_references").fetchone()[0]
				n_snps = db.query("SELECT count(*) FROM snp_annotation").fetchone()[0]
				n_typed = db.query("SELECT count(*) FROM snp_annotation WHERE typed = 1").fetchone()[0]
				orphans = db.query('''SELECT snp_id, node_id FROM snp_annotation
						WHERE node_id NOT IN (SELECT id FROM nodes) OR node_id NOT IN (SELECT child FROM tree)
						ORDER BY node_id''').fetchall()
			finally:
				db.disconnect()
		except Exception as e:
			raise CanSNPer2Error("The database {db} could not be read: {e}".format(db=database, e=e))
		if n_references == 0:
			raise CanSNPer2Error("The database {db} lists no reference genomes, nothing can be typed against it; build it "
				"again with its references file".format(db=database))
		if n_snps == 0:
			raise CanSNPer2Error("The database {db} has no canSNPs, nothing can be typed against it; build it again with "
				"its canSNP file".format(db=database))
		if n_typed == 0:
			raise CanSNPer2Error("The database {db} has {n} canSNP(s) but none of them is marked typed, nothing can be "
				"called against it; build it again from a canSNP file that marks one canSNP of each node typed".format(
				db=database, n=n_snps))
		if orphans:
			raise CanSNPer2Error("The database {db} has {n} SNP(s) on node ids that are not in the nodes or tree table: "
				"{snps}; export it with cansnper2 db export, correct the files and build it again".format(
				db=database, n=len(orphans),
				snps=listing("{snp} (node {id})".format(snp=snp, id=node_id) for snp, node_id in orphans)))
		return self.read_reference_snps()

	def check_database_file(self):
		'''sqlite creates an empty file for a path that does not exist'''
		if not self.database or not os.path.isfile(self.database):
			raise CanSNPer2Error("The database {database} is not a file".format(database=self.database))

	def read_reference_snps(self):
		'''The SNPs of each reference genome, which the reference files are checked against and
			which give the genomes that need a reference file. SNPs that cannot be read stop the
			run: the reference check used to be skipped with a warning'''
		self.check_database_file()
		try:
			db = CanSNPdbFunctions(self.database,verbose=self.verbose)
			try:
				self.reference_snps = reference_snps(db)
			finally:
				db.disconnect()
		except Exception as e:
			raise CanSNPer2Error("The SNP positions of {db} could not be read ({e}), so the references cannot be checked "
				"against them; export it with cansnper2 db export, correct the files and build it again".format(
				db=self.database, e=e))
		return self.reference_snps

	def check_aligner(self):
		'''Stop before the first sample when the aligner program cannot be found: this used to show
			up after the run had started, once for every sample, as a raw "No such file" error.
			--skip_alignment reads alignments that exist and needs no aligner.'''
		if self.skip_alignment or self.aligner.available():
			return
		raise CanSNPer2Error("The aligner {name} was not found (looked for {program}); install it, give its location with "
			"--aligner_path, or choose another --aligner. CanSNPer2-test lists the aligners that are installed".format(
			name=self.aligner.name, program=self.aligner.executable()))

	def check_start(self):
		'''Check the database, and the reference folder when the query is aligned, before any
			sample is typed. Returns the reference files that are aligned against'''
		self.read_database()
		if self.skip_alignment:
			return []
		self.references = None
		self.references = self.get_references()
		return self.references

	def get_references(self,selected=False):
		'''The reference files in refdir, <genome>.fna, of the reference genomes of the database
			that carry SNPs. Only these are aligned against: a shared reference folder otherwise
			makes mauve align a genome against unrelated species, which wastes hours and fails
			with a non zero exitcode (issue #124). A genome of the database without its file
			stops the run: its SNPs were silently never compared, and when no file matched
			the database every file in the folder was used. With allow_missing_references
			the run continues against the genomes that have a file, with one warning that
			names the others; it still stops when no genome has a file.
		'''
		if getattr(self, "references", None):
			return self.references
		if not hasattr(self, "reference_snps"):
			self.read_reference_snps()
		try:
			refs = sorted(ref for ref in os.listdir(self.refdir) if ref.endswith(".fna"))
		except OSError:
			refs = []
		download = "cansnper2 ref get -d {database} -o {refdir}".format(database=self.database,refdir=self.refdir)
		if not refs:
			'''without references there is nothing to align against, say so instead of failing later
				as "Mauve error skip" for every sample (issue #124)'''
			raise CanSNPer2Error("No reference .fna files were found in {refdir}, run {download} first".format(refdir=self.refdir,download=download))
		genomes = sorted(self.reference_snps)
		missing = [genome for genome in genomes if genome + ".fna" not in refs]
		if missing and (not self.allow_missing_references or len(missing) == len(genomes)):
			raise CanSNPer2Error("{n} reference genome(s) of {database} have no file <genome>.fna in {refdir}: {names}; "
				"the SNPs on them would not be compared. Fetch them with {download}{option}".format(
				n=len(missing),database=self.database,refdir=self.refdir,names=listing(missing),download=download,
				option="" if self.allow_missing_references else ", or give --allow_missing_references to type without them"))
		if missing:
			logger.warning("{n} reference genome(s) of {database} have no file <genome>.fna in {refdir}: {names}; "
				"the SNPs on them are not compared (--allow_missing_references)".format(
				n=len(missing),database=self.database,refdir=self.refdir,names=listing(missing)))
		matched = [genome + ".fna" for genome in genomes if genome not in missing]
		skipped = [ref.rsplit(".",1)[0] for ref in refs if ref not in matched]
		if skipped:
			logger.warning("Skipping {n} reference file(s) in {refdir} that are not reference genomes with SNPs in the database: {names}".format(n=len(skipped),refdir=self.refdir,names=", ".join(skipped)))
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
		with open(path,"rb") as raw:
			start = raw.read(3)
		if start[:2] == b"\x1f\x8b":
			raise CanSNPer2Error("{path} is gzip-compressed, decompress it first".format(path=path))
		if start == b"\xef\xbb\xbf":
			dirty = True ## a byte order mark is removed, aligners read it as part of the first header
		'''utf-8-sig drops a byte order mark, and bytes that are not UTF-8 (a latin-1 character in a header)
			become U+FFFD instead of stopping the run with a decoding error; sequence characters
			that are not ACGT become N below in any case'''
		with open(path,encoding="utf-8-sig",errors="replace") as seq_in:
			for line in seq_in:
				if "\ufffd" in line:
					dirty = True
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
			'''an empty file, a file without a header and a header without sequence gave a silent NA result'''
			if not lines:
				raise CanSNPer2Error("query {path} is empty; give a FASTA file".format(path=path))
			if not any(line.startswith(">") for line in lines):
				raise CanSNPer2Error("query {path} has no FASTA header (a line starting with \">\"); give a FASTA file, "
					"or a list of file paths as the only query argument".format(path=path))
			if not any(not line.startswith(">") for line in lines):
				raise CanSNPer2Error("query {path} has a header but no sequence".format(path=path))
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
				self.read_reference_snps()
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
		query_file = self.not_an_option(self.stage_input(self.clean_input_fasta(query,"query"),"query"))
		for ref in references:	  ## For each reference in the reference folder align to query
			ref_name = ref.rsplit(".",1)[0] ## remove file ending
			output = self.aligner.alignment_path(self.tmpdir.rstrip("/"), ref_name, self.query_name)
			ref_path = "{refdir}/{ref}".format(refdir=self.refdir, ref=ref)
			'''the reference is also checked with --skip_alignment: a renamed reference gave
				NA calls and nothing in the log. The references are listed from refdir and
				opened by clean_input_fasta below, so the file is always present here and the
				isfile test does not leave the check out; --skip_alignment without the
				reference folder fails before this point (a known defect, tracked separately).'''
			if not self.skip_alignment or os.path.isfile(ref_path):
				self.check_reference_file(ref_path, ref_name)
			ref_file = self.not_an_option(self.stage_input(self.clean_input_fasta("{refdir}/{ref}".format(refdir=self.refdir, ref=ref),"ref"),"ref"))
			commands.append(self.aligner.command(ref_file, query_file, output))
			if not self.skip_alignment:
				self.aligner.prepare(ref_file, query_file, output)
			logs.append(self.aligner.log_path(self.logdir, ref_name, self.query_name))
			names.append((ref_name, self.query_name))
			self.alignment_files.append(output)        ## the alignment files are read later
			self.alignment_references.append(ref_name)
		return commands,logs,names

	def not_an_option(self, path):
		'''a relative path that starts with "-" (a query file named -x.fasta) was read by the aligner
			as an option; it is passed as ./-x.fasta'''
		return os.path.join(".", path) if path.startswith("-") else path

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

	def create_tree(self,SNPS,name,called_snps,save_tree,min_required_hits,strictness=0.7, summary=False, snp_info=None, outdir=None):
		'''This function colors the SNP tree of the database with the SNPS found in the query and
			writes the self-contained html viewer of it (NewickTree.draw_tree)
		'''
		outdir = outdir or self.outdir
		if save_tree and not os.path.isdir(outdir):
			os.makedirs(outdir, exist_ok=True)
		newickTree = NewickTree(self.database,name,outdir,min_required_hits=min_required_hits, strictness=strictness)
		try:
			final_snp = newickTree.draw_tree(SNPS,called_snps,save_tree,summary=summary,annotations=self.annotations,
					snp_info=snp_info,run_info=self.run_info())
		finally:
			newickTree.close()
		if save_tree:
			logger.info("{outdir}/{name}_tree.html".format(outdir=outdir, name=name))
		return final_snp

	def is_list_file(self,path):
		'''A query list is an existing text file whose first non-blank line is not a FASTA header
			(">"); a FASTA file is told from it by its content, not by its file name. Empty and
			gzip-compressed files are not lists, the FASTA checks name them.'''
		if not os.path.isfile(path):
			return False
		try:
			with open(path,"rb") as handle:
				head = handle.read(4096)
		except OSError:
			return False
		head = head[3:] if head.startswith(b"\xef\xbb\xbf") else head
		head = head.lstrip()
		return bool(head) and not head.startswith((b">", b"\x1f\x8b"))

	def resolve_queries(self,queries):
		'''The query files of the run. A file whose first non-blank line is not a FASTA header is a
			query list, one file path per line; it must be the only query argument. All other
			arguments are query FASTA files, whatever their file name.'''
		lists = [q for q in queries if self.is_list_file(q)]
		if not lists:
			return list(queries)
		if len(queries) > 1:
			raise CanSNPer2Error("{lst} is a query list (it does not start with a FASTA header) and a list must be the only query argument, "
				"not combined with {n} other query argument(s); put all files in the list, or give the FASTA files without a list".format(
				lst=lists[0], n=len(queries)-1))
		logger.info("A query list was found, reading file paths from {q}".format(q=lists[0]))
		paths = self.read_query_textfile_input(lists[0])
		'''a FASTA file without its header line is read as a list too; its lines are sequence, not
			file paths, and each was reported as a missing query file'''
		if all(SEQUENCE_LINE.match(line) for line in paths[:20]) and not any(os.path.exists(line) for line in paths[:20]):
			raise CanSNPer2Error("{lst} has no FASTA header: its lines are sequence, not file paths of a query list; "
				"add a header line (>name) before the sequence".format(lst=lists[0]))
		return paths

	def read_query_textfile_input(self,query_file):
		'''Read a query list: one file path per line.

			query_file is the query list from the command line, so the text file is its
			first entry; passing the list itself to open() raised TypeError and the
			documented "one fasta per line" input never worked (issue #181). Blank lines
			are skipped, they would otherwise fail the file existence check.'''
		if isinstance(query_file,list):
			query_file = query_file[0]
		with open(query_file, "r") as f:
			query = [query.strip() for query in f if query.strip()]
		if not query:
			raise CanSNPer2Error("The query input file {query_file} is empty, it contained no filepaths".format(query_file=query_file))
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

	def sample_output(self, name, suffix):
		'''path of a per-sample output inside the samples subfolder of --outdir; the folder
			is created here, when the first file actually goes in'''
		path = "{results}/{name}{suffix}".format(results=self.results_dir, name=name, suffix=suffix)
		if not os.path.isdir(self.results_dir):
			os.makedirs(self.results_dir, exist_ok=True)
		return path

	def sample_output_exists(self, name, suffix):
		'''existence check without creating anything (reuse detection)'''
		return os.path.exists("{results}/{name}{suffix}".format(results=self.results_dir, name=name, suffix=suffix))

	def check_sample_names(self, queries):
		'''Every result file, alignment and summary row is named after the sample, so two
			different query files with the same file name (a/x.fasta and b/x.fasta) were typed as
			one sample: the second was skipped as already done and reported with the call of the
			first. Different files with the same name stop the run before anything is aligned;
			the same file given twice is typed once.'''
		unique, seen, clashes = [], {}, {}
		for q in queries:
			name = self.sample_name(q)
			if not name:
				raise CanSNPer2Error("{q} gives an empty sample name (samples are named after the file name without its extension); rename the file".format(q=q))
			if name == "summary" and self.summary and self.save_tree:
				raise CanSNPer2Error("{q} gives the sample name summary, whose tree would overwrite the summary tree of --summary; rename the file".format(q=q))
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

	def query_record(self):
		'''the file next to the alignments of the current sample that records the sha256 of its query'''
		return os.path.join(self.tmpdir, "{name}.query_sha256".format(name=self.query_name))

	def check_alignment_query(self, q):
		'''--skip_alignment finds the alignments by sample name only, so a different file with the same
			name was typed from them. Alignments from versions that did not record their query are
			read as before.'''
		record = self.query_record()
		if not os.path.isfile(record):
			logger.info("{tmpdir} does not record which query {name} was aligned from (older version)".format(tmpdir=self.tmpdir, name=self.query_name))
			return
		with open(record) as handle:
			recorded = handle.read().strip()
		if recorded != self.query_sha256:
			raise CanSNPer2Error("the alignments of {name} in {tmpdir} were made from a different query file; align {q} again without --skip_alignment".format(
				name=self.query_name, tmpdir=self.tmpdir, q=q))

	def file_sha256(self, path):
		digest = hashlib.sha256()
		with open(path, "rb") as fin:
			for block in iter(lambda: fin.read(1 << 20), b""):
				digest.update(block)
		return digest.hexdigest()

	def result_trailer(self, resultfile):
		'''The lines that follow the table of a result file, without their "# " prefix. They are
			written as # comment lines, so a table reader that skips comments reads the table;
			result files of earlier versions wrote them without the #, and are read the same way'''
		with open(resultfile) as handle:
			for row in handle:
				if row.startswith("#"):
					row = row[1:].lstrip(" ")
				yield row

	def read_result_source(self, resultfile):
		'''(path, sha256) of the query a result file was made from, (None, None) for result
			files from versions that did not record it'''
		for row in self.result_trailer(resultfile):
			if row.startswith("Query file:"):
				words = row.split(":",1)[1].split()
				sha = next((w.split(":",1)[1] for w in words if w.startswith("sha256:")), None)
				path = " ".join(w for w in words if not w.startswith("sha256:")) or None
				return path, sha
		return None, None

	def query_source_line(self, query):
		'''the trailer line that records which input a result file was made from'''
		return "Query file: {path} sha256:{sha}".format(path=os.path.abspath(query), sha=self.query_sha256)

	def settings_line(self):
		'''the trailer line that records the database and the settings a result was made with:
			a result file from another database or another threshold was reused without a word'''
		if getattr(self, "_database_sha256", None) is None:
			self._database_sha256 = self.file_sha256(self.database)
		return "Typing settings: database sha256:{sha} min_required_hits={hits} strictness={strict} mask={mask}".format(
			sha=self._database_sha256, hits=self.min_required_hits, strict=self.strictness, mask=self.mask)

	def run_info(self):
		'''What the result pages and snps_summary.txt record about the run: the CanSNPer2
			version, the date, the database file and its sha256, the scheme and reference genomes
			of the database, the aligner and the typing settings'''
		if getattr(self, "_run_info", None) is None:
			info = database_provenance(self.database)
			if getattr(self, "_database_sha256", None) is None:
				self._database_sha256 = self.file_sha256(self.database)
			info.update({"version": __version__,
				"date": getattr(self, "run_date", None) or datetime.datetime.now().strftime("%Y-%m-%d %H:%M"),
				"database": os.path.basename(self.database), "database_sha256": self._database_sha256,
				"aligner": self.aligner.name, "min_required_hits": self.min_required_hits,
				"strictness": self.strictness, "mask": self.mask, "save_tree": bool(self.save_tree),
				"snp_files": bool(self.export)})
			self._run_info = info
		return self._run_info

	def read_result_settings(self, resultfile):
		'''the settings line of a result file, None for files from versions that did not record it'''
		for row in self.result_trailer(resultfile):
			if row.startswith("Typing settings:"):
				return row.rstrip("\n")
		return None

	def read_result_support(self, resultfile):
		'''the support of the call of a result file, "found/depth" as its "Final SNP:" line
			gives it, NA when the line has none (a sample that was not called)'''
		if os.path.isfile(resultfile):
			for row in self.result_trailer(resultfile):
				if row.startswith("Final SNP:") and "found/depth:" in row:
					if row.split(":",1)[1].split()[0] == "NA":
						return "NA"
					return row.split("found/depth:",1)[1].split()[0]
		return "NA"

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
			for row in self.result_trailer(resultfile):
				if row.startswith("SNP path:"):
					path = row.split(":",1)[1].strip()
					if path and path != "NA":
						snp_path = path
				elif row.startswith("Final SNP:"):
					final_snp = row.split(":",1)[1].strip().split()[0]
		return final_snp,snp_path

	def print_summary(self):
		'''Create summary file and tree. Every sample of the run is listed, NA included, matching
			what the run reported on stdout: the old summary only contained samples that had a
			final SNP, silently omitting the rest of the run (issue #24). The file starts with #
			lines on the run (version, date, database, scheme, settings), followed by the header
			row "sample final_snp support snp_path" and one tab separated row per sample; support
			is found/n, the derived SNPs on the SNP path out of the SNPs on it, NA when the
			sample was not called'''
		summarypath = "{outdir}/snps_summary.txt".format(outdir=self.outdir)
		final_snps = set()
		with open(summarypath,"w") as summaryout:
			for line in provenance_lines(self.run_info()):
				print("# " + line,file=summaryout)
			print("\t".join(["sample","final_snp","support","snp_path"]),file=summaryout)
			for q in self.query:
				name = self.sample_name(q)
				if name in self.reused:
					resultfile = self.sample_output(name, "_snps.txt")
					final_snp,snp_path = self.read_result_file(resultfile)
					support = self.read_result_support(resultfile)
				else:
					final_snp = self.typed.get(name,"NA")
					snp_path = self.result_paths.get(name,"NA")
					support = self.sample_supports.get(name,"NA") if final_snp != "NA" else "NA"
				print("{name}\t{final_snp}\t{support}\t{snp_path}".format(name=name,final_snp=final_snp,support=support,snp_path=snp_path),file=summaryout)
				if final_snp != "NA":
					final_snps.add(final_snp)

		'''Print summary tree showing all unique SNPs in the final tree'''
		self.create_tree([],"summary",final_snps,True,min_required_hits=self.min_required_hits,summary=True)

	def add_sample_report(self,name,severity,final_snp,snp_path,reasons,states=None,called=None,edge=None,no_data=None,
			support="NA",path_no_data=None,path_snps=None):
		'''one row of the results dashboard: what the sample became, why a manual look may be
			wanted, and the SNP states the embedded tree paints when it is selected. support is
			found/n, the derived SNPs on the SNP path out of the SNPs on it (NA when not called),
			path_no_data the SNPs on the SNP path without alignment data and path_snps the SNP
			path with the bases of each SNP (NewickTree.path_details)'''
		self.sample_reports[name] = {"name":name,"severity":severity,"final":final_snp,"path":snp_path,
				"reasons":list(reasons),"states":dict(states or {}),"called":list(called or []),
				"edge":sorted(edge or []),"no_data":sorted(no_data or []),"support":support,
				"path_no_data":list(path_no_data or [])}
		if path_snps is not None:
			self.sample_reports[name]["path_snps"] = list(path_snps)

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
		try:
			dashboard = newickTree.write_dashboard([self.sample_reports[name] for name in self.sample_reports],
					"{outdir}/summary.html".format(outdir=self.outdir), self.annotations, run_info=self.run_info())
		finally:
			newickTree.close()
		logger.info("Results overview: open {dashboard} in a browser".format(dashboard=dashboard))

	def run(self,database):
		'''Run CanSNPer2. On a failure the run directory is removed unless --keep_temp is set; an
			unexpected error used to leave it behind.'''
		try:
			return self._run(database)
		except BaseException:
			if not self.keep_temp:
				self.cleanup()
			raise

	def _run(self,database):
		logger.info("Running CanSNPer2 version-{version}".format(version=__version__))
		self.run_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
		if len(self.query) > 0:
			'''A query list file (file paths, one per line) is read, see resolve_queries'''
			self.query = self.resolve_queries(self.query)
			self.query = self.check_sample_names(self.query)
			'''a missing query used to be found only after the queries before it were aligned'''
			unusable = [q for q in self.query if not os.path.isfile(q)]
			if unusable and not self.keep_going:
				shown = ", ".join(unusable[:5]) + (" and {n} more".format(n=len(unusable)-5) if len(unusable) > 5 else "")
				raise CanSNPer2Error("{n} query file(s) were not found or are not files: {files}".format(n=len(unusable), files=shown))
			'''the database and the reference folder are checked once, before the first sample'''
			self.check_start()
			'''the aligner is checked only when a sample will be aligned: a run whose samples all have a
				result from an earlier run (made again for a summary, for example) needs no aligner'''
			if self.rerun or any(not self.sample_output_exists(self.sample_name(q), "_snps.txt") for q in self.query):
				self.check_aligner()

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
					if os.path.isdir(q):
						raise CanSNPer2Error("Input file: {qfile} is a directory, not a FASTA file".format(qfile=q))
					if not os.path.exists(q):
						raise FileNotFoundError("Input file: {qfile} was not found!".format(qfile=q))

					outputfile = self.sample_output(self.query_name, "_snps.txt")
					self.query_sha256 = self.file_sha256(q)
					if os.path.exists(outputfile) and not self.rerun:
						'''a result is only used again when it was made from the same input: the
							name alone used to let a different genome with the same file name be
							reported with the old call'''
						old_path, old_sha = self.read_result_source(outputfile)
						if old_sha is not None and old_sha != self.query_sha256:
							raise CanSNPer2Error("{outputfile} exists for a different input ({old}); use --rerun to replace it or choose another --outdir".format(
									outputfile=outputfile, old=old_path or "unknown file"))
						old_settings = self.read_result_settings(outputfile)
						if old_settings is not None and old_settings != self.settings_line():
							raise CanSNPer2Error("{outputfile} was made with a different database or settings ({old}); use --rerun to replace it or choose another --outdir".format(
									outputfile=outputfile, old=old_settings.split(":",1)[1].strip()))
						if old_sha is None:
							logger.info("{outputfile} does not record its input (older version), it is used again as it is".format(outputfile=outputfile))
						logger.debug("{outputfile} already exits, skip!".format(outputfile=outputfile))
						'''the skip used to be a debug line, so a run that typed nothing said nothing
							after announcing the alignments (issue #78); the call of the earlier run is
							printed like a computed one, with what to do to recompute it'''
						self.skipped[self.query_name] = "{outputfile} already exists, use --rerun".format(outputfile=outputfile)
						self.reused.add(self.query_name)
						'''the dashboard lists reused samples too, but without states: the tree cannot
							paint a result it did not compute in this run'''
						old_final,old_path = self.read_result_file(outputfile)
						print("{query}: {SNP} (from an earlier run, use --rerun to type it again)".format(
								query=self.query_name, SNP=old_final))
						self.add_sample_report(self.query_name,"reused",old_final,old_path,
								["result from an earlier run, use --rerun to recompute it and see it in the tree"],
								support=self.read_result_support(outputfile))
						continue
					logger.info("Running CanSNPer2 on {query}".format(query=qfile))
					if not self.skip_alignment: ### alignments made before (--skip_alignment) are only read
						logger.info("Run {program} alignments".format(program=self.aligner.name))

					'''For each query fasta align to all CanSNP references the reference folder
						with --skip_alignment the align function only lists the alignment file paths
					'''
					if self.skip_alignment:
						self.check_alignment_query(q)
					alignment_files = self.align(q)
					logger.debug(alignment_files)
					if alignment_files and not self.skip_alignment:
						with open(self.query_record(), "w") as record:
							print(self.query_sha256, file=record)
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
							warnings, and the folder that was searched is what the user needs (issue #176, #177).
							A sample without its alignments is a failed sample: it stops the run unless
							--keep_going is set, as other failures do, instead of being skipped quietly.'''
						self.alignment_files = []
						raise CanSNPer2Error("alignment files of {qfile} are missing: {error}".format(qfile=qfile,error=e)) from e
					'''If file export is requested the result of each SNP position is written to the two
						result files after the call; the rows of the called SNPs are collected here'''
					if self.export:
						outputfile = self.sample_output(self.query_name, "_not_called.txt")
						outputfile2 = self.sample_output(self.query_name, "_snps.txt")
						self.csnpdict = {snp[0]: snp for snp in SNP_info if snp[0] in called_snps}

					'''Color the tree with the states of the sample and call the final SNP'''
					SNP = "NA" ## Default message if SNP cannot be confirmed
					tree_result = self.create_tree(SNPS,self.query_name,called_snps,self.save_tree,min_required_hits=self.min_required_hits,strictness=self.strictness,snp_info=SNP_info,outdir=self.results_dir)
					if not isinstance(tree_result, tuple):  ## create_tree returns None when it only draws a summary tree, no SNP call to consume
						tree_result = (False,"No SNP call was produced.",[],[])
					final_snp,message,called,no_data_snps = tree_result
					'''The SNP path of a confirmed call is the chain of the final SNP: the called SNPs from the
						root of the tree down to the final SNP. Derived SNPs that were called outside that chain
						(a competing branch) are not part of the path, they are listed in the dashboard
						report for review. The table of called SNPs keeps one row for every called SNP.
						Without a confirmed call (NA) the path lists all called SNPs, as before.'''
					called_names = [snp[1] for snp in called] if called else []
					path_names = called_names
					off_chain = []
					if final_snp and final_snp[2][0] and len(final_snp) > 5:
						chain = final_snp[5]
						path_names = [name for name in chain if name in set(called_names)]
						off_chain = sorted(name for name in called_names if name not in set(chain))
					self.result_paths[self.query_name] = ";".join(path_names) if path_names else "NA"
					'''the support of a call: the derived SNPs on the SNP path out of the SNPs on it'''
					support = "NA"
					if final_snp:
						SNP = final_snp[1]
						if not final_snp[2][0]:  ## if snp was never confirmed print NA
							SNP = "NA"
						else:
							support = "{found}/{n}".format(found=final_snp[2][1], n=final_snp[3] if len(final_snp) > 3 else int(final_snp[0]))
						called_rows = [self.csnpdict[snp[1]] for snp in called] if self.export else []
						final_line = "Final SNP: {snp} found/depth: {found}/{depth}".format(snp=SNP,depth=int(final_snp[0]),found=final_snp[2][1])
						path_text = ";".join(path_names)
						logger.info(final_line)
					else:
						'''an NA sample used to write a snps file with nothing but the final
							line: the SNPs that WERE called on the way to the failing point are
							real information, and with them in this file the two result files
							together list every evaluated position of the tree (issue #35)'''
						called_ids = set(called_snps)
						called_rows = [row for row in SNP_info if row[0] in called_ids]
						final_line = "Final SNP: {snp}".format(snp=SNP)
						path_text = ";".join(path_names) if path_names else "NA"
						logger.info(message)
					self.sample_supports[self.query_name] = support
					if self.export:
						'''Both files start with the table, a header row and the rows; the lines on the
							call and the input follow as # comment lines, so that a table reader that skips
							comments (pandas.read_csv(comment="#")) reads the table. SNP positions of the tree
							that no alignment carried a base for are listed in the not called file with an
							N/A target base (issue #27)'''
						with open(outputfile2, "w") as called_out:
							print("\t".join(SNP_FILE_COLUMNS),file=called_out)
							for row in called_rows:
								print("\t".join(row),file=called_out)
							for line in ("SNP path: {path}".format(path=path_text), final_line, self.query_source_line(q), self.settings_line()):
								print("# " + line,file=called_out)
						logger.info("Printing SNP info of non called SNPs to {file}".format(file=outputfile))
						if no_data_snps:
							logger.info("Listing {n} SNP positions without any alignment data in {file}".format(n=len(no_data_snps),file=outputfile))
						with open(outputfile,"w") as snplist_out:
							print("\t".join(SNP_FILE_COLUMNS),file=snplist_out)
							for snp in SNP_info:
								if snp[0] not in called_snps:
									print("\t".join(snp),file=snplist_out)
							for snp in no_data_snps:
								print("\t".join(snp),file=snplist_out)
							for line in (self.query_source_line(q), self.settings_line()):
								print("# " + line,file=snplist_out)
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
					called_path = path_names
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
						reasons.append("{k} SNP positions on the SNP path had no alignment data (the call leans on --strictness): {ids}".format(k=len(path_missing),ids=", ".join(path_missing)))
					if edge_on_path:
						if severity == "clean":
							severity = "review"
						reasons.append("SNPs flagged at an alignment edge on the SNP path: {ids}".format(ids=", ".join(sorted(edge_on_path))))
					if off_chain and SNP != "NA":
						if severity == "clean":
							severity = "review"
						reasons.append("{k} derived SNP(s) called outside the SNP path, on another branch: {ids}".format(k=len(off_chain),ids=", ".join(off_chain)))
					self.add_sample_report(self.query_name,severity,SNP,";".join(called_path) if called_path else "NA",reasons,
							states=SNPS,called=called_snps,edge=sorted(s for s in edge_flags if s in set(called_snps)),
							no_data=[row[0] for row in no_data_snps],support=support,
							path_no_data=path_missing if SNP != "NA" else [],path_snps=path_details(called_path,SNP_info))
					print("{query}: {SNP}".format(query=self.query_name, SNP=SNP)) ## the call is printed whether or not result files are written
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
		if not self.keep_temp: ## a run without queries also removes its (empty) run directory
			self.cleanup()
		elif self.own_tmpdir:
			logger.info("Alignments were kept in {tmpdir}, reuse them with --skip_alignment --tmpdir {tmpdir}".format(tmpdir=self.tmpdir))

		if self.typed or self.reused or not self.query:
			logger.info("CanSNPer2 finished successfully, files can be found in {outdir}".format(outdir=self.outdir.rstrip("/")+"/")) ## a trailing slash was already there would print dir// (issue #76)
