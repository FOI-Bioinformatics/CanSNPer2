'''
The pairwise aligners CanSNPer2 can use to align a query genome against a reference.

An aligner describes how to run the program (command line, input staging, exit codes,
retry) and which AlignmentReader reads its output. CanSNPer2 itself only runs the commands
and hands the alignment files to the reader. To add an aligner, subclass Aligner, give it
a reader (a subclass of AlignmentReader) and add it to ALIGNERS.

What a reader must report is narrow: for every SNP position of the reference that the
alignment covers, the query base on the forward strand of the reference, in the reference
coordinates the database uses (sequences of a multi-sequence reference joined end to end
in file order).
'''

import os
import shutil
import logging

from CanSNPer2.modules.ParseXMFA import ParseXMFA

logger = logging.getLogger(__name__)


class Aligner(object):
	'''Base class of the aligners'''
	name = None            ## name used on the command line (--aligner) and in messages
	program = None         ## executable name, looked up in PATH unless a path is given
	suffix = None          ## file ending of the alignment files
	log_tag = None         ## alignment logs are <reference>_<query>.<log_tag>.log
	stage_inputs = False   ## link the input genomes into the run directory before aligning
	retry_codes = ()       ## return codes of run_alignments that ask align() to retry once

	def __init__(self, path=""):
		'''path is the program itself, or the folder that holds it (a folder may be given
			with a trailing "/", the form the old mauve_path prefix had)'''
		self.path = path or ""

	def executable(self):
		if not self.path:
			return self.program
		if self.path.endswith("/") or os.path.isdir(self.path):
			return os.path.join(self.path, self.program)
		return self.path

	def available(self):
		'''True when the program can be found'''
		return shutil.which(self.executable()) is not None

	def alignment_path(self, tmpdir, reference, query):
		return os.path.join(tmpdir, "{ref}_{query}.{suffix}".format(ref=reference, query=query, suffix=self.suffix))

	def log_path(self, logdir, reference, query):
		return os.path.join(logdir, "{ref}_{query}.{tag}.log".format(ref=reference, query=query, tag=self.log_tag or self.name))

	def command(self, reference_file, query_file, output):
		'''the command as an argument list, so paths with spaces stay single arguments (issue #3)'''
		raise NotImplementedError

	def reader(self, **kwargs):
		'''an AlignmentReader for the alignment files of this aligner'''
		raise NotImplementedError

	def exit_status(self, exitcode):
		'''how a finished alignment process is judged: ("ok", 0), ("retry", code) for a
			failure that a retry with prepare_retry() may fix, or ("error", exitcode)'''
		if exitcode == 0:
			return "ok", 0
		return "error", exitcode

	def retry_message(self, exitcode):
		'''warning logged for a "retry" exit code'''
		return "{name} finished with exitcode {code}".format(name=self.name, code=exitcode)

	def prepare_retry(self, query, tmpdir):
		'''write the query for the retry and return its path, or None when there is no retry'''
		return None


class ProgressiveMauve(Aligner):
	'''progressiveMauve (Darling et al. 2010), writes xmfa'''
	name = "progressiveMauve"
	program = "progressiveMauve"
	suffix = "xmfa"
	log_tag = "mauve"
	'''progressiveMauve writes a multi megabyte "<input>.sslist" progress file next to every input
		genome, which polluted the reference and query folders and made them required to be writable.
		There is no mauve option that turns this off, so the sequences are linked into the run
		directory and the progress files land there and are removed by cleanup() (issue #148)'''
	stage_inputs = True
	retry_codes = (6, 11)

	def command(self, reference_file, query_file, output):
		return [self.executable(), "--output", output, reference_file, query_file]

	def reader(self, **kwargs):
		return ParseXMFA(**kwargs)

	def exit_status(self, exitcode):
		if exitcode == 11:
			return "retry", 11
		if exitcode == -6:
			return "retry", 6
		return super(ProgressiveMauve, self).exit_status(exitcode)

	def retry_message(self, exitcode):
		if exitcode == 6:
			return "Input sequence is not free of gaps, replace gaps with N and retry!!"
		return "WARNING progressiveMauve finished with a exitcode: {code} (seen for genomes with short repetitive contigs or dashes in the sequence)".format(code=exitcode)

	def prepare_retry(self, query, tmpdir):
		'''exitcode 11 is often mauve not handling dashes(-) in the sequence, and exitcode 6 is
			"not free of gaps", so both are retried with every "-" replaced by N (issue #127)'''
		new_query = os.path.join(tmpdir, os.path.basename(query) + ".tmp")
		logger.debug("Replacing all '-' with 'N' in {query}, writing copy to {tmp}".format(query=query,tmp=new_query))
		with open(query, "r") as seq_in, open(new_query, "w") as seq_out:
			for line in seq_in:
				seq_out.write(line.replace("-", "N"))
		return new_query


'''the aligners that can be selected, by name'''
ALIGNERS = {aligner.name: aligner for aligner in (ProgressiveMauve,)}
DEFAULT_ALIGNER = ProgressiveMauve.name


def get_aligner(name=None, path=""):
	'''the aligner called name (case does not matter), set up with the program path'''
	name = name or DEFAULT_ALIGNER
	for known, cls in ALIGNERS.items():
		if known.lower() == name.lower():
			return cls(path)
	raise ValueError("Unknown aligner {name}, choose one of: {known}".format(name=name, known=", ".join(sorted(ALIGNERS))))
