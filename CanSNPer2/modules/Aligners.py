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
import shlex
import shutil
import logging

from CanSNPer2.modules.ParseXMFA import ParseXMFA
from CanSNPer2.modules.ParsePAF import ParsePAF, SEQUENCES_SUFFIX, write_sequence_table
from CanSNPer2.modules.ParseDelta import ParseDelta

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

	def prepare(self, reference_file, query_file, output):
		'''called before an alignment is run (not with --skip_alignment), to write what the
			reader needs next to the alignment besides the alignment itself'''
		return None

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
			'''progressiveMauve aborts (signal 6) on gaps in the sequence, but also when it cannot
				write its .sslist files next to the inputs (an unwritable folder) or runs out of
				memory; the message used to blame gaps only, also for files without any'''
			return "progressiveMauve was aborted (signal 6: gaps in the sequence, an unwritable folder for its .sslist files or too little memory), retrying once with gaps replaced by N"
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


class Minimap2(Aligner):
	'''minimap2 (Li 2018, 2021) with an assembly preset, writes PAF with the long cs tag.
		asm20 allows for up to about 5 % sequence divergence between query and reference.
		Two settings differ from the preset, both to align the short stretches between the
		IS elements and rearrangements that separate Francisella genomes of different
		subspecies:
		- the long-join bandwidth is lowered from 100 kb to 10 kb (-r 1k,10k): with 100 kb
		  about 28 % of a reference stayed unaligned and 153 of 1703 canSNP positions in
		  8 genomes were left without data that progressiveMauve aligned, against 6 with 10 kb;
		- secondary chains are kept down to 10 % of the best score (-N 50 -p 0.1, from 5 and
		  0.8): a short block at a rearrangement breakpoint that overlaps a longer chain on the
		  query was otherwise dropped before it could become a primary alignment. In 208
		  F. tularensis genomes this recovered 60 of the 351 positions only progressiveMauve
		  had aligned, with no position changing its base.
		Only primary alignments with a mapping quality of at least min_mapq are read, see
		ParsePAF, so the extra secondary alignments are written but not used.'''
	name = "minimap2"
	program = "minimap2"
	suffix = "paf"
	log_tag = "minimap2"
	preset = "asm20"
	options = ["-r", "1k,10k", "-N", "50", "-p", "0.1"]
	min_mapq = 20

	def command(self, reference_file, query_file, output):
		return [self.executable(), "-c", "--cs=long", "-x", self.preset] + list(self.options) + ["-o", output, reference_file, query_file]

	def prepare(self, reference_file, query_file, output):
		'''the reference sequence table that maps PAF positions to the joined reference'''
		write_sequence_table(reference_file, output + SEQUENCES_SUFFIX)

	def reader(self, **kwargs):
		kwargs.setdefault("min_mapq", self.min_mapq)
		return ParsePAF(**kwargs)


'''the aligners that can be selected, by name'''
class Mummer4(Aligner):
	'''MUMmer4 nucmer (Marcais et al. 2018) followed by delta-filter -1, writes delta. The
		1-to-1 filter keeps, for every reference and query position, only the best alignment,
		allowing rearrangements; it is the filter MUMmer recommends for SNP finding and plays
		the part of the primary alignment and mapping quality rule of minimap2. nucmer runs with
		its default settings. The two programs run as one shell command, so a failure of either
		fails the alignment; delta-filter is taken from the folder of nucmer when a path is given.'''
	name = "mummer4"
	program = "nucmer"
	suffix = "delta"
	log_tag = "mummer4"
	options = []

	def filter_program(self):
		folder = os.path.dirname(self.executable())
		return os.path.join(folder, "delta-filter") if folder else "delta-filter"

	def available(self):
		return shutil.which(self.executable()) is not None and shutil.which(self.filter_program()) is not None

	def command(self, reference_file, query_file, output):
		'''the delta file names its input files, which the reader opens, so they are given as
			absolute paths: the alignment can then be read from any folder (--skip_alignment)'''
		raw = output + ".unfiltered"
		nucmer = [self.executable()] + list(self.options) + ["--delta=" + raw, os.path.abspath(reference_file), os.path.abspath(query_file)]
		script = "{nucmer} && {filt} -1 {raw} > {out} && rm -f {raw}".format(
				nucmer=" ".join(shlex.quote(a) for a in nucmer), filt=shlex.quote(self.filter_program()),
				raw=shlex.quote(raw), out=shlex.quote(output))
		return ["sh", "-c", script]

	def reader(self, **kwargs):
		return ParseDelta(**kwargs)


ALIGNERS = {aligner.name: aligner for aligner in (ProgressiveMauve, Minimap2, Mummer4)}
DEFAULT_ALIGNER = ProgressiveMauve.name


def get_aligner(name=None, path=""):
	'''the aligner called name (case does not matter), set up with the program path'''
	name = name or DEFAULT_ALIGNER
	for known, cls in ALIGNERS.items():
		if known.lower() == name.lower():
			return cls(path)
	raise ValueError("Unknown aligner {name}, choose one of: {known}".format(name=name, known=", ".join(sorted(ALIGNERS))))
