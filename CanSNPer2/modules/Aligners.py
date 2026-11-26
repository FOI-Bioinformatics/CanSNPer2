'''
The pairwise aligners CanSNPer2 can use to align a query genome against a reference.

An aligner describes how to run the program (command line, input staging, exit codes,
retry) and which AlignmentReader reads its output. CanSNPer2 itself only runs the commands
and hands the alignment files to the reader. To add an aligner, subclass Aligner, give it
a reader (a subclass of AlignmentReader) and add it to ALIGNERS.

What a reader must report is narrow: for every SNP position of the reference that the
alignment covers, the query base on the forward strand of the reference, in the reference
coordinates the database uses (the reference sequence and the 1-based position on it).
'''

import os
import shlex
import shutil
import logging

from CanSNPer2.modules.AlignmentReader import SEQUENCES_SUFFIX, write_sequence_table
from CanSNPer2.modules.ParseXMFA import ParseXMFA
from CanSNPer2.modules.ParsePAF import ParsePAF
from CanSNPer2.modules.ParseDelta import ParseDelta
from CanSNPer2.modules.ParsePileup import PileupReader, PLATFORMS, DEFAULT_PLATFORM

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

	def prepare(self, reference_file, query_file, output):
		'''the xmfa counts reference positions along the reference sequences joined end to
			end; the names and lengths of the sequences, in the order mauve joins them, are
			written next to the alignment for ParseXMFA to convert them back'''
		write_sequence_table(reference_file, output + SEQUENCES_SUFFIX)

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
		'''read as bytes: a header that is not UTF-8 stopped the retry with a UnicodeDecodeError'''
		with open(query, "rb") as seq_in, open(new_query, "wb") as seq_out:
			for line in seq_in:
				seq_out.write(line.replace(b"-", b"N"))
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

	def reader(self, **kwargs):
		kwargs.setdefault("min_mapq", self.min_mapq)
		return ParsePAF(**kwargs)


class Mummer4(Aligner):
	'''MUMmer4 nucmer (Marcais et al. 2018) followed by delta-filter -1, writes delta. The
		1-to-1 filter keeps, for every reference and query position, only the best alignment,
		allowing rearrangements; it is the filter MUMmer recommends for SNP finding and plays
		the part of the primary alignment and mapping quality rule of minimap2. nucmer runs with
		its default settings except for one thread (-t 1): with its default of two threads, two
		runs on the same genomes gave different alignments after the 1-to-1 filter (in 2 of 10
		Y. pestis drafts, up to 20 of 970 canSNP positions changed between no data and a base),
		while one thread gave the same output every time. The two programs run as one shell
		command, so a failure of either fails the alignment; delta-filter is taken from the
		folder of nucmer when a path is given.'''
	name = "mummer4"
	program = "nucmer"
	suffix = "delta"
	log_tag = "mummer4"
	options = ["-t", "1"]

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


def is_gzip(path):
	'''True when the file starts with the gzip magic bytes; False for a file that cannot be
		read, which the mapping then reports'''
	try:
		with open(path, "rb") as handle:
			return handle.read(2) == b"\x1f\x8b"
	except OSError:
		return False


def threads_per_mapping(threads, mappings):
	'''the threads of one mapping when the threads of a sample are divided between mappings
		that run at the same time; at least one each, so more mappings than threads use more
		threads than given'''
	return max(1, int(threads) // max(1, int(mappings)))


def write_sites_bed(path, snps):
	'''the canSNP positions of one reference as a BED file for samtools mpileup -l; snps are
		(snp_id, sequence, position, ancestral, derived) with 1-based positions'''
	with open(path, "w") as out:
		for sequence, position in sorted({(row[1], int(row[2])) for row in snps}):
			out.write("{seq}\t{start}\t{end}\n".format(seq=sequence, start=position - 1, end=position))


class ReadMapper(Aligner):
	'''Sequencing reads instead of an assembly: minimap2 maps the reads to the reference with
		the preset of the platform, samtools sorts and indexes them and samtools mpileup lists
		the read bases at the canSNP positions (the BED file <output>.sites.bed, written by
		CanSNPer2 before the command runs). The pileup is the "alignment" PileupReader reads;
		the sorted reads are kept next to it as <output>.bam.

		mpileup is run with -B (no base alignment quality: it lowers base qualities next to
		indels in a way that depends on the reference, and a canSNP is a single base), -A
		(pairs that are not properly paired are counted: with a distant reference many pairs
		are not), -d 0 (no depth limit) and -a (a line for every site, also without reads).
		Reads below min_mapq and bases below min_baseq are left out; secondary, QC failed and
		duplicate-flagged reads are left out by mpileup's defaults, and a base covered by both
		reads of a pair is counted once.

		The three programs run as one bash command with pipefail, so a failure of minimap2
		inside the pipe fails the alignment instead of leaving an empty pileup.

		Only the first max_reads reads of each read file are mapped (read pairs when there are
		two files; 0 maps all of them): head cuts each file in a process substitution, after
		gzip -cd for a compressed one, so no shortened copy is written. The first reads of a
		file from a sequencer come from the whole genome; a file sorted by position (converted
		from a sorted BAM) does not, and needs max_reads 0.

		Not one of the --aligner choices: CanSNPer2 uses it when the query files are reads.'''
	name = "minimap2-reads"
	program = "minimap2"
	suffix = "pileup"
	log_tag = "reads"
	needs_sites = True

	def __init__(self, path="", platform=DEFAULT_PLATFORM, threads=4, samtools="samtools",
			min_mapq=None, min_baseq=None, min_depth=None, min_allele_fraction=None, max_reads=None):
		super(ReadMapper, self).__init__(path)
		if platform not in PLATFORMS:
			raise ValueError("Unknown platform {p}, choose one of: {known}".format(p=platform, known=", ".join(sorted(PLATFORMS))))
		self.platform = platform
		'''the threads of one sample; the run divides them between the references the sample is
			mapped to at the same time (threads_per_mapping)'''
		self.threads = int(threads)
		self.samtools = samtools or "samtools"
		self.given = {"min_mapq": min_mapq, "min_baseq": min_baseq, "min_depth": min_depth,
				"min_allele_fraction": min_allele_fraction, "max_reads": max_reads}

	def setting(self, key):
		'''a threshold: the one given, else the default of the platform'''
		value = self.given.get(key)
		return PLATFORMS[self.platform][key] if value is None else value

	def settings(self):
		'''the settings that change a call, for the settings line of the result files'''
		return "platform={p} min_depth={d} min_allele_fraction={f} min_mapq={q} min_baseq={b} max_reads={n}".format(
				p=self.platform, d=self.setting("min_depth"), f=self.setting("min_allele_fraction"),
				q=self.setting("min_mapq"), b=self.setting("min_baseq"), n=self.setting("max_reads"))

	def available(self):
		return shutil.which(self.executable()) is not None and shutil.which(self.samtools) is not None

	def read_input(self, path):
		'''path as minimap2 is given it: the file itself, or the first max_reads FASTQ records of it
			(4 lines each) through a bash process substitution'''
		n = int(self.setting("max_reads"))
		if n <= 0:
			return shlex.quote(path)
		if is_gzip(path):
			return "<(gzip -cd {p} | head -n {lines})".format(p=shlex.quote(path), lines=4 * n)
		return "<(head -n {lines} {p})".format(p=shlex.quote(path), lines=4 * n)

	def command(self, reference_file, query_file, output, threads=None):
		'''query_file is the read file, or the two files of paired-end reads; threads, when given,
			replaces the threads of the sample for this one mapping'''
		reads = [query_file] if isinstance(query_file, str) else list(query_file)
		bam = output + ".bam"
		q = shlex.quote
		t = str(self.threads if threads is None else int(threads))
		mapping = " ".join([q(a) for a in [self.executable(), "-a", "-x", PLATFORMS[self.platform]["preset"], "-t", t, reference_file]]
				+ [self.read_input(r) for r in reads])
		sort = " ".join(q(a) for a in [self.samtools, "sort", "-@", t, "-o", bam, "-"])
		index = " ".join(q(a) for a in [self.samtools, "index", bam])
		pileup = " ".join(q(a) for a in [self.samtools, "mpileup", "-l", output + ".sites.bed", "-a", "-B", "-A", "-d", "0",
				"-q", str(self.setting("min_mapq")), "-Q", str(self.setting("min_baseq")), "-o", output, bam])
		return ["bash", "-o", "pipefail", "-c", "{m} | {s} && {i} && {p}".format(m=mapping, s=sort, i=index, p=pileup)]

	def reader(self, **kwargs):
		kwargs.setdefault("min_depth", self.setting("min_depth"))
		kwargs.setdefault("min_allele_fraction", self.setting("min_allele_fraction"))
		return PileupReader(**kwargs)


'''the aligners that can be selected, by name'''
ALIGNERS = {aligner.name: aligner for aligner in (ProgressiveMauve, Minimap2, Mummer4)}
DEFAULT_ALIGNER = Minimap2.name


def get_aligner(name=None, path=""):
	'''the aligner called name (case does not matter), set up with the program path'''
	name = name or DEFAULT_ALIGNER
	for known, cls in ALIGNERS.items():
		if known.lower() == name.lower():
			return cls(path)
	raise ValueError("Unknown aligner {name}, choose one of: {known}".format(name=name, known=", ".join(sorted(ALIGNERS))))
