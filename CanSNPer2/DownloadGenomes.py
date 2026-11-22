#!/usr/bin/env python3
import sys, os
import textwrap
import gzip
import zlib
import queue
import shutil
import tempfile
import http.client
from urllib.request import urlopen
from urllib.error import HTTPError, URLError

import argparse,os
import multiprocessing
'''python 3.14 made forkserver the default start method on linux, where the worker target has
	to be picklable and the bound self.download (with the whole DownloadGenomes state) is not,
	so the fork context is requested explicitly, as in modules/CanSNPer2.py (issue #186).
	Queues must come from the same context as the processes.'''
_MP_CTX = multiprocessing.get_context("fork")
from CanSNPer2.modules.TopLevel import run_main
from CanSNPer2.modules.DatabaseConnection import CanSNPdbFunctions
from CanSNPer2.modules.ReferenceCheck import check_reference, reference_snps
from CanSNPer2.modules.LogSetup import setup_logging
from CanSNPer2.modules.FileMode import set_default_mode
from CanSNPer2 import COMMANDS_HELP
import logging
logger = logging.getLogger(__name__)

DOWNLOAD_TIMEOUT = 120 ## seconds per connection, urlretrieve had no timeout at all (issue #218)

def assembly_report_url(genome_url):
	'''the assembly report in the same NCBI folder as the genomic fasta file'''
	suffix = "_genomic.fna.gz"
	base = genome_url[:-len(suffix)] if genome_url.endswith(suffix) else genome_url
	return base + "_assembly_report.txt"

def genbank_names(report_text):
	'''{refseq_accn: genbank_accn} from the text of an NCBI assembly report. The columns are
		found by name in the last comment line; rows without a GenBank or RefSeq accession
		(value na) are skipped. Rows whose Relationship is <> are skipped as well, with one
		warning that names them: the RefSeq and the GenBank sequence differ, so a position on
		one need not hold on the other. Such a sequence keeps its RefSeq name, and the
		reference check then refuses the genome when a SNP needs it, as CanSNPCreator does.'''
	header, names, different = None, {}, []
	for line in report_text.splitlines():
		if line.startswith("#"):
			fields = line.lstrip("#").strip().split("\t")
			if "GenBank-Accn" in fields and "RefSeq-Accn" in fields:
				header = fields
			continue
		if header is None or not line.strip():
			continue
		row = line.split("\t")
		if len(row) < len(header):
			continue
		genbank, refseq = row[header.index("GenBank-Accn")].strip(), row[header.index("RefSeq-Accn")].strip()
		if genbank and refseq and genbank != "na" and refseq != "na":
			if "Relationship" in header and row[header.index("Relationship")].strip() == "<>":
				different.append(refseq)
				continue
			names[refseq] = genbank
	if different:
		logger.warning("{n} RefSeq sequence(s) differ from their GenBank sequence (Relationship <>) and keep "
			"their RefSeq name: {names}; SNPs on them cannot be placed by their GenBank name".format(
			n=len(different), names=", ".join(different)))
	return names

def rename_sequences(path, names):
	'''rewrite the first word of each fasta header found in names, keeping the rest of the line.
		The file is replaced through a temporary file in the same folder, which is removed
		again if anything fails'''
	fd, tmp = tempfile.mkstemp(dir=os.path.dirname(os.path.abspath(path)), suffix=".tmp")
	try:
		set_default_mode(tmp) ## the renamed genome gets the mode of open(), not the 0600 of mkstemp
		## surrogateescape keeps bytes that are not UTF-8 as they are, so a header is not altered or lost
		with os.fdopen(fd, "w", encoding="utf-8", errors="surrogateescape") as out, open(path, encoding="utf-8", errors="surrogateescape") as fin:
			for line in fin:
				if line.startswith(">"):
					words = line[1:].rstrip("\n").split(None, 1)
					if words and words[0] in names:
						rest = (" " + words[1]) if len(words) > 1 else ""
						line = ">" + names[words[0]] + rest + "\n"
				out.write(line)
		os.replace(tmp, path)
	finally:
		if os.path.exists(tmp):
			os.remove(tmp)

def missing_accession(source):
	'''the reason given for a genome without the accession of source'''
	column = "refseq_id" if source == "refseq" else "genbank_id"
	return ("the database has no {column} for this genome; give its assembly accession in the "
		"references file and build the database again").format(column=column)

class DownloadGenomes(object):
	"""docstring for DownloadGenomes."""
	def __init__(self, database, source="genbank",directory="references",verbose=False):
		super(DownloadGenomes, self).__init__()
		'''database None gives a downloader for single genomes through download(), as the
			conversion of SNP source files uses it; run() needs the database'''
		self.database = CanSNPdbFunctions(database,verbose=verbose) if database is not None else None
		self.database_path = database
		self.source = source
		self.directory = directory
		if os.path.exists(self.directory) and not os.path.isdir(self.directory):
			raise NotADirectoryError("--outdir {d} is not a directory, give a folder for the references".format(d=self.directory))
		os.makedirs(self.directory, exist_ok=True)

	def get_genomes(self, database=False,source="genbank"):
		'''The reference genomes of the database, one (accession, assembly name, genome) record per
			row of snp_references; the accession is the refseq_id with source refseq, otherwise the
			genbank_id. Accessions and assembly names used to share one dictionary, so two genomes
			with the same assembly name were downloaded as one, and a genome without accession
			and assembly name stopped the command with a KeyError'''
		records = []
		if not database:
			raise ValueError("get_genomes needs a database connection, none was provided")
		QUERY = '''SELECT id,genome,strain,refseq_id,genbank_id,assembly_name FROM {table}'''.format(table="snp_references")
		logger.debug(QUERY)
		logger.debug("Selected source: {source}".format(source=source))
		for id,genome_id,strain,refseq_id,genbank_id,assembly_name in database.query(QUERY).fetchall():
			logging.debug([id,genome_id,strain,refseq_id,genbank_id,assembly_name])
			accession = refseq_id if source == "refseq" else genbank_id
			records.append((accession, assembly_name, genome_id))
		logging.info(records)
		return records

	def retrieve(self, link, target, timeout=DOWNLOAD_TIMEOUT):
		'''Fetch link into target; urlopen takes a timeout, urlretrieve never did (issue #218)'''
		with urlopen(link, timeout=timeout) as response, open(target,"wb") as out:
			shutil.copyfileobj(response, out)

	def download(self,genome_id, assembly, refid ,source="genbank", directory="references", results=None, snps=None):
		'''Download genomes from refseq or genbank on request
			https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/00x/xxx/GCF_00000xxxx.x_ASMxxxv1/GCF_00000xxxx.x_ASMxxxv1_genomic.fna.gz
			(the ftp:// protocol was retired by NCBI in 2021, the same tree is served over https (issue #120))
			Each genome runs in its own process, so it reports (refid, ok, reason) back through
			results: an exception in a child never reaches the log or the exit code otherwise (issue #218)
			snps: [(snp_id, sequence, position, ancestral, derived), ...] of this reference; a genome
			without a sequence named there, with a sequence shorter than its SNP positions or
			without the ancestral or the derived base at many SNPs is reported as failed
		'''
		def report(success, reason=""):
			if results is not None:
				results.put((refid, success, reason))
		sourceD = {"refseq":"F", "genbank": "A"}
		if not refid or refid in (".", "..") or any(c in str(refid) for c in ("/", "\\", "\0")) or ".." in str(refid):
			'''the name becomes the file name <refid>.fna in the reference folder'''
			reason = "the genome name '{refid}' contains a path separator or '..', rename the genome in the references file and build the database again".format(refid=refid)
			logger.error("Could not download {refid}: {reason}".format(refid=refid, reason=reason))
			report(False, reason)
			return
		if not genome_id:
			'''a NULL accession raised an AttributeError in the worker, which was reported as a
				process that ended without reporting, with a traceback'''
			report(False, missing_accession(source))
			return
		try:
			n1,n2,n3 = textwrap.wrap(genome_id.split("_")[-1].split(".")[0],3)  ## Get the three number parts
		except ValueError:
			logger.error("Could not download {refid}: {genome_id} is not a valid accession".format(refid=refid,genome_id=genome_id))
			report(False, "{genome_id} is not a valid accession".format(genome_id=genome_id))
			return
		link_url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GC{source}/{n1}/{n2}/{n3}/{genome_id}_{assembly}/{genome_id}_{assembly}_genomic.fna.gz".format(
						source=sourceD[source],
						n1=n1,
						n2=n2,
						n3=n3,
						genome_id=genome_id,
						assembly=assembly
						)
		retrieved_file = "{directory}/source/{genome_id}_{assembly}_genomic.fna.gz".format(genome_id=genome_id,assembly=assembly,directory=directory)
		ungzipped = retrieved_file[:-3] if retrieved_file.endswith(".gz") else retrieved_file ## strip(".gz") also eats letters of the name, cut the suffix instead (issue #73)
		link = os.path.join(directory, "{refid}.fna".format(refid=refid))
		if os.path.exists(ungzipped) and not os.path.exists(retrieved_file):
			'''a genome downloaded before is used again, but only when it is long enough for the
				SNPs of the database: older versions could leave a half unpacked file behind and
				took it for complete on the next run'''
			problem = check_reference(ungzipped, refid, snps)
			if problem is None:
				logger.debug("Reference file for {refid} already exists!".format(refid=refid))
				problem = self.link_reference(ungzipped, link, refid)
				report(problem is None, problem or "")
				return
			logger.warning("Downloading {refid} again: {problem}".format(refid=refid, problem=problem))
			os.remove(ungzipped)
		logging.debug("Downloading: "+ link_url +" <-> "+ retrieved_file)
		'''downloaded into a temp file first, so a failed or interrupted fetch leaves no half
			files and no empty source folder behind (issue #218). The genome is unpacked into a
			temp file next to its final name and only renamed once gzip has read the whole stream:
			the gzip trailer carries the CRC and the length of the data, so a transfer that was cut
			short fails here (EOFError) instead of leaving a shorter genome behind'''
		fd, tmpfile = tempfile.mkstemp(suffix=".fna.gz")
		os.close(fd) ## only the name is used, the descriptor was never closed
		partial = ungzipped + ".part"
		try:
			self.retrieve(link_url, tmpfile)
			os.makedirs("{directory}/source".format(directory=directory), exist_ok=True)
			### unzip file (mauve cant handle zipped sources) with the gzip module, no external gunzip needed (issue #74)
			with gzip.open(tmpfile, "rb") as f_in, open(partial, "wb") as f_out:
				shutil.copyfileobj(f_in, f_out)
			os.replace(partial, ungzipped)
			if source == "refseq":
				'''the database names sequences by GenBank accession, the RefSeq assembly by NC_ or NZ_'''
				self.rename_to_genbank(link_url, ungzipped)
		except Exception as e:
			'''any failure of the worker is reported with its reason: IncompleteRead and InvalidURL
				are not OSError, and an exception that is not reported shows only as a process that
				ended without reporting'''
			if isinstance(e, (EOFError, zlib.error, http.client.IncompleteRead)):
				reason = "the download was incomplete ({error})".format(error=e)
			elif isinstance(e, (HTTPError, URLError, OSError, http.client.HTTPException, ValueError)):
				reason = str(e)
			else:
				reason = "{kind}: {error}".format(kind=type(e).__name__, error=e)
			logger.error("Could not download {refid} from {link}: {reason}".format(refid=refid,link=link_url,reason=reason))
			report(False, reason)
			return
		finally:
			os.remove(tmpfile)
			if os.path.exists(partial):
				os.remove(partial) ## a half unpacked genome would look complete to the next run (issue #218)
		problem = check_reference(ungzipped, refid, snps)
		if problem is not None:
			logger.error("Downloaded {refid}, but {problem}".format(refid=refid, problem=problem))
			report(False, problem)
			return
		problem = self.link_reference(ungzipped, link, refid)
		report(problem is None, problem or "")

	def rename_to_genbank(self, link_url, path):
		'''fetch the assembly report of a RefSeq genome and rename its sequences to the GenBank accessions'''
		fd, report_file = tempfile.mkstemp(suffix=".txt")
		os.close(fd)
		try:
			self.retrieve(assembly_report_url(link_url), report_file)
			with open(report_file, encoding="utf-8", errors="replace") as fin: ## a report that is not UTF-8 is read, not fatal
				names = genbank_names(fin.read())
			rename_sequences(path, names)
		except BaseException:
			os.remove(path) ## a genome that keeps the RefSeq names would be reused as renamed by the next run
			raise
		finally:
			os.remove(report_file)

	def link_reference(self, ungzipped, link, refid):
		'''the <refid>.fna link CanSNPer2 reads; a run that found the genome already unpacked
			used to report success without making it'''
		target = "source/{bname}".format(bname=os.path.basename(ungzipped))
		if os.path.lexists(link):
			if not os.path.islink(link):
				return "{link} is a file, not a link to {target}; it was kept, so CanSNPer2 reads that file. Remove it to use the downloaded genome".format(link=link, target=target)
			if os.readlink(link) == target and os.path.exists(link):
				logger.debug("{link} already exists, skip!".format(link=link))
				return
			'''a link to another genome or to nothing is replaced: the new link is made under a
				temporary name and renamed over the old one, so the name never goes missing'''
			logger.warning("{link} pointed to {old}, now {target}".format(link=link, old=os.readlink(link), target=target))
			tmp_link = "{link}.{pid}.tmp".format(link=link, pid=os.getpid())
			if os.path.lexists(tmp_link):
				os.remove(tmp_link)
			os.symlink(target, tmp_link)
			os.replace(tmp_link, link)
			return
		### symlink target is relative to the link location, so no chdir needed (issue #8)
		os.symlink(target, link)

	def run(self):
		'''Start download and collect what each worker actually managed to fetch (issue #218)'''
		database = getattr(self, "database_path", None) or "the database"
		records = self.get_genomes(self.database,self.source)
		if not records:
			'''this was a warning with exit status 0, so a script went on without references'''
			logger.error("{db} lists no reference genomes, there is nothing to download; build the database "
				"again with its references file".format(db=database))
			return 1
		targets = []
		failed = []
		for genome_id, assembly, refid in records:
			if not genome_id:
				failed.append("{refid}: {reason}".format(refid=refid, reason=missing_accession(self.source)))
				continue
			targets.append((genome_id, assembly, refid))
		'''the SNPs of each reference, which its downloaded genome is checked against; without
			them the check of sequence names, lengths and bases was skipped with only a warning'''
		try:
			snps = reference_snps(self.database)
		except Exception as e:
			logger.error("The SNP positions of {db} could not be read ({e}), so the downloaded references "
				"cannot be checked; build the database again with CanSNPer2-database --create".format(db=database, e=e))
			return 1
		results = _MP_CTX.Queue()
		threads = []
		logging.info("Downloading {n} references from {source} into {directory}".format(
				n=len(targets), source=self.source, directory=self.directory))
		for genome_id, assembly, refid in targets:
			kwargs = {"source":self.source, "directory":self.directory, "results":results, "snps":snps.get(refid)}
			p = _MP_CTX.Process(target=self.download, args=(genome_id, assembly, refid), kwargs=kwargs)
			p.start()
			threads.append(p)
		for p in threads:
			p.join() # this blocks until the process terminates
		'''the workers are separate processes: this is the only way their outcome
			reaches the log and the exit code (issue #218)'''
		reports = {}
		for _ in range(len(targets)):
			try:
				refid, success, reason = results.get(timeout=5)
			except queue.Empty: ## Queue from multiprocessing is a function and has no Empty, the worker that died without reporting crashed the whole command
				break
			reports[refid] = (success, reason)
		for genome_id, assembly, refid in targets:
			success, reason = reports.get(refid, (False, "the download process ended without reporting"))
			if not success:
				failed.append("{refid} ({genome_id} / {assembly}): {reason}".format(
						refid=refid, genome_id=genome_id, assembly=assembly, reason=reason))
		total = len(records)
		logging.info("{n} of {m} references are now in {directory}".format(
				n=total-len(failed), m=total, directory=self.directory))
		if failed:
			logger.error("{n} of {m} references could not be downloaded:".format(n=len(failed), m=total))
			for line in failed:
				logger.error("  "+line)
			return 1
		logging.info("Done!")
		return 0

def command():
	parser = argparse.ArgumentParser(
		description='CanSNPer2-download - fetch the reference genomes a database needs from NCBI',
		epilog="Give a database with --database, or name a scheme with --scheme to fetch it from "
			"CanSNPer2-data (or --repo, or --data_dir), build its database and download its references in one step."
			+ COMMANDS_HELP,
		formatter_class=argparse.RawDescriptionHelpFormatter)

	baseopts = parser.add_argument_group('Database')
	baseopts.add_argument('-db',  '--database', 	metavar='',													help='CanSNP database (with --scheme: the database to build, default <scheme>.db)')

	schemeopts = parser.add_argument_group('CanSNPer2-data scheme options')
	schemeopts.add_argument('--scheme', 			metavar='',													help="Fetch this scheme from CanSNPer2-data and build its database")
	schemeopts.add_argument('--branch', 			metavar='',													help="Branch, tag or commit of the repository; in CanSNPer2-data master holds verified schemes, dev also draft schemes (default master, with --repo the default branch of that repository)")
	schemeopts.add_argument('--repo', 				metavar='',													help="GitHub repository laid out as CanSNPer2-data, as owner/name (default FOI-Bioinformatics/CanSNPer2-data)")
	schemeopts.add_argument('--data_dir', 			metavar='',													help="Local directory laid out as CanSNPer2-data (a clone, or a folder of custom schemes) to read schemes from instead of GitHub")
	schemeopts.add_argument('--list_schemes', 		action='store_true',										help="List the schemes of the repository or --data_dir and exit")
	schemeopts.add_argument('--overwrite', 		action='store_true',										help="Replace an existing database when building from --scheme")
	schemeopts.add_argument('--no_references', 	action='store_true',										help="Build the database from --scheme without downloading its references")
	schemeopts.add_argument('--allow_unverified', 	action='store_true',										help="Build a verified scheme from GitHub although its files cannot be checked against a manifest (no manifest, or no hash for a file)")

	downlopts = parser.add_argument_group('Download options')
	downlopts.add_argument('-s', '--source', 			metavar='', default="genbank",	choices=["genbank","refseq"], 	help="Source for download (genbank/refseq)")
	downlopts.add_argument('-o', '--outdir', 			metavar='', default="references",								help="reference genomes folder")

	debugopts = parser.add_argument_group('Logging and debug options')
	debugopts.add_argument('--logs', metavar='', default='logs', 				help='Specify log directory')
	debugopts.add_argument('--verbose',	action='store_const', const=logging.DEBUG, help='Verbose logging')
	debugopts.add_argument('--debug',	action='store_const', const=logging.DEBUG, help='Debug logging, and the traceback of an unexpected error')

	args = parser.parse_args()
	if not (args.database or args.scheme or args.list_schemes):
		parser.error("give a database with --database, a scheme with --scheme, or --list_schemes")
	if args.no_references and not args.scheme:
		'''with only a database there is nothing to build, the references are all the command does'''
		parser.error("--no_references skips the download after building a database from --scheme, give --scheme or leave it out")
	if args.data_dir and (args.repo or args.branch):
		parser.error("--data_dir reads the directory as it is on disk, it cannot be combined with --repo or --branch")

	def scheme_source():
		'''where the schemes are read from: a local directory or a GitHub repository'''
		from CanSNPer2.modules.FetchScheme import GitHubSchemeSource, LocalSchemeSource, DATA_REPO
		if args.data_dir:
			return LocalSchemeSource(args.data_dir)
		if args.repo:
			return GitHubSchemeSource(ref=args.branch, repo=args.repo) ## default branch of that repository unless --branch is given
		return GitHubSchemeSource(ref=args.branch or "master")

	logval = logging.WARNING
	if args.verbose or args.debug:
		logval = logging.DEBUG

	if args.list_schemes:
		'''only prints the schemes, so no log file is written'''
		logging.basicConfig(level=logval, format="%(levelname)s: %(message)s")
		from CanSNPer2.modules.FetchScheme import SchemeError
		try:
			source = scheme_source()
			schemes = source.list_schemes()
		except SchemeError as e:
			print(e)
			exit(1)
		print("SNP schemes in {loc}:".format(loc=source.location()))
		for name, status in schemes:
			print("  {name:<30} {status}".format(name=name, status=status or "status not given"))
		exit()

	args.logs = args.logs.rstrip("/")+"/"
	logpath = setup_logging(args.logs, "CanSNPer2-download", level=logval) ## own tag: it used to append to the CanSNPer2-database log (issue #217)
	logger = logging.getLogger(__name__)

	if args.scheme:
		'''fetch the scheme files, check them against the manifest of CanSNPer2-data and
			build the database with this CanSNPer2 version'''
		from CanSNPer2.modules.FetchScheme import build_scheme, SchemeError
		from CanSNPer2.modules.ModifyTree import InputError
		from CanSNPer2.modules.ReadTaxonomyCanSNPer import InputError as TreeFileError
		from CanSNPer2.modules.DatabaseConnection import SourceFileError, DatabaseWriteError
		from CanSNPer2.modules.DatabaseConnection import ConnectionError as DatabaseConnectionError
		if not args.database:
			args.database = "{scheme}.db".format(scheme=args.scheme)
		try:
			if os.path.exists(args.database) and not args.overwrite:
				raise SchemeError("{db} already exists, use --overwrite to replace it".format(db=args.database))
			source = scheme_source()
			build_scheme(args.scheme, args.database, overwrite=args.overwrite, source=source,
				allow_unverified=args.allow_unverified)
		except (SchemeError, InputError, TreeFileError, SourceFileError, DatabaseWriteError, DatabaseConnectionError) as e:
			logger.error(str(e))
			exit(1)
		print("Built {db} from the scheme {scheme} in {loc}".format(db=args.database, scheme=args.scheme, loc=source.location()))
		if args.no_references:
			exit()

	'''fail with a clear message on a bad database path instead of creating an empty database (issue #125)'''
	from CanSNPer2.modules.DatabaseConnection import check_database, ConnectionError
	try:
		check_database(args.database)
	except ConnectionError as e:
		print(e, file=sys.stderr) ## an error message, so a script that reads stdout does not take it as output
		exit(1)
	if not args.scheme:
		'''which scheme the references are fetched for; a database just built from --scheme
			was already reported by the build (draft status, known issues)'''
		from CanSNPer2.modules.SchemeInfo import log_scheme_notice
		log_scheme_notice(args.database, logger)

	'''the folder is checked before anything is downloaded: a file in its place was found only when the first genome was unpacked'''
	if os.path.exists(args.outdir) and not os.path.isdir(args.outdir):
		logger.error("--outdir {d} exists and is not a directory, give a folder for the references".format(d=args.outdir))
		print("--outdir {d} is not a directory, give a folder for the references".format(d=args.outdir), file=sys.stderr)
		exit(1)
	# create a database connection
	DG = DownloadGenomes(args.database,args.source,args.outdir)
	if DG.run():
		logger.error("CanSNPer2-download did not fetch every reference, see the messages above")
		exit(1)
	print("Downloaded the reference genomes of {db} from {source} into {d}".format(db=args.database, source=args.source, d=args.outdir))


def main():
	'''the entry point of the command: command() with the handler of unexpected errors'''
	return run_main(command)


if __name__ == '__main__':
	main()
