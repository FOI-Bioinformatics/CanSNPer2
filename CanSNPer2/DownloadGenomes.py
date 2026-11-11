#!/usr/bin/env python3
import sys, os
import textwrap
import gzip
import shutil
import tempfile
from urllib.request import urlopen
from urllib.error import HTTPError, URLError

import argparse,os
import multiprocessing
from multiprocessing import Queue
'''python 3.14 made forkserver the default start method on linux, where the worker target has
	to be picklable and the bound self.download (with the whole DownloadGenomes state) is not,
	so the fork context is requested explicitly, as in modules/CanSNPer2.py (issue #186).
	Queues must come from the same context as the processes.'''
_MP_CTX = multiprocessing.get_context("fork")
from CanSNPer2.modules.DatabaseConnection import CanSNPdbFunctions
from CanSNPer2.modules.LogSetup import setup_logging
import logging
logger = logging.getLogger(__name__)

DOWNLOAD_TIMEOUT = 120 ## seconds per connection, urlretrieve had no timeout at all (issue #218)

class DownloadGenomes(object):
	"""docstring for DownloadGenomes."""
	def __init__(self, database, source="genbank",directory="references",verbose=False):
		super(DownloadGenomes, self).__init__()
		self.database = CanSNPdbFunctions(database,verbose=verbose)
		self.source = source
		self.directory = directory
		if not os.path.exists(self.directory):
			os.makedirs(self.directory)

	def get_genomes(self, database=False,source="genbank"):
		'''Get the list of genomes in the database'''
		## This is a many to many relation, so all genomes has to be put in a set for each taxonomy id
		genomeDict = {}
		keys = []
		if not database:
			raise ValueError("get_genomes needs a database connection, none was provided")
		QUERY = '''SELECT id,genome,strain,refseq_id,genbank_id,assembly_name FROM {table}'''.format(table="snp_references")
		logger.debug(QUERY)
		logger.debug("Selected source: {source}".format(source=source))
		for id,genome_id,strain,refseq_id,genbank_id,assembly_name in database.query(QUERY).fetchall():
			logging.debug([id,genome_id,strain,refseq_id,genbank_id,assembly_name])
			logging.debug("Downloading {genome}".format(genome=genome_id))
			if source == "refseq":
				genomeDict[refseq_id] = assembly_name
				keys.append(refseq_id)
			else:
				genomeDict[genbank_id] = assembly_name
				keys.append(genbank_id)
			genomeDict[assembly_name] = genome_id
		logging.info(genomeDict)
		return genomeDict,keys

	def retrieve(self, link, target, timeout=DOWNLOAD_TIMEOUT):
		'''Fetch link into target; urlopen takes a timeout, urlretrieve never did (issue #218)'''
		with urlopen(link, timeout=timeout) as response, open(target,"wb") as out:
			shutil.copyfileobj(response, out)

	def download(self,genome_id, assembly, refid ,source="genbank", directory="references", results=None):
		'''Download genomes from refseq or genbank on request
			https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/00x/xxx/GCF_00000xxxx.x_ASMxxxv1/GCF_00000xxxx.x_ASMxxxv1_genomic.fna.gz
			(the ftp:// protocol was retired by NCBI in 2021, the same tree is served over https (issue #120))
			Each genome runs in its own process, so it reports (refid, ok, reason) back through
			results: an exception in a child never reaches the log or the exit code otherwise (issue #218)
		'''
		def report(success, reason=""):
			if results is not None:
				results.put((refid, success, reason))
		sourceD = {"refseq":"F", "genbank": "A"}
		try:
			n1,n2,n3 = textwrap.wrap(genome_id.split("_")[-1].split(".")[0],3)  ## Get the three number parts
		except ValueError:
			logger.error("Could not download {refid}: {genome_id} is not a valid accession".format(refid=refid,genome_id=genome_id))
			report(False, "{genome_id} is not a valid accession".format(genome_id=genome_id))
			return
		link = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GC{source}/{n1}/{n2}/{n3}/{genome_id}_{assembly}/{genome_id}_{assembly}_genomic.fna.gz".format(
						source=sourceD[source],
						n1=n1,
						n2=n2,
						n3=n3,
						genome_id=genome_id,
						assembly=assembly
						)
		retrieved_file = "{directory}/source/{genome_id}_{assembly}_genomic.fna.gz".format(genome_id=genome_id,assembly=assembly,directory=directory)
		ungzipped = retrieved_file[:-3] if retrieved_file.endswith(".gz") else retrieved_file ## strip(".gz") also eats letters of the name, cut the suffix instead (issue #73)
		if not os.path.exists(ungzipped) or os.path.exists(retrieved_file):
			logging.debug("Downloading: "+ link +" <-> "+ retrieved_file)
			'''downloaded into a temp file first, so a failed or interrupted fetch leaves no half
				files and no empty source folder behind (issue #218)'''
			tmpfile = tempfile.mkstemp(suffix=".fna.gz")[1]
			unpacking = False
			try:
				self.retrieve(link, tmpfile)
				os.makedirs("{directory}/source".format(directory=directory), exist_ok=True)
				### unzip file (mauve cant handle zipped sources) with the gzip module, no external gunzip needed (issue #74)
				with gzip.open(tmpfile, "rb") as f_in, open(ungzipped, "wb") as f_out:
					unpacking = True
					shutil.copyfileobj(f_in, f_out)
			except (HTTPError,URLError,OSError) as e:
				logger.error("Could not download {refid} from {link}: {error}".format(refid=refid,link=link,error=e))
				if unpacking and os.path.exists(ungzipped):
					os.remove(ungzipped) ## a half unpacked genome would look complete to the next run (issue #218)
				report(False, str(e))
				return
			finally:
				os.remove(tmpfile)
			bname = os.path.basename(ungzipped)
			try:
				### symlink target is relative to the link location, so no chdir needed (issue #8)
				os.symlink("source/{bname}".format(bname=bname), os.path.join(directory, "{refid}.fna".format(refid=refid)))
			except FileExistsError:
				logger.debug("symlink already exists, skip!")
			report(True)
		else:
			logger.debug("Reference file for {f} already exists!".format(f=retrieved_file.split("/")[-1].strip("_genomic.fna.gz")))
			report(True)
		return

	def run(self):
		'''Start download and collect what each worker actually managed to fetch (issue #218)'''
		self.genomes,keys = self.get_genomes(self.database,self.source)
		if not keys:
			logger.warning("The database lists no genome references, there is nothing to download!")
			return 0
		targets = []
		for genome_id in keys:
			assembly = self.genomes[genome_id]
			targets.append((genome_id, assembly, self.genomes[assembly]))
		results = _MP_CTX.Queue()
		threads = []
		logging.info("Downloading {n} references from {source} into {directory}".format(
				n=len(targets), source=self.source, directory=self.directory))
		for genome_id, assembly, refid in targets:
			kwargs = {"source":self.source, "directory":self.directory, "results":results}
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
			except Queue.Empty:
				break
			reports[refid] = (success, reason)
		failed = []
		for genome_id, assembly, refid in targets:
			success, reason = reports.get(refid, (False, "the download process ended without reporting"))
			if not success:
				failed.append("{refid} ({genome_id} / {assembly}): {reason}".format(
						refid=refid, genome_id=genome_id, assembly=assembly, reason=reason))
		logging.info("{n} of {m} references are now in {directory}".format(
				n=len(targets)-len(failed), m=len(targets), directory=self.directory))
		if failed:
			logger.error("{n} of {m} references could not be downloaded:".format(n=len(failed), m=len(targets)))
			for line in failed:
				logger.error("  "+line)
			return 1
		logging.info("Done!")
		return 0

def main():
	parser = argparse.ArgumentParser(description='CanSNPer2-download')

	baseopts = parser.add_argument_group('Required')
	baseopts.add_argument('-db',  '--database', 	metavar='', required=True,									help='CanSNP database')

	downlopts = parser.add_argument_group('Download options')
	downlopts.add_argument('-s', '--source', 			metavar='', default="genbank",	choices=["genbank","refseq"], 	help="Source for download (genbank/refseq)")
	downlopts.add_argument('-o', '--outdir', 			metavar='', default="references",								help="reference genomes folder")

	debugopts = parser.add_argument_group('Logging and debug options')
	debugopts.add_argument('--logs', metavar='', default='logs', 				help='Specify log directory')
	debugopts.add_argument('--verbose',	action='store_const', const=logging.DEBUG, help='Verbose logging')

	args = parser.parse_args()

	'''fail with a clear message on a bad database path instead of creating an empty database (issue #125)'''
	from CanSNPer2.modules.DatabaseConnection import check_database, ConnectionError
	try:
		check_database(args.database)
	except ConnectionError as e:
		print(e)
		exit(1)

	logval = logging.WARNING
	if args.verbose:
		logval = args.verbose

	args.logs = args.logs.rstrip("/")+"/"
	logpath = setup_logging(args.logs, "CanSNPer2-download", level=logval) ## own tag: it used to append to the CanSNPer2-database log (issue #217)
	logger = logging.getLogger(__name__)
	# create a database connection
	DG = DownloadGenomes(args.database,args.source,args.outdir)
	if DG.run():
		logger.error("CanSNPer2-download finished with missing references, check the failures listed above")
		exit(1)

if __name__ == '__main__':
	main()
