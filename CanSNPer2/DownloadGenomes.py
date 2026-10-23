#!/usr/bin/env python3
import sys, os
import textwrap
import gzip
import shutil
from urllib.request import urlretrieve

import argparse,os
from multiprocessing import Process
from CanSNPer2.modules.DatabaseConnection import CanSNPdbFunctions
import logging
logger = logging.getLogger(__name__)

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

	def download(self,genome_id, assembly, refid ,source="genbank", directory="references"):
		'''Download genomes from refseq or genbank on request
			https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/00x/xxx/GCF_00000xxxx.x_ASMxxxv1/GCF_00000xxxx.x_ASMxxxv1_genomic.fna.gz
			(the ftp:// protocol was retired by NCBI in 2021, the same tree is served over https (issue #120))
		'''
		sourceD = {"refseq":"F", "genbank": "A"}
		try:
			n1,n2,n3 = textwrap.wrap(genome_id.split("_")[-1].split(".")[0],3)  ## Get the three number parts
		except ValueError:
			logger.warning("Could not download {refid}".format(refid=refid))
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
			urlretrieve(link, retrieved_file)
			### unzip file (mauve cant handle zipped sources) with the gzip module, no external gunzip needed (issue #74)
			with gzip.open(retrieved_file, "rb") as f_in, open(ungzipped, "wb") as f_out:
				shutil.copyfileobj(f_in, f_out)
			os.remove(retrieved_file)
			bname = os.path.basename(ungzipped)
			try:
				### symlink target is relative to the link location, so no chdir needed (issue #8)
				os.symlink("source/{bname}".format(bname=bname), os.path.join(directory, "{refid}.fna".format(refid=refid)))
			except FileExistsError:
				logger.debug("symlink already exists, skip!")
				return
		else:
			logger.debug("Reference file for {f} already exists!".format(f=retrieved_file.split("/")[-1].strip("_genomic.fna.gz")))
		return

	def run(self):
		'''Start download'''
		self.genomes,keys = self.get_genomes(self.database,self.source)
		threads = []
		logging.info("Downloading references")
		if not os.path.isdir("{directory}/source".format(directory=self.directory)): os.mkdir("{directory}/source".format(directory=self.directory))
		for genome_id in keys:
			assembly = self.genomes[genome_id]
			refid = self.genomes[assembly]
			kwargs = {"source":self.source, "directory":self.directory}
			p = Process(target=self.download, args=(genome_id, assembly, refid), kwargs=kwargs)
			p.start()
			threads.append(p)
		for p in threads:
			p.join() # this blocks until the process terminates
		logging.info("Done!")
		return

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

	from datetime import date,time
	t = time()
	today = date.today()
	args.logs = args.logs.rstrip("/")+"/"
	logpath = args.logs+"CanSNPer2-database"+today.strftime("%b-%d-%Y")+"{}.log"
	if not os.path.exists(args.logs):
		os.makedirs(args.logs, exist_ok=True)
	if os.path.exists(logpath):
		logpath=logpath.format("-{:%H:%M}".format(t))
	else: logpath = logpath.format("")
	logging.basicConfig(
			#filename=logpath,
			level=logval,
			format="%(asctime)s %(module)s [%(levelname)-5.5s]  %(message)s",
		    handlers=[
		        logging.FileHandler(logpath),
		        logging.StreamHandler()
		    ])
	logger = logging.getLogger(__name__)
	# create a database connection
	DG = DownloadGenomes(args.database,args.source,args.outdir)
	DG.run()

if __name__ == '__main__':
	main()
