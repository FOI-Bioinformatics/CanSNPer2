import sys, os
import textwrap
from urllib.request import urlretrieve

import argparse,os
from subprocess import Popen
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
			ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/00x/xxx/GCF_00000xxxx.x_ASMxxxv1/GCF_00000xxxx.x_ASMxxxv1_genomic.fna.gz
		'''
		sourceD = {"refseq":"F", "genbank": "A"}
		try:
			n1,n2,n3 = textwrap.wrap(genome_id.split("_")[-1].split(".")[0],3)  ## Get the three number parts
		except ValueError:
			logger.warning("Could not download {refid}".format(refid=refid))
			return
		link = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GC{source}/{n1}/{n2}/{n3}/{genome_id}_{assembly}/{genome_id}_{assembly}_genomic.fna.gz".format(
						source=sourceD[source],
						n1=n1,
						n2=n2,
						n3=n3,
						genome_id=genome_id,
						assembly=assembly
						)
		retrieved_file = "{directory}/source/{genome_id}_{assembly}_genomic.fna.gz".format(genome_id=genome_id,assembly=assembly,directory=directory)
		if not os.path.exists(retrieved_file.strip(".gz")) or os.path.exists(retrieved_file):
			logging.debug("Downloading: "+ link +" <-> "+ retrieved_file)
			urlretrieve(link, retrieved_file)
			### unzip file (mauve cant handle zipped sources)
			p = Popen(["gunzip", "-f", retrieved_file])
			os.chdir(directory)
			bname = os.path.basename(retrieved_file.rstrip(".gz"))
			try:
				os.symlink("source/{bname}".format(bname=bname), "{refid}.fna".format(refid=refid))
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

	logval = logging.WARNING
	if args.verbose:
		logval = args.verbose

	from datetime import date,time
	t = time()
	today = date.today()
	args.logs = args.logs.rstrip("/")+"/"
	logpath = args.logs+"CanSNPer2-database"+today.strftime("%b-%d-%Y")+"{}.log"
	if not os.path.exists(args.logs):
		os.system("mkdir -p {path}".format(path=args.logs))
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
