"""Reference downloads that are cut short, references too short for their SNPs, blank lines
and whitespace in the inputs, and alignments that never return.

- A transfer that ended early left a half unpacked genome: gzip raised EOFError, which was not
  caught, the worker died without reporting, and run() itself crashed on `Queue.Empty`
  (multiprocessing.Queue is a function). The next run took the half file for complete,
  reported success and never made the <reference>.fna link.
- A reference shorter than the last SNP position on it gave missing or wrong genotypes and no
  error.
- progressiveMauve never returns on a query that starts with a blank line, and CanSNPer2 had
  no time limit on the alignments."""
import gzip
import os
import sqlite3
import sys
import time

import pytest

from CanSNPer2.DownloadGenomes import DownloadGenomes
from CanSNPer2.modules.Aligners import ProgressiveMauve
from CanSNPer2.modules.CanSNPer2 import AlignmentError, CanSNPer2, CanSNPer2Error
from CanSNPer2.modules.ReferenceCheck import check_reference, last_snp_positions, sequence_length

GENOME = ">chr1 test\n" + "ACGT" * 2500 + "\n"         ## 10,000 bases
GZ = gzip.compress(GENOME.encode())


def downloader(tmp_path, payload):
	'''a DownloadGenomes whose retrieve() writes payload instead of fetching from NCBI'''
	obj = DownloadGenomes.__new__(DownloadGenomes)
	obj.source, obj.directory, obj.database = "genbank", str(tmp_path / "refs"), None
	os.makedirs(obj.directory, exist_ok=True)
	obj.fetched = 0
	def retrieve(link, target, timeout=None):
		obj.fetched += 1
		with open(target, "wb") as out:
			out.write(payload["data"])
	obj.retrieve = retrieve
	return obj


class Results(list):
	def put(self, item):
		self.append(item)


def fetch(obj, last_position=None):
	results = Results()
	obj.download("GCA_000007845.1", "ASM784v1", "Ames", directory=obj.directory, results=results, last_position=last_position)
	return results[0]


def files(obj):
	return sorted(os.path.relpath(os.path.join(r, f), obj.directory) for r, _, fs in os.walk(obj.directory) for f in fs)


def test_a_transfer_cut_short_leaves_nothing_behind(tmp_path):
	payload = {"data": GZ[: len(GZ) // 2]}
	obj = downloader(tmp_path, payload)
	refid, ok, reason = fetch(obj)
	assert not ok and "incomplete" in reason
	assert files(obj) == [], "no half genome, no .part file and no link"
	# the next run fetches again and succeeds
	payload["data"] = GZ
	assert fetch(obj)[1]
	assert files(obj) == ["Ames.fna", "source/GCA_000007845.1_ASM784v1_genomic.fna"]
	assert sequence_length(os.path.join(obj.directory, "Ames.fna")) == 10000


def test_a_half_genome_from_an_older_run_is_fetched_again_and_linked(tmp_path):
	obj = downloader(tmp_path, {"data": GZ})
	source = os.path.join(obj.directory, "source")
	os.makedirs(source)
	with open(os.path.join(source, "GCA_000007845.1_ASM784v1_genomic.fna"), "w") as half:
		half.write(">chr1\n" + "ACGT" * 1000 + "\n")    ## 4,000 bases, no link
	refid, ok, reason = fetch(obj, last_position=9000)
	assert ok and obj.fetched == 1
	assert sequence_length(os.path.join(obj.directory, "Ames.fna")) == 10000


def test_a_complete_genome_without_its_link_is_linked_without_fetching(tmp_path):
	obj = downloader(tmp_path, {"data": GZ})
	source = os.path.join(obj.directory, "source")
	os.makedirs(source)
	with open(os.path.join(source, "GCA_000007845.1_ASM784v1_genomic.fna"), "w") as done:
		done.write(GENOME)
	assert fetch(obj, last_position=9000)[1] and obj.fetched == 0
	assert os.path.islink(os.path.join(obj.directory, "Ames.fna"))


def test_a_genome_shorter_than_its_snps_is_reported(tmp_path):
	obj = downloader(tmp_path, {"data": GZ})
	refid, ok, reason = fetch(obj, last_position=20000)
	assert not ok and "10000 bases" in reason and "20000" in reason
	assert "Ames.fna" not in files(obj), "not linked, so CanSNPer2 does not use it"


def test_a_worker_that_dies_without_reporting_is_a_failure_not_a_crash(tmp_path, monkeypatch):
	obj = downloader(tmp_path, {"data": GZ})
	obj.get_genomes = lambda *args: ({"GCA_000007845.1": "ASM784v1", "ASM784v1": "Ames"}, ["GCA_000007845.1"])
	monkeypatch.setattr(DownloadGenomes, "download", lambda self, *args, **kwargs: os._exit(3))
	assert obj.run() == 1


def make_db(path, positions, genome="REF"):
	conn = sqlite3.connect(path)
	conn.executescript("""
		CREATE TABLE nodes (id integer PRIMARY KEY, name text NOT NULL);
		CREATE TABLE tree (parent integer, child integer, rank_i integer);
		CREATE TABLE snp_annotation (node_id integer PRIMARY KEY, snp_id varchar(6), position integer,
			ancestral_base varchar(1), derived_base varchar(1), reference varchar(20),
			date datetime, genome_i integer);
		CREATE TABLE snp_references (id integer PRIMARY KEY, genome varchar(6), strain varchar(6));
	""")
	conn.execute("INSERT INTO snp_references (id, genome, strain) VALUES (5, ?, ?)", (genome, genome))
	for i, pos in enumerate(positions):
		conn.execute("INSERT INTO nodes (id, name) VALUES (?, ?)", (i + 1, "S.{0}".format(i)))
		conn.execute("INSERT INTO snp_annotation (node_id, snp_id, position, ancestral_base, derived_base, reference, date, genome_i)"
			" VALUES (?, ?, ?, 'A', 'T', 'test', '06/11/26', 5)", (i + 1, "S.{0}".format(i), pos))
	conn.commit()
	conn.close()


def test_typing_stops_on_a_reference_too_short_for_its_snps(tmp_path):
	db = str(tmp_path / "t.db")
	make_db(db, [500, 20000])
	refdir = tmp_path / "refs"
	refdir.mkdir()
	(refdir / "REF.fna").write_text(GENOME)
	query = tmp_path / "sample.fna"
	query.write_text(GENOME)
	obj = CanSNPer2([str(query)], workdir=str(tmp_path), refdir="refs", tmpdir="shared", logdir="logs", outdir="out", database=db)
	with pytest.raises(CanSNPer2Error, match="10000 bases, but the database places SNPs up to position 20000"):
		obj.align(str(query))
	from CanSNPer2.modules.DatabaseConnection import CanSNPdbFunctions
	assert last_snp_positions(CanSNPdbFunctions(db)) == {"REF": 20000}
	assert check_reference(str(refdir / "REF.fna"), "REF", 10000) is None


def test_whitespace_is_removed_from_the_inputs(tmp_path):
	# a blank first line makes progressiveMauve run forever; spaces inside a line used to
	# become N and shift every position after them
	obj = CanSNPer2([], workdir=str(tmp_path))
	src = tmp_path / "genome.fna"
	src.write_text("\n>seq one\nACGT ACGT\n   \nAC\tGT\n\n>seq two\nACGT\n")
	cleaned = obj.clean_input_fasta(str(src), "query")
	assert cleaned != str(src)
	assert open(cleaned).read() == ">seq one\nACGTACGT\nACGT\n>seq two\nACGT\n"
	assert sequence_length(str(src)) == sequence_length(cleaned) == 16
	# a file whose only fault is a blank first line must be cleaned as well
	src2 = tmp_path / "lead.fna"
	src2.write_text("\n>seq\nACGT\n")
	cleaned2 = obj.clean_input_fasta(str(src2), "query")
	assert cleaned2 != str(src2) and open(cleaned2).read() == ">seq\nACGT\n"


def test_alignments_that_never_return_are_stopped(tmp_path):
	logs = tmp_path / "logs"
	logs.mkdir()
	slow = [sys.executable, "-c", "import time; time.sleep(60)", "ref_REF.fna", "query_q.fna"]
	fast = [sys.executable, "-c", "pass", "ref_OTHER.fna", "query_q.fna"]
	obj = CanSNPer2([], workdir=str(tmp_path), alignment_timeout=0.02)   ## 1.2 seconds
	start = time.time()
	with pytest.raises(AlignmentError, match="1 progressiveMauve processes did not run correctly"):
		obj.run_alignments([slow, fast], [str(logs / "a.log"), str(logs / "b.log")])
	assert time.time() - start < 20
	obj = CanSNPer2([], workdir=str(tmp_path), alignment_timeout=0.02, keep_going=True)
	assert obj.run_alignments([slow], [str(logs / "c.log")]) == 1


def test_the_abort_message_does_not_blame_gaps_only():
	message = ProgressiveMauve().retry_message(6)
	assert "unwritable folder" in message and "not free of gaps" not in message
