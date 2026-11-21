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
import sys
import time

import pytest

from helpers import make_db
from CanSNPer2.DownloadGenomes import DownloadGenomes
from CanSNPer2.modules.Aligners import ProgressiveMauve
from CanSNPer2.modules.CanSNPer2 import AlignmentError, CanSNPer2, CanSNPer2Error
from CanSNPer2.modules.ReferenceCheck import check_reference, read_sequences, reference_snps

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


def snp_at(position):
	'''the SNPs of a reference with one SNP at position of chr1 (ancestral A, derived T)'''
	return [("S.0", "chr1", position, "A", "T")]


def length(path):
	return sum(len(seq) for _name, seq in read_sequences(path))


def fetch(obj, snps=None):
	results = Results()
	obj.download("GCA_000007845.1", "ASM784v1", "Ames", directory=obj.directory, results=results, snps=snps)
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
	assert length(os.path.join(obj.directory, "Ames.fna")) == 10000


def test_a_half_genome_from_an_older_run_is_fetched_again_and_linked(tmp_path):
	obj = downloader(tmp_path, {"data": GZ})
	source = os.path.join(obj.directory, "source")
	os.makedirs(source)
	with open(os.path.join(source, "GCA_000007845.1_ASM784v1_genomic.fna"), "w") as half:
		half.write(">chr1\n" + "ACGT" * 1000 + "\n")    ## 4,000 bases, no link
	refid, ok, reason = fetch(obj, snp_at(9001))
	assert ok and obj.fetched == 1
	assert length(os.path.join(obj.directory, "Ames.fna")) == 10000


def test_a_complete_genome_without_its_link_is_linked_without_fetching(tmp_path):
	obj = downloader(tmp_path, {"data": GZ})
	source = os.path.join(obj.directory, "source")
	os.makedirs(source)
	with open(os.path.join(source, "GCA_000007845.1_ASM784v1_genomic.fna"), "w") as done:
		done.write(GENOME)
	assert fetch(obj, snp_at(9001))[1] and obj.fetched == 0
	assert os.path.islink(os.path.join(obj.directory, "Ames.fna"))


def test_a_genome_shorter_than_its_snps_is_reported(tmp_path):
	obj = downloader(tmp_path, {"data": GZ})
	refid, ok, reason = fetch(obj, snp_at(20001))
	assert not ok and "10000 bases" in reason and "20001" in reason
	assert "Ames.fna" not in files(obj), "not linked, so CanSNPer2 does not use it"


def test_a_worker_that_dies_without_reporting_is_a_failure_not_a_crash(tmp_path, monkeypatch):
	obj = downloader(tmp_path, {"data": GZ})
	obj.get_genomes = lambda *args: ({"GCA_000007845.1": "ASM784v1", "ASM784v1": "Ames"}, ["GCA_000007845.1"])
	monkeypatch.setattr(DownloadGenomes, "download", lambda self, *args, **kwargs: os._exit(3))
	assert obj.run() == 1


def test_typing_stops_on_a_reference_too_short_for_its_snps(tmp_path):
	db = str(tmp_path / "t.db")
	make_db(db, [("S.0", "chr1", 501, "A", "T"), ("S.1", "chr1", 20001, "A", "T")])
	refdir = tmp_path / "refs"
	refdir.mkdir()
	(refdir / "REF.fna").write_text(GENOME)
	query = tmp_path / "sample.fna"
	query.write_text(GENOME)
	obj = CanSNPer2([str(query)], workdir=str(tmp_path), refdir="refs", tmpdir="shared", logdir="logs", outdir="out", database=db)
	with pytest.raises(CanSNPer2Error, match="chr1 of reference REF .* has 10000 bases, but the database places SNPs on it up to position 20001"):
		obj.align(str(query))
	from CanSNPer2.modules.DatabaseConnection import CanSNPdbFunctions
	assert reference_snps(CanSNPdbFunctions(db))["REF"][1] == ("S.1", "chr1", 20001, "A", "T")
	assert check_reference(str(refdir / "REF.fna"), "REF", [("S.1", "chr1", 10000, "T", "A")]) is None


def test_whitespace_is_removed_from_the_inputs(tmp_path):
	# a blank first line makes progressiveMauve run forever; spaces inside a line used to
	# become N and shift every position after them
	obj = CanSNPer2([], workdir=str(tmp_path))
	src = tmp_path / "genome.fna"
	src.write_text("\n>seq1 one\nACGT ACGT\n   \nAC\tGT\n\n>seq2 two\nACGT\n")
	cleaned = obj.clean_input_fasta(str(src), "query")
	assert cleaned != str(src)
	assert open(cleaned).read() == ">seq1 one\nACGTACGT\nACGT\n>seq2 two\nACGT\n"
	assert length(str(src)) == length(cleaned) == 16
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
