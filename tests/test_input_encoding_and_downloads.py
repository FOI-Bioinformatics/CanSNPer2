"""Encoding of reference and SNP files, and the failure paths of the reference download.

- A UTF-8 byte order mark on the first header dropped the first record of a reference, and a
  non-UTF-8 byte or a gzipped file stopped with a bare UnicodeDecodeError that names no file.
- A <refid>.fna link to another file, or to nothing, was kept and the download reported as done.
- IncompleteRead and ValueError were not caught in a download worker, and a non-UTF-8 assembly
  report stopped it without a report; the descriptors of mkstemp were not closed.
- --outdir naming a file was found after the genome had been downloaded, and a genome name
  with a path separator was used as a file name."""
import gzip
import http.client
import os
import subprocess
import sys

import pytest

from CanSNPer2.DownloadGenomes import DownloadGenomes
from CanSNPer2.modules.Aligners import ProgressiveMauve
from CanSNPer2.modules.ConvertPositions import (ConversionError, JoinedReference, convert_table,
	references_from_file, table_genomes)
from CanSNPer2.modules.ReferenceCheck import check_reference, read_sequences, sequence_names

FASTA = b">chr1 first\nACGTACGT\n>chr2 second\nGGGG\n"
GZ = gzip.compress(FASTA)


def write(path, data):
	with open(str(path), "wb") as out:
		out.write(data)
	return str(path)


def test_bom_on_first_header_keeps_the_first_record(tmp_path):
	path = write(tmp_path / "bom.fna", b"\xef\xbb\xbf" + FASTA)
	assert [name for name, _seq in read_sequences(path)] == ["chr1", "chr2"]
	assert sequence_names(path) == ["chr1", "chr2"]


def test_non_utf8_header_is_read(tmp_path):
	path = write(tmp_path / "latin.fna", b">chr1 caf\xe9\nACGT\n")
	assert read_sequences(path) == [("chr1", "ACGT")]
	assert sequence_names(path) == ["chr1"]


def test_gzipped_reference_is_refused_with_its_name(tmp_path):
	path = write(tmp_path / "ref.fna", GZ)
	with pytest.raises(ValueError, match="ref.fna is gzip-compressed, decompress it first"):
		read_sequences(path)
	with pytest.raises(ValueError, match="gzip-compressed"):
		sequence_names(path)
	problem = check_reference(path, "REF", [("S.1", "chr1", 2, "C", "T")])
	assert "gzip-compressed" in problem and "ref.fna" in problem


def test_convert_names_a_gzipped_or_non_utf8_snp_file(tmp_path):
	ref = write(tmp_path / "REF.fna", FASTA)
	for name, data, words in (("gz.txt", GZ, "gzip-compressed"), ("bad.txt", b"snp_id\tcaf\xe9\n", "not UTF-8")):
		source = write(tmp_path / name, data)
		with pytest.raises(ConversionError, match=name) as info:
			convert_table(source, str(tmp_path / "out" / name), {"REF": ref})
		assert words in str(info.value)
		with pytest.raises(ConversionError, match=name):
			table_genomes(source)
		with pytest.raises(ConversionError, match=name):
			references_from_file(source)


def test_convert_strips_a_bom_from_the_snp_and_references_files(tmp_path):
	ref = write(tmp_path / "REF.fna", FASTA)
	source = write(tmp_path / "snp.txt", b"\xef\xbb\xbfsnp_id\tgenome\tposition\tderived_base\tancestral_base\nA.1\tREF\t2\tT\tC\n")
	assert table_genomes(source) == ["REF"]
	target = str(tmp_path / "out" / "snp.txt")
	assert convert_table(source, target, {"REF": ref}) == 1
	with open(target, "rb") as fin:
		assert fin.read().startswith(b"snp_id\t")
	refs = write(tmp_path / "refs.txt", b"\xef\xbb\xbfgenome\tgenbank_id\tassembly_name\nREF\tGCA_1.1\tASM1\n")
	assert references_from_file(refs) == {"REF": ("GCA_1.1", "ASM1")}


def test_convert_names_a_gzipped_reference(tmp_path):
	path = write(tmp_path / "REF.fna", GZ)
	with pytest.raises(ConversionError, match="REF.fna is gzip-compressed"):
		JoinedReference(path)


def test_mauve_retry_copy_of_a_non_utf8_query(tmp_path):
	query = write(tmp_path / "q.fna", b">q caf\xe9\nAC-GT\n")
	copy = ProgressiveMauve().prepare_retry(query, str(tmp_path))
	with open(copy, "rb") as fin:
		assert fin.read() == b">q caf\xe9\nACNGT\n"


def downloader(tmp_path, payload=None, retrieve=None):
	obj = DownloadGenomes.__new__(DownloadGenomes)
	obj.source, obj.directory, obj.database = "genbank", str(tmp_path / "refs"), None
	os.makedirs(obj.directory, exist_ok=True)
	def fetch(link, target, timeout=None):
		with open(target, "wb") as out:
			out.write(payload)
	obj.retrieve = retrieve or fetch
	return obj


class Results(list):
	def put(self, item):
		self.append(item)


def fetch(obj, refid="REF", source="genbank"):
	results = Results()
	obj.download("GCA_000007845.1", "ASM784v1", refid, source=source, directory=obj.directory, results=results)
	return results[0]


def test_stale_link_is_replaced(tmp_path):
	obj = downloader(tmp_path, GZ)
	os.makedirs(os.path.join(obj.directory, "source"))
	other = write(os.path.join(obj.directory, "source", "other.fna"), FASTA)
	link = os.path.join(obj.directory, "REF.fna")
	os.symlink("source/other.fna", link)
	assert fetch(obj) == ("REF", True, "")
	assert os.readlink(link) == "source/GCA_000007845.1_ASM784v1_genomic.fna"
	assert os.path.exists(other)


def test_dangling_link_is_replaced(tmp_path):
	obj = downloader(tmp_path, GZ)
	link = os.path.join(obj.directory, "REF.fna")
	os.symlink("source/gone.fna", link)
	assert fetch(obj)[1] is True
	assert os.path.exists(link)
	assert os.readlink(link) == "source/GCA_000007845.1_ASM784v1_genomic.fna"


def test_real_file_in_place_of_the_link_is_kept_and_reported(tmp_path):
	obj = downloader(tmp_path, GZ)
	link = write(os.path.join(obj.directory, "REF.fna"), b">x\nAAAA\n")
	refid, ok, reason = fetch(obj)
	assert ok is False and "REF.fna" in reason and "not a link" in reason
	with open(link, "rb") as fin:
		assert fin.read() == b">x\nAAAA\n"


@pytest.mark.parametrize("error", [http.client.IncompleteRead(b"ab", 10), http.client.InvalidURL("bad url"),
	ValueError("bad value"), RuntimeError("unexpected")])
def test_worker_reports_any_download_exception(tmp_path, error):
	def retrieve(link, target, timeout=None):
		raise error
	refid, ok, reason = fetch(downloader(tmp_path, retrieve=retrieve))
	assert ok is False and reason


def test_non_utf8_assembly_report_is_reported(tmp_path):
	def retrieve(link, target, timeout=None):
		with open(target, "wb") as out:
			out.write(GZ if link.endswith(".fna.gz") else b"# caf\xe9\n")
	refid, ok, reason = fetch(downloader(tmp_path, retrieve=retrieve), source="refseq")
	assert ok is True or reason ## reported, not an exception


def open_descriptors():
	count = 0
	for fd in range(256):
		try:
			os.fstat(fd)
			count += 1
		except OSError:
			pass
	return count


def test_download_closes_its_temporary_descriptors(tmp_path):
	obj = downloader(tmp_path, GZ)
	genome = os.path.join(obj.directory, "source", "GCA_000007845.1_ASM784v1_genomic.fna")
	fetch(obj)
	before = open_descriptors()
	for i in range(5):
		os.remove(genome) ## a genome that is present is reused without a download
		fetch(obj, refid="R" + str(i))
	assert open_descriptors() == before


@pytest.mark.parametrize("refid", ["../esc", "a/b", "..", "x\\y"])
def test_genome_name_with_a_path_separator_is_refused(tmp_path, refid):
	obj = downloader(tmp_path, GZ)
	results = Results()
	obj.download("GCA_000007845.1", "ASM784v1", refid, directory=obj.directory, results=results)
	name, ok, reason = results[0]
	assert ok is False and refid in reason
	assert not os.path.exists(str(tmp_path / "esc.fna"))


def test_outdir_that_is_a_file_is_refused_before_any_download(tmp_path):
	from helpers import make_db
	db = make_db(tmp_path / "y.db", [("S.1", "chr1", 2, "C", "T")])
	afile = write(tmp_path / "afile", b"x")
	env = dict(os.environ, PYTHONPATH=os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
	result = subprocess.run([sys.executable, "-m", "CanSNPer2.DownloadGenomes", "-db", db, "-o", afile,
		"--logs", str(tmp_path / "logs")], capture_output=True, text=True, cwd=str(tmp_path), env=env)
	assert result.returncode == 1
	assert "afile" in result.stderr + result.stdout and "not a directory" in result.stderr + result.stdout
	assert "Traceback" not in result.stderr
