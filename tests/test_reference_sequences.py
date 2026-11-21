"""References of several sequences (chromosomes, plasmids).

Database positions run along the reference sequences joined end to end in file order. The
database does not record the sequences, and the reference check only compared the total
length, so a reference with the same sequences in another order (or another genome of
about the same length) passed the check. Every aligner then read the bases of the wrong
positions: on Brucella (two chromosomes) most calls became NA without an error that
pointed at the reference.

- The reference base at each SNP position must be the ancestral or the derived base of
  that SNP. A reference that fails this at many SNPs is refused, a few mismatches are
  reported as a warning.
- The result tables give the reference sequence and the position on that sequence next
  to the position on the joined reference, so the positions can be compared with
  publications that give them per chromosome."""
import logging
import random
import sqlite3

import pytest

from CanSNPer2.DownloadGenomes import DownloadGenomes
from CanSNPer2.modules.CanSNPer2 import CanSNPer2, CanSNPer2Error, SNP_FILE_COLUMNS, add_reference_sequence
from CanSNPer2.modules.DatabaseConnection import CanSNPdbFunctions
from CanSNPer2.modules.ReferenceCheck import base_mismatches, check_reference, locate, reference_snps, sequence_table

rng = random.Random(7)
CHR1 = "".join(rng.choice("ACGT") for _ in range(600))
CHR2 = "".join(rng.choice("ACGT") for _ in range(400))
JOINED = CHR1 + CHR2
POSITIONS = [50, 120, 333, 480, 599, 601, 700, 850, 999, 1000]   ## 5 on chr1, 5 on chr2
OTHER = {"A": "G", "C": "T", "G": "A", "T": "C"}


def fasta(*records):
	return "".join(">{0} test\n{1}\n".format(name, seq) for name, seq in records)


def snp_rows():
	'''(snp_id, position, ancestral, derived): the reference carries the ancestral base, the
		derived base is another one'''
	return [("S.{0}".format(i), pos, JOINED[pos - 1], OTHER[JOINED[pos - 1]]) for i, pos in enumerate(POSITIONS)]


def make_db(path, rows=None):
	conn = sqlite3.connect(path)
	conn.executescript("""
		CREATE TABLE nodes (id integer PRIMARY KEY, name text NOT NULL);
		CREATE TABLE tree (parent integer, child integer, rank_i integer);
		CREATE TABLE snp_annotation (node_id integer PRIMARY KEY, snp_id varchar(6), position integer,
			ancestral_base varchar(1), derived_base varchar(1), reference varchar(20),
			date datetime, genome_i integer);
		CREATE TABLE snp_references (id integer PRIMARY KEY, genome varchar(6), strain varchar(6));
		INSERT INTO snp_references (id, genome, strain) VALUES (5, 'REF', 'REF');
	""")
	for i, (snp_id, pos, anc, der) in enumerate(rows or snp_rows()):
		conn.execute("INSERT INTO nodes (id, name) VALUES (?, ?)", (i + 1, snp_id))
		conn.execute("INSERT INTO snp_annotation VALUES (?, ?, ?, ?, ?, 'test', '06/11/26', 5)", (i + 1, snp_id, pos, anc, der))
	conn.commit()
	conn.close()
	return path


@pytest.fixture
def references(tmp_path):
	right = tmp_path / "right.fna"
	right.write_text(fasta(("chr1", CHR1), ("chr2", CHR2)))
	swapped = tmp_path / "swapped.fna"
	swapped.write_text(fasta(("chr2", CHR2), ("chr1", CHR1)))
	return str(right), str(swapped)


def test_the_sequence_table_gives_names_and_lengths_in_file_order(references):
	right, swapped = references
	assert sequence_table(right) == [("chr1", 600), ("chr2", 400)]
	assert sequence_table(swapped) == [("chr2", 400), ("chr1", 600)]


def test_a_joined_position_is_located_on_its_sequence(references):
	table = sequence_table(references[0])
	assert locate(table, 1) == ("chr1", 1)
	assert locate(table, 600) == ("chr1", 600)
	assert locate(table, 601) == ("chr2", 1)
	assert locate(table, 1000) == ("chr2", 400)
	assert locate(table, 1001) == ("-", "-")
	assert locate(table, "-") == ("-", "-")


def test_the_snps_of_each_reference_are_read_from_the_database(tmp_path):
	db = make_db(str(tmp_path / "t.db"))
	snps = reference_snps(CanSNPdbFunctions(db))
	assert sorted(snps) == ["REF"]
	assert sorted(snps["REF"]) == sorted(snp_rows())


def test_the_right_reference_has_the_ancestral_or_derived_base_at_every_snp(references):
	right, swapped = references
	assert base_mismatches(right, snp_rows()) == []
	assert check_reference(right, "REF", 1000, snp_rows()) is None


def test_a_reference_with_its_sequences_in_another_order_is_refused(references):
	'''same sequences, same total length: the length check alone let this through'''
	right, swapped = references
	assert check_reference(swapped, "REF", 1000) is None
	mismatches = base_mismatches(swapped, snp_rows())
	assert len(mismatches) >= 5
	problem = check_reference(swapped, "REF", 1000, snp_rows())
	assert problem is not None
	assert "chr2, chr1" in problem and "order" in problem
	assert "{n} of 10 SNPs".format(n=len(mismatches)) in problem


def test_a_single_mismatching_snp_is_a_warning_only(references, caplog):
	'''one wrong base in a scheme is a scheme error to report, not a reason to stop typing'''
	rows = snp_rows()
	snp_id, pos, anc, der = rows[3]
	neither = [b for b in "ACGT" if b not in (anc, der)]
	rows[3] = (snp_id, pos, neither[0], neither[1])
	with caplog.at_level(logging.WARNING):
		assert check_reference(references[0], "REF", 1000, rows) is None
	assert "S.3" in caplog.text and "1 of 10 SNPs" in caplog.text


def test_typing_stops_on_a_reference_with_its_sequences_in_another_order(tmp_path, references):
	db = make_db(str(tmp_path / "t.db"))
	refdir = tmp_path / "refs"
	refdir.mkdir()
	(refdir / "REF.fna").write_text(open(references[1]).read())
	query = tmp_path / "sample.fna"
	query.write_text(fasta(("chr1", CHR1), ("chr2", CHR2)))
	obj = CanSNPer2([str(query)], workdir=str(tmp_path), refdir="refs", tmpdir="shared", logdir="logs", outdir="out", database=db)
	with pytest.raises(CanSNPer2Error, match="reference base"):
		obj.align(str(query))


class Results(list):
	def put(self, item):
		self.append(item)


def test_a_download_with_its_sequences_in_another_order_is_refused(tmp_path, references):
	import gzip
	obj = DownloadGenomes.__new__(DownloadGenomes)
	obj.source, obj.directory, obj.database = "genbank", str(tmp_path / "refs"), None
	data = gzip.compress(open(references[1]).read().encode())
	def retrieve(link, target, timeout=None):
		with open(target, "wb") as out:
			out.write(data)
	obj.retrieve = retrieve
	results = Results()
	obj.download("GCA_000007125.1", "ASM712v1", "REF", directory=obj.directory, results=results, last_position=1000, snps=snp_rows())
	assert results[0][1] is False and "reference base" in results[0][2]


def test_the_result_rows_carry_the_reference_sequence_and_its_position(references):
	assert SNP_FILE_COLUMNS[-2:] == ["ReferenceSeq", "ReferenceSeqPosition"]
	tables = {"REF": sequence_table(references[0])}
	rows = [["S.6", "REF", "700", "q2", "100", "A", "G", "A", ""],
			["S.0", "REF", "50", "q1", "50", "A", "G", "A", "edge"],
			["X", "-", "-", "N/A", "N/A", "-", "-", "N/A", ""],
			["S.9", "OTHER", "5", "N/A", "N/A", "A", "G", "N/A", ""]]
	out = add_reference_sequence(rows, tables)
	assert out[0][-2:] == ["chr2", "100"]
	assert out[1][-2:] == ["chr1", "50"]
	assert out[2][-2:] == ["-", "-"]
	assert out[3][-2:] == ["-", "-"]          ## reference file not at hand (--skip_alignment)
	assert all(len(row) == len(SNP_FILE_COLUMNS) for row in out)
	assert rows[0][-1] == ""                  ## the input rows are not changed
