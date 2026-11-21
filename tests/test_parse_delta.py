"""The MUMmer4 aligner and its delta reader. Delta positions are per reference sequence, the
coordinates the database stores, and carry no bases: the reader takes the query bases from the
query file named in the first line of the delta file. The reference file is not read."""
import os
import random
import shutil

import pytest

from helpers import make_db
from CanSNPer2.modules.Aligners import Mummer4, get_aligner
from CanSNPer2.modules.CanSNPer2 import CanSNPer2
from CanSNPer2.modules.ParseDelta import ParseDelta


# the reference: chr1 (100 bases) and chr2 (50 bases); the SNP B.1 sits at chr2 position 15,
# where chr2 has C
CHR2 = "AAAAAAAAAA" + "GGGGCAAAAA" + "A" * 30
ROWS = [("B.1", "chr2", 15, "C", "T")]


def read(tmp_path, query, record, rows=ROWS, name="a", records=None, keep_reference=True):
	'''record: (s1, e1, s2, e2, [indels]) of one alignment of chr2 against the query ctg1;
		records: (reference sequence, s1, e1, s2, e2, [indels]) against ctg1, in place of record'''
	ref = tmp_path / (name + "_ref.fna")
	ref.write_text(">chr1 first\n" + "A" * 60 + "\n" + "A" * 40 + "\n>chr2\n" + CHR2 + "\n")
	qry = tmp_path / (name + "_qry.fna")
	qry.write_text(">ctg0\nCCCC\n>ctg1 the one\n" + query + "\n")
	if records is None:
		records = [("chr2",) + tuple(record)]
	lengths = {"chr1": 100, "chr2": 50}
	text = "{0} {1}\nNUCMER\n".format(ref, qry)
	for rname, s1, e1, s2, e2, indels in records:
		text += ">{0} ctg1 {1} {2}\n{3} {4} {5} {6} 1 1 0\n{7}0\n".format(
			rname, lengths[rname], len(query), s1, e1, s2, e2, "".join("{0}\n".format(d) for d in indels))
	delta = tmp_path / ("REF_" + name + ".delta")
	delta.write_text(text)
	if not keep_reference:
		ref.unlink()
	db = tmp_path / (name + ".db")
	make_db(db, rows)
	parser = ParseDelta()
	snps = parser.run(str(delta), reference="REF", database=str(db))
	return snps, {row[0]: row for row in parser.get_snp_info()}


def test_forward_alignment_in_reference_sequence_coordinates(tmp_path):
	# chr2 11..20 = GGGGCAAAAA against ctg1 1..10 = GGGGTAAAAA: the SNP is the 5th pair
	snps, info = read(tmp_path, "GGGGTAAAAA", (11, 20, 1, 10, []))
	assert snps == {"B.1": 1}
	assert info["B.1"][1:6] == ["REF", "chr2", "15", "ctg1", "5"] and info["B.1"][8] == "T"


def test_delta_reads_reference_coordinates_without_the_reference_file(tmp_path):
	# the reference file named in the delta header is gone: the query base is still found
	snps, info = read(tmp_path, "GGGGTAAAAA", (11, 20, 1, 10, []), name="noref", keep_reference=False)
	assert not (tmp_path / "noref_ref.fna").exists()
	assert snps == {"B.1": 1}
	assert info["B.1"][2:6] == ["chr2", "15", "ctg1", "5"]


def test_same_position_on_two_sequences_both_called(tmp_path):
	# chr1 41..50 against ctg1 1..10 and chr2 41..50 against ctg1 11..20; position 50 on chr1
	# and position 50 on chr2 are different SNPs
	rows = [("A.1", "chr1", 50, "A", "T"), ("B.2", "chr2", 50, "A", "G")]
	snps, info = read(tmp_path, "AAAAAAAAAT" + "AAAAAAAAAG", None, rows=rows, name="two",
		records=[("chr1", 41, 50, 1, 10, []), ("chr2", 41, 50, 11, 20, [])])
	assert snps == {"A.1": 1, "B.2": 1}
	assert info["A.1"][2:6] == ["chr1", "50", "ctg1", "10"] and info["A.1"][8] == "T"
	assert info["B.2"][2:6] == ["chr2", "50", "ctg1", "20"] and info["B.2"][8] == "G"


def test_query_insertion_moves_only_the_query_position(tmp_path):
	# -3: two pairs, then a query base with a reference gap; the SNP pairs with query 6
	snps, info = read(tmp_path, "GGXGGTAAAAA".replace("X", "C"), (11, 20, 1, 11, [-3]), name="ins")
	assert snps == {"B.1": 1} and info["B.1"][5] == "6"


def test_deletion_at_the_snp_is_no_data(tmp_path):
	# +5: four pairs, then reference position 15 with a gap in the query
	snps, info = read(tmp_path, "GGGGAAAAA", (11, 20, 1, 9, [5]), name="del")
	assert snps == {"B.1": 0} and info["B.1"][5] == "N/A" and info["B.1"][8] == "-"


def test_reverse_alignment_is_complemented(tmp_path):
	# S2 > E2: chr2 11 pairs with ctg1 10, so chr2 15 pairs with ctg1 6, whose base A is T on
	# the reference strand
	snps, info = read(tmp_path, "TTTTTACCCC", (11, 20, 10, 1, []), name="rev")
	assert snps == {"B.1": 1}
	assert info["B.1"][5] == "6" and info["B.1"][8] == "T"


def test_missing_query_file_is_an_error(tmp_path):
	delta = tmp_path / "REF_x.delta"
	delta.write_text("{0}/gone_ref.fna {0}/gone_qry.fna\nNUCMER\n".format(tmp_path))
	db = tmp_path / "x.db"
	make_db(db, ROWS)
	with pytest.raises(FileNotFoundError, match="query sequence file named in"):
		ParseDelta().run(str(delta), reference="REF", database=str(db))


def test_mummer4_definition(tmp_path):
	mum = get_aligner("mummer4", "/opt/mummer/bin/")
	assert isinstance(mum, Mummer4)
	command = mum.command("refs/r.fna", "q.fna", "/tmp/run/R_q.delta")
	assert command[:2] == ["sh", "-c"]
	script = command[2]
	assert script.startswith("/opt/mummer/bin/nucmer -t 1 --delta=/tmp/run/R_q.delta.unfiltered " + os.path.abspath("refs/r.fna"))
	assert "/opt/mummer/bin/delta-filter -1 /tmp/run/R_q.delta.unfiltered > /tmp/run/R_q.delta" in script
	assert mum.alignment_path("/tmp/run", "CO92", "s") == "/tmp/run/CO92_s.delta"
	assert isinstance(mum.reader(), ParseDelta)
	assert Mummer4().filter_program() == "delta-filter"


def revcomp(seq):
	return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


@pytest.mark.skipif(shutil.which("nucmer") is None or shutil.which("delta-filter") is None, reason="MUMmer4 is not installed")
def test_mummer4_through_the_whole_alignment_step(tmp_path):
	# a two-sequence reference and a query split into three contigs, one of them reverse
	# complemented and the order shuffled; SNPs on both reference sequences
	rng = random.Random(7)
	chr1 = "".join(rng.choice("ACGT") for _ in range(6000))
	chr2 = "".join(rng.choice("ACGT") for _ in range(4000))
	joined = chr1 + chr2
	## the query is built from the joined sequence; the SNPs are stored per reference sequence
	snps = {1500: ("B.1", "chr1", 1500), 4500: ("B.2", "chr1", 4500), 6000 + 2000: ("B.3", "chr2", 2000)}
	rows, query = [], list(joined)
	for pos, (snp_id, sequence, seq_pos) in snps.items():
		ref_base = joined[pos - 1]
		alt = "A" if ref_base != "A" else "C"
		rows.append((snp_id, sequence, seq_pos, ref_base, alt))
		if snp_id != "B.2":
			query[pos - 1] = alt
	query = "".join(query)
	contigs = [("q1", query[:3000]), ("q2", revcomp(query[3000:7000])), ("q3", query[7000:])]
	refdir = tmp_path / "refs"
	refdir.mkdir()
	(refdir / "REF.fna").write_text(">chr1\n{0}\n>chr2\n{1}\n".format(chr1, chr2))
	qfile = tmp_path / "sample.fna"
	qfile.write_text("".join(">{0}\n{1}\n".format(name, seq) for name, seq in (contigs[2], contigs[0], contigs[1])))
	db = str(tmp_path / "t.db")
	make_db(db, rows)
	obj = CanSNPer2([str(qfile)], workdir=str(tmp_path), refdir="refs", tmpdir="shared", logdir="logs", outdir="out",
		database=db, aligner="mummer4")
	files = obj.align(str(qfile))
	assert files and all(os.path.isfile(f) for f in files)
	assert not [f for f in os.listdir(obj.tmpdir) if f.endswith(".unfiltered")], "the unfiltered delta is removed"
	states, info, called = obj.find_snps_multiproc(xmfa_obj=obj.aligner.reader(database=db), xmfa_files=files,
		references=obj.alignment_references)
	assert states == {"B.1": 1, "B.2": 2, "B.3": 1}
	assert {row[0]: (row[2], row[3]) for row in info} == {"B.1": ("chr1", "1500"), "B.2": ("chr1", "4500"), "B.3": ("chr2", "2000")}
	where = {row[0]: (row[4], int(row[5])) for row in info}
	assert where["B.1"] == ("q1", 1500)
	assert where["B.2"] == ("q2", 7000 - 4500 + 1)
	assert where["B.3"] == ("q3", 8000 - 7000)


def test_query_path_with_a_space(tmp_path):
	# the header line holds two paths; the query path contains a space (and so does the reference path)
	snps, info = read(tmp_path, "GGGGTAAAAA", (11, 20, 1, 10, []), name="with space")
	assert snps == {"B.1": 1}
	assert info["B.1"][4:6] == ["ctg1", "5"]


def test_query_path_from_the_query_option_when_the_header_is_ambiguous(tmp_path):
	# the header names a query file that does not exist; the query file the caller gave is used
	ref = tmp_path / "ref.fna"
	ref.write_text(">chr2\n" + CHR2 + "\n")
	qry = tmp_path / "real query.fna"
	qry.write_text(">ctg1\nGGGGTAAAAA\n")
	delta = tmp_path / "REF_q.delta"
	delta.write_text("{0} /gone/some other/qry.fna\nNUCMER\n>chr2 ctg1 50 10\n11 20 1 10 1 1 0\n0\n".format(ref))
	db = tmp_path / "q.db"
	make_db(db, ROWS)
	parser = ParseDelta(query=str(qry))
	assert parser.run(str(delta), reference="REF", database=str(db)) == {"B.1": 1}


def test_query_name_missing_from_the_fasta_is_a_clear_error(tmp_path):
	ref = tmp_path / "ref.fna"
	ref.write_text(">chr2\n" + CHR2 + "\n")
	qry = tmp_path / "qry.fna"
	qry.write_text(">other\nGGGGTAAAAA\n")
	delta = tmp_path / "REF_m.delta"
	delta.write_text("{0} {1}\nNUCMER\n>chr2 ctg1 50 10\n11 20 1 10 1 1 0\n0\n".format(ref, qry))
	db = tmp_path / "m.db"
	make_db(db, ROWS)
	with pytest.raises(ValueError) as e:
		ParseDelta().run(str(delta), reference="REF", database=str(db))
	assert "ctg1" in str(e.value) and "qry.fna" in str(e.value)


def test_position_past_the_end_of_the_query_is_no_data(tmp_path):
	# the query has 3 bases, the alignment claims 10: the SNP at the 5th pair is not covered
	snps, info = read(tmp_path, "GGG", (11, 20, 1, 10, []), name="short")
	assert snps == {} and info == {}


def test_negative_query_position_does_not_read_the_last_base(tmp_path):
	# reverse alignment running below position 1 of the query; the SNP (5th pair) would be at
	# query position 0, which must not be read as the last base of the query
	snps, info = read(tmp_path, "GGGGTAAAAT", (11, 20, 4, -5, []), name="neg")
	assert snps == {} and info == {}
