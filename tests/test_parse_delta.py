"""The MUMmer4 aligner and its delta reader. Delta positions are per sequence and carry no
bases: the reader takes the bases from the query file and the reference sequence order and
lengths from the reference file, both named in the first line of the delta file."""
import os
import random
import shutil
import sqlite3

import pytest

from CanSNPer2.modules.Aligners import Mummer4, get_aligner
from CanSNPer2.modules.CanSNPer2 import CanSNPer2
from CanSNPer2.modules.ParseDelta import ParseDelta


def make_db(path, rows, genome="REF"):
	'''rows: (snp_id, position, ancestral_base, derived_base) against reference genome'''
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
	for i, (snp_id, pos, anc, der) in enumerate(rows):
		conn.execute("INSERT INTO nodes (id, name) VALUES (?, ?)", (i + 1, snp_id))
		conn.execute("INSERT INTO snp_annotation (node_id, snp_id, position, ancestral_base, derived_base, reference, date, genome_i)"
			" VALUES (?, ?, ?, ?, ?, 'test', '06/11/26', 5)", (i + 1, snp_id, pos, anc, der))
	conn.commit()
	conn.close()


# the reference: chr1 (100 bases) and chr2 (50 bases), so chr2 base n is joined 100 + n;
# the SNP B.1 sits at chr2 position 15, joined 115, where chr2 has C
CHR2 = "AAAAAAAAAA" + "GGGGCAAAAA" + "A" * 30
ROWS = [("B.1", 115, "C", "T")]


def read(tmp_path, query, record, rows=ROWS, name="a"):
	'''record: (s1, e1, s2, e2, [indels]) of one alignment of chr2 against the query ctg1'''
	ref = tmp_path / (name + "_ref.fna")
	ref.write_text(">chr1 first\n" + "A" * 60 + "\n" + "A" * 40 + "\n>chr2\n" + CHR2 + "\n")
	qry = tmp_path / (name + "_qry.fna")
	qry.write_text(">ctg0\nCCCC\n>ctg1 the one\n" + query + "\n")
	s1, e1, s2, e2, indels = record
	delta = tmp_path / ("REF_" + name + ".delta")
	delta.write_text("{0} {1}\nNUCMER\n>chr2 ctg1 50 {2}\n{3} {4} {5} {6} 1 1 0\n{7}0\n".format(
		ref, qry, len(query), s1, e1, s2, e2, "".join("{0}\n".format(d) for d in indels)))
	db = tmp_path / (name + ".db")
	make_db(db, rows)
	parser = ParseDelta()
	snps = parser.run(str(delta), reference="REF", database=str(db))
	return snps, {row[0]: row for row in parser.get_snp_info()}


def test_forward_alignment_in_joined_coordinates(tmp_path):
	# chr2 11..20 = GGGGCAAAAA against ctg1 1..10 = GGGGTAAAAA: the SNP is the 5th pair
	snps, info = read(tmp_path, "GGGGTAAAAA", (11, 20, 1, 10, []))
	assert snps == {"B.1": 1}
	assert info["B.1"][1:5] == ["REF", "115", "ctg1", "5"] and info["B.1"][7] == "T"


def test_query_insertion_moves_only_the_query_position(tmp_path):
	# -3: two pairs, then a query base with a reference gap; the SNP pairs with query 6
	snps, info = read(tmp_path, "GGXGGTAAAAA".replace("X", "C"), (11, 20, 1, 11, [-3]), name="ins")
	assert snps == {"B.1": 1} and info["B.1"][4] == "6"


def test_deletion_at_the_snp_is_no_data(tmp_path):
	# +5: four pairs, then reference position 15 with a gap in the query
	snps, info = read(tmp_path, "GGGGAAAAA", (11, 20, 1, 9, [5]), name="del")
	assert snps == {"B.1": 0} and info["B.1"][4] == "N/A" and info["B.1"][7] == "-"


def test_reverse_alignment_is_complemented(tmp_path):
	# S2 > E2: chr2 11 pairs with ctg1 10, so chr2 15 pairs with ctg1 6, whose base A is T on
	# the reference strand
	snps, info = read(tmp_path, "TTTTTACCCC", (11, 20, 10, 1, []), name="rev")
	assert snps == {"B.1": 1}
	assert info["B.1"][4] == "6" and info["B.1"][7] == "T"


def test_missing_sequence_files_are_an_error(tmp_path):
	delta = tmp_path / "REF_x.delta"
	delta.write_text("{0}/gone_ref.fna {0}/gone_qry.fna\nNUCMER\n".format(tmp_path))
	db = tmp_path / "x.db"
	make_db(db, ROWS)
	with pytest.raises(FileNotFoundError, match="sequence files named in"):
		ParseDelta().run(str(delta), reference="REF", database=str(db))


def test_mummer4_definition(tmp_path):
	mum = get_aligner("mummer4", "/opt/mummer/bin/")
	assert isinstance(mum, Mummer4)
	command = mum.command("refs/r.fna", "q.fna", "/tmp/run/R_q.delta")
	assert command[:2] == ["sh", "-c"]
	script = command[2]
	assert script.startswith("/opt/mummer/bin/nucmer --delta=/tmp/run/R_q.delta.unfiltered " + os.path.abspath("refs/r.fna"))
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
	snps = {1500: "B.1", 4500: "B.2", 6000 + 2000: "B.3"}
	rows, query = [], list(joined)
	for pos, snp_id in snps.items():
		ref_base = joined[pos - 1]
		alt = "A" if ref_base != "A" else "C"
		rows.append((snp_id, pos, ref_base, alt))
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
	where = {row[0]: (row[3], int(row[4])) for row in info}
	assert where["B.1"] == ("q1", 1500)
	assert where["B.2"] == ("q2", 7000 - 4500 + 1)
	assert where["B.3"] == ("q3", 8000 - 7000)
