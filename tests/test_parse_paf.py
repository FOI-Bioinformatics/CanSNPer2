"""The minimap2 aligner and its PAF reader. Database positions run along the reference
sequences joined end to end; PAF positions are per sequence, so the reader maps them with
the sequence table the aligner writes next to the alignment. The cs tag is on the reference
forward strand, also for "-" alignments, and the query position counts down on those."""
import os
import random
import shutil
import sqlite3

import pytest

from CanSNPer2.modules.Aligners import Minimap2, get_aligner
from CanSNPer2.modules.CanSNPer2 import CanSNPer2
from CanSNPer2.modules.ParsePAF import ParsePAF, write_sequence_table


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


def paf_line(cs, strand="+", tname="chr2", tstart=9, qstart=20, mapq=60, tp="P", qname="ctg1"):
	'''one PAF record; the reference span and the query span follow from the cs string'''
	import re
	ref_len = sum(len(m) for m in re.findall(r"=([A-Za-z]+)", cs)) + cs.count("*") + sum(len(m) for m in re.findall(r"-([a-z]+)", cs))
	query_len = sum(len(m) for m in re.findall(r"=([A-Za-z]+)", cs)) + cs.count("*") + sum(len(m) for m in re.findall(r"\+([a-z]+)", cs))
	fields = [qname, "1000", qstart, qstart + query_len, strand, tname, "50", tstart, tstart + ref_len,
		ref_len, ref_len, mapq, "tp:A:" + tp, "cs:Z:" + cs]
	return "\t".join(str(f) for f in fields)


def read(tmp_path, lines, rows, name="a", table=True, **kwargs):
	'''the reference has chr1 (100 bases) and chr2 (50 bases): chr2 base n is joined 100 + n'''
	db, paf = tmp_path / (name + ".db"), tmp_path / ("REF_" + name + ".paf")
	make_db(db, rows)
	paf.write_text("\n".join(lines) + "\n")
	if table:
		fasta = tmp_path / (name + ".fna")
		fasta.write_text(">chr1 first\n" + "A" * 60 + "\n" + "A" * 40 + "\n>chr2\n" + "C" * 50 + "\n")
		write_sequence_table(str(fasta), str(paf) + ".sequences")
	parser = ParsePAF(**kwargs)
	snps = parser.run(str(paf), reference="REF", database=str(db))
	return snps, {row[0]: row for row in parser.get_snp_info()}


# chr2 from its 10th base (joined 110): ACGT, then joined 114 is c in the reference, g in the query
ROWS = [("B.1", 114, "C", "G")]


def test_sequence_table_maps_to_the_joined_reference(tmp_path):
	snps, info = read(tmp_path, [paf_line("=ACGT*cg=AAAA")], ROWS)
	assert snps == {"B.1": 1}
	assert info["B.1"][:5] == ["B.1", "REF", "114", "ctg1", "25"], "qstart 20 + 4 bases before it, 1-based"
	assert info["B.1"][7] == "G"


def test_query_insertion_and_reverse_strand(tmp_path):
	snps, info = read(tmp_path, [paf_line("=ACGT+tt*cg=AAAA")], ROWS, name="ins")
	assert info["B.1"][4] == "27", "two inserted query bases move the query position, not the reference one"
	# "-": the cs string is still on the reference strand, the query runs down from qend 31
	snps, info = read(tmp_path, [paf_line("=ACGT+tt*cg=AAAA", strand="-")], ROWS, name="rev")
	assert snps == {"B.1": 1}
	assert info["B.1"][4] == "25", "qend 31 - 6 query bases before it"


def test_deletion_at_a_snp_is_no_data(tmp_path):
	snps, info = read(tmp_path, [paf_line("=ACGT-c=AAAA")], ROWS)
	assert snps == {"B.1": 0}
	assert info["B.1"][4] == "N/A" and info["B.1"][7] == "-"


def test_only_primary_alignments_with_enough_mapping_quality_are_read(tmp_path):
	lines = [paf_line("=ACGT*cg=AAAA"),
		paf_line("=ACGT*ca=AAAA", tp="S", qname="paralog"),        ## secondary
		paf_line("=ACGT*ct=AAAA", mapq=3, qname="repeat")]         ## primary, ambiguous placement
	snps, info = read(tmp_path, lines, ROWS)
	assert snps == {"B.1": 1} and info["B.1"][3] == "ctg1"
	snps, info = read(tmp_path, lines, ROWS, name="mapq0", min_mapq=0)
	assert snps == {"B.1": 0} and info["B.1"][7] == "N", "two primary alignments disagree, no evidence either way"


def test_missing_sequence_table_or_cs_tag_is_an_error(tmp_path):
	with pytest.raises(FileNotFoundError, match="sequence table"):
		read(tmp_path, [paf_line("=ACGT*cg=AAAA")], ROWS, table=False)
	with pytest.raises(ValueError, match="cs tag"):
		read(tmp_path, [paf_line("=ACGT*cg=AAAA").rsplit("\t", 1)[0]], ROWS, name="nocs")


def test_minimap2_definition(tmp_path):
	mm2 = get_aligner("minimap2", "/opt/bin/")
	assert isinstance(mm2, Minimap2)
	assert mm2.command("r.fna", "q.fna", "o.paf") == ["/opt/bin/minimap2", "-c", "--cs=long", "-x", "asm20", "-r", "1k,10k", "-o", "o.paf", "r.fna", "q.fna"]
	assert mm2.alignment_path("/tmp/run", "CO92", "s") == "/tmp/run/CO92_s.paf"
	assert isinstance(mm2.reader(), ParsePAF) and mm2.reader().min_mapq == 20
	fasta = tmp_path / "r.fna"
	fasta.write_text(">c1 x\nACGT\nAC\n>c2\nGG\n")
	mm2.prepare(str(fasta), "q.fna", str(tmp_path / "o.paf"))
	assert (tmp_path / "o.paf.sequences").read_text() == "c1\t6\nc2\t2\n"


def revcomp(seq):
	return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


@pytest.mark.skipif(shutil.which("minimap2") is None, reason="minimap2 is not installed")
def test_minimap2_through_the_whole_alignment_step(tmp_path):
	# a two-sequence reference and a query split into three contigs, one of them reverse
	# complemented and the order shuffled; SNPs on both reference sequences
	rng = random.Random(7)
	chr1 = "".join(rng.choice("ACGT") for _ in range(6000))
	chr2 = "".join(rng.choice("ACGT") for _ in range(4000))
	joined = chr1 + chr2
	snps = {1500: "B.1", 4500: "B.2", 6000 + 2000: "B.3"}  ## joined 1-based positions
	rows, query = [], list(joined)
	for pos, snp_id in snps.items():
		ref_base = joined[pos - 1]
		alt = "A" if ref_base != "A" else "C"
		derived = snp_id != "B.2"                 ## B.2 stays ancestral in the query
		rows.append((snp_id, pos, ref_base, alt))
		if derived:
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
		database=db, aligner="minimap2")
	files = obj.align(str(qfile))
	assert files and all(os.path.isfile(f) and os.path.isfile(f + ".sequences") for f in files)
	states, info, called = obj.find_snps_multiproc(xmfa_obj=obj.aligner.reader(database=db), xmfa_files=files,
		references=obj.alignment_references)
	assert states == {"B.1": 1, "B.2": 2, "B.3": 1}
	assert sorted(called) == ["B.1", "B.3"]
	where = {row[0]: (row[3], int(row[4])) for row in info}
	assert where["B.1"] == ("q1", 1500)
	assert where["B.2"] == ("q2", 7000 - 4500 + 1), "q2 is reverse complemented: query 4500 is its base 2501"
	assert where["B.3"] == ("q3", 8000 - 7000)
