"""Regression tests for issue #39: SNPs added to the database as custom SNPs
were not called because the scan of an alignment block could stop early.

Each test feeds ParseXMFA a small database (SNP rows inserted the way
ModifyCanSNPer2Database.add_annotation does it) and a hand-made xmfa, and
checks which SNPs end up evaluated and called."""
import os
import random
import shutil

import pytest

from helpers import make_db
from CanSNPer2.modules.AlignmentReader import SEQUENCES_SUFFIX, read_sequence_table
from CanSNPer2.modules.Aligners import ProgressiveMauve
from CanSNPer2.modules.CanSNPer2 import CanSNPer2
from CanSNPer2.modules.ParseXMFA import ParseXMFA

## the single reference sequence of the hand-made alignments, longer than every block
SEQUENCE = "chr"
SEQUENCE_LENGTH = 100000


def write_table(xmfa, sequences=((SEQUENCE, SEQUENCE_LENGTH),)):
	'''the reference sequence table that ProgressiveMauve.prepare writes next to an alignment'''
	with open(str(xmfa) + SEQUENCES_SUFFIX, "w") as out:
		for name, length in sequences:
			out.write("{0}\t{1}\n".format(name, length))


def on_sequence(rows):
	'''(snp_id, position, ancestral, derived) rows on the single reference sequence'''
	return [(snp_id, SEQUENCE, pos, anc, der) for snp_id, pos, anc, der in rows]


def make_xmfa(path, blocks):
	'''blocks: (start, end, ref_seq, target_seq) or (start, end, ref_seq, target_seq, sign)
		using the mauve xmfa layout'''
	with open(path, "w") as f:
		f.write("#FormatVersion Mauve1\n")
		for block in blocks:
			start, end, ref, target = block[:4]
			sign = block[4] if len(block) > 4 else "+"
			f.write("============\n> 1:{s}-{e} {sign} ref\n{r}\n> 1:{s}-{e} {sign} query\n{t}\n".format(s=start, e=end, sign=sign, r=ref, t=target))


def pair(length, overrides):
	'''identical ref/target strings, overrides: column -> (ref_base, target_base)'''
	ref = list("C" * length)
	target = list("C" * length)
	for col, (r, t) in overrides.items():
		ref[col] = r
		target[col] = t
	return "".join(ref), "".join(target)


def run(tmp_path, rows, blocks, name="t", sequences=((SEQUENCE, SEQUENCE_LENGTH),), table=True):
	'''rows: (snp_id, position, ancestral, derived) on the single reference sequence, or
		(snp_id, sequence, position, ancestral, derived) with sequences given'''
	db, xmfa = tmp_path / (name + ".db"), tmp_path / ("REF_" + name + ".xmfa")
	make_db(db, rows if rows and len(rows[0]) == 5 else on_sequence(rows))
	make_xmfa(xmfa, blocks)
	if table:
		write_table(xmfa, sequences)
	parser = ParseXMFA()
	snps = parser.run(xmfa=str(xmfa), reference="REF", database=str(db))
	return snps, parser.get_called_snps(), parser


def test_custom_snp_at_block_start_is_called(tmp_path):
	# the next SNP sits exactly at the block start: the whole block used to be
	# skipped by the strict range check, dropping every SNP of the block (#39)
	rows = [("C.1", 10, "A", "T"), ("B.1", 25, "C", "T"), ("C.2", 40, "G", "T")]
	ref, target = pair(31, {0: ("A", "T"), 15: ("C", "T"), 30: ("G", "T")})
	snps, called, parser = run(tmp_path, rows, [(10, 40, ref, target)])
	assert called == ["C.1", "B.1", "C.2"]
	assert snps == {"C.1": 1, "B.1": 1, "C.2": 1}


def test_snp_beyond_current_block_is_found_in_a_later_block(tmp_path):
	# advancing past the end of one block must not stop the scan of the next one
	rows = [("B.1", 15, "C", "T"), ("C.gap", 25, "A", "T"), ("C.1", 45, "G", "T")]
	ref1, target1 = pair(11, {5: ("C", "T")})
	ref2, target2 = pair(21, {15: ("G", "T")})
	snps, called, parser = run(tmp_path, rows, [(10, 20, ref1, target1), (30, 50, ref2, target2)])
	assert called == ["B.1", "C.1"]
	assert "C.gap" not in snps ## unaligned region, the SNP was never comparable

def test_block_without_snps_does_not_end_the_scan(tmp_path):
	# a block before the first registered SNP is skipped, later SNPs still called
	rows = [("C.1", 45, "G", "T")]
	ref1, target1 = pair(11, {})
	ref2, target2 = pair(21, {15: ("G", "T")})
	snps, called, parser = run(tmp_path, rows, [(10, 20, ref1, target1), (30, 50, ref2, target2)])
	assert called == ["C.1"]

def test_short_target_only_skips_the_uncovered_tail(tmp_path):
	# issue #97 break: the target string ends before the block does, SNPs in the
	# covered part are still evaluated, SNPs in the uncovered tail are not
	rows = [("B.1", 15, "C", "T"), ("C.1", 38, "G", "T")]
	ref, target = pair(31, {5: ("C", "T"), 28: ("G", "T")})
	snps, called, parser = run(tmp_path, rows, [(10, 40, ref, target[:20])])
	assert called == ["B.1"]
	assert "C.1" not in snps


def test_snps_at_the_alignment_edges_are_flagged_but_still_called(tmp_path):
	# issue #48: a SNP sitting within EDGE_MARGIN bases of either end of an alignment block
	# is the most likely misalignment, it must be evaluated AND flagged "edge" in the info
	# rows so it can be reviewed, not silently trusted (mask) and not silently dropped
	rows = [("B.1", 10, "C", "T"), ("M.1", 100, "C", "T"), ("C.1", 180, "G", "T")]
	ref, target = pair(171, {0: ("C", "T"), 90: ("C", "T"), 170: ("G", "T")})
	snps, called, parser = run(tmp_path, rows, [(10, 180, ref, target)])
	assert called == ["B.1", "M.1", "C.1"], "the edge SNPs are still called"
	flags = {row[0]: row[9] for row in parser.get_snp_info()}
	assert flags == {"B.1": "edge", "M.1": "", "C.1": "edge"}, flags
	assert parser.edge_snps == {"B.1", "C.1"}


def test_blocks_out_of_reference_order_are_all_read(tmp_path):
	# the SNP list used to be consumed as a queue while walking the blocks, so a SNP in a
	# block that came after a block further along the reference was never evaluated and was
	# missing from the result entirely; positions are looked up per block now
	rows = [("B.1", 15, "C", "T"), ("C.1", 45, "G", "T")]
	ref1, target1 = pair(11, {5: ("C", "T")})
	ref2, target2 = pair(21, {15: ("G", "T")})
	snps, called, parser = run(tmp_path, rows, [(30, 50, ref2, target2), (10, 20, ref1, target1)])
	assert snps == {"B.1": 1, "C.1": 1}
	assert called == ["B.1", "C.1"], "rows and calls come in reference position order"
	assert [row[0] for row in parser.get_snp_info()] == ["B.1", "C.1"]


def test_reverse_block_reports_bases_on_the_reference_strand(tmp_path):
	# a "-" block shows the reverse complement of reference 10..20: position p sits in
	# column 20 - p, and the query base there must be complemented before it is compared
	rows = [("B.1", 12, "C", "T"), ("B.2", 18, "A", "G")]
	ref, target = pair(11, {20 - 12: ("G", "A"), 20 - 18: ("T", "T")})
	snps, called, parser = run(tmp_path, rows, [(10, 20, ref, target, "-")])
	assert snps == {"B.1": 1, "B.2": 2}
	info = {row[0]: row for row in parser.get_snp_info()}
	assert info["B.1"][8] == "T" and info["B.2"][8] == "A"
	assert [row[0] for row in parser.get_snp_info()] == ["B.1", "B.2"]


def test_position_covered_twice_with_different_bases_is_no_data(tmp_path, caplog):
	# progressiveMauve never covers a reference position twice, other aligners can; two
	# different query bases at one SNP are no evidence either way
	rows = [("B.1", 15, "C", "T")]
	ref1, target1 = pair(11, {5: ("C", "T")})
	ref2, target2 = pair(11, {5: ("C", "C")})
	snps, called, parser = run(tmp_path, rows, [(10, 20, ref1, target1), (10, 20, ref2, target2)])
	assert snps == {"B.1": 0}
	assert called == []
	assert [row[8] for row in parser.get_snp_info()] == ["N"]
	assert "more than one alignment" in caplog.text
	# the same base twice is one observation
	snps, called, parser = run(tmp_path, rows, [(10, 20, ref1, target1), (10, 20, ref1, target1)], name="same")
	assert snps == {"B.1": 1} and len(parser.get_snp_info()) == 1


QUERY_FASTA = ">ctgA first contig\n" + "A" * 30 + "\n>ctgB\n" + "A" * 20 + "\n" + "A" * 20 + "\n"


def query_case(tmp_path, ref_block, query_block, rows, name, header_file=None, query_file=True):
	'''one block: ref_block (start, end, sign, seq), query_block (start, end, sign, seq); the
		query fasta has ctgA (30 bases) and ctgB (40 bases), so ctgB starts at joined 31'''
	db, xmfa = tmp_path / (name + ".db"), tmp_path / ("REF_" + name + ".xmfa")
	fasta = tmp_path / (name + ".fasta")
	fasta.write_text(QUERY_FASTA)
	make_db(db, on_sequence(rows))
	write_table(xmfa)
	with open(xmfa, "w") as f:
		f.write("#FormatVersion Mauve1\n")
		if header_file:
			f.write("#Sequence2File\t{0}\n#Sequence2Format\tFastA\n".format(header_file if header_file is not True else fasta))
		f.write("> 1:{0}-{1} {2} ref.fna\n{3}\n> 2:{4}-{5} {6} {7}\n{8}\n=\n".format(
			*ref_block[:3], ref_block[3], *query_block[:3], fasta, query_block[3]))
	parser = ParseXMFA(query=str(fasta) if query_file else None)
	parser.run(xmfa=str(xmfa), reference="REF", database=str(db))
	return {row[0]: row for row in parser.get_snp_info()}, parser


def test_query_contig_and_position_are_reported(tmp_path):
	# issue #23: every called position must say where it sits in the submitted assembly. The
	# xmfa coordinates of the query run along all its contigs joined end to end, they are
	# mapped to the contig name and the position in that contig; target gaps are skipped
	rows = [("B.1", 15, "C", "T")]
	ref, target = pair(11, {5: ("C", "T")})
	target = target[:3] + "-" + target[4:]              ## a gap before the SNP
	info, _ = query_case(tmp_path, (10, 20, "+", ref), (41, 50, "+", target), rows, "plus")
	assert info["B.1"][4:6] == ["ctgB", "15"], "joined 41 + 4 query bases before the SNP = 45 = ctgB 15"
	assert info["B.1"][8] == "T"
	info, _ = query_case(tmp_path, (10, 20, "+", ref), (5, 14, "+", target), rows, "first")
	assert info["B.1"][4:6] == ["ctgA", "9"]


def test_query_position_on_every_strand_combination(tmp_path):
	# the walk follows the reference strand and the query row follows the query strand;
	# the position used to be right only when both strands were the same
	rows = [("B.1", 12, "C", "T")]
	# reference "-": reference position 12 sits in column 20 - 12 = 8 of the reversed block,
	# the query "+" row has 8 query bases left of it: joined 41 + 8 = 49 = ctgB 19
	ref, target = pair(11, {8: ("G", "A")})
	info, _ = query_case(tmp_path, (10, 20, "-", ref), (41, 51, "+", target), rows, "rminus_qplus")
	assert info["B.1"][4:6] == ["ctgB", "19"] and info["B.1"][8] == "T"
	# reference "+", query "-": column 2, the "-" query row runs down from 51: 51 - 2 = 49
	ref, target = pair(11, {2: ("C", "T")})
	info, _ = query_case(tmp_path, (10, 20, "+", ref), (41, 51, "-", target), rows, "rplus_qminus")
	assert info["B.1"][4:6] == ["ctgB", "19"] and info["B.1"][8] == "T"
	# both "-": column 8, 51 - 8 = 43 = ctgB 13
	ref, target = pair(11, {8: ("G", "A")})
	info, _ = query_case(tmp_path, (10, 20, "-", ref), (41, 51, "-", target), rows, "both_minus")
	assert info["B.1"][4:6] == ["ctgB", "13"]


def test_query_sequences_come_from_the_alignment_header_or_are_not_reported(tmp_path, caplog):
	rows = [("B.1", 15, "C", "T")]
	ref, target = pair(11, {5: ("C", "T")})
	info, _ = query_case(tmp_path, (10, 20, "+", ref), (41, 51, "+", target), rows, "header", header_file=True, query_file=False)
	assert info["B.1"][4:6] == ["ctgB", "16"], "the file named by #Sequence2File is read"
	info, _ = query_case(tmp_path, (10, 20, "+", ref), (41, 51, "+", target), rows, "none", query_file=False)
	assert info["B.1"][4:6] == ["-", "N/A"], "without the query sequences no position is invented"
	assert info["B.1"][8] == "T", "the SNP itself is still called"
	assert "could not be read" in caplog.text


TWO_SEQUENCES = (("chr1", 600), ("chr2", 400))


def test_xmfa_maps_joined_coordinates_to_sequences(tmp_path):
	# progressiveMauve counts along the reference sequences joined end to end: with chr1 of
	# 600 bases, joined 595 is chr1:595 and joined 615 is chr2:15. One block 590-620 spans
	# the end of chr1 and the start of chr2 and holds both SNPs. The last base of chr1
	# (joined 600) and the first base of chr2 (joined 601) test the edge between them
	rows = [("B.1", "chr1", 595, "C", "T"), ("B.2", "chr2", 15, "G", "A"),
		("E.1", "chr1", 600, "C", "A"), ("E.2", "chr2", 1, "C", "G")]
	ref, target = pair(31, {595 - 590: ("C", "T"), 615 - 590: ("G", "A"), 600 - 590: ("C", "A"), 601 - 590: ("C", "G")})
	snps, called, parser = run(tmp_path, rows, [(590, 620, ref, target)], sequences=TWO_SEQUENCES)
	assert snps == {"B.1": 1, "B.2": 1, "E.1": 1, "E.2": 1}
	assert called == ["B.1", "E.1", "E.2", "B.2"]
	assert [(row[0], row[2], row[3], row[8]) for row in parser.get_snp_info()] == [
		("B.1", "chr1", "595", "T"), ("E.1", "chr1", "600", "A"), ("E.2", "chr2", "1", "G"), ("B.2", "chr2", "15", "A")]
	assert parser.joined == [595, 600, 601, 615]
	assert parser.joined_keys == {595: ("chr1", 595), 600: ("chr1", 600), 601: ("chr2", 1), 615: ("chr2", 15)}
	assert parser.joined_snps_in(590, 600) == [595, 600]
	assert parser.joined_snps_in(601, 620) == [601, 615]


def test_xmfa_reverse_block_on_second_sequence(tmp_path):
	# a "-" block over joined 610-630 is chr2:10-30 reverse complemented; chr2:12 is joined
	# 612 and sits in column 630 - 612 = 18, its query base is complemented before the call
	rows = [("B.1", "chr2", 12, "C", "T")]
	ref, target = pair(21, {630 - 612: ("G", "A")})
	snps, called, parser = run(tmp_path, rows, [(610, 630, ref, target, "-")], sequences=TWO_SEQUENCES)
	assert snps == {"B.1": 1}
	row = parser.get_snp_info()[0]
	assert (row[2], row[3], row[8]) == ("chr2", "12", "T")


def test_xmfa_without_sequence_table_fails_with_realign_hint(tmp_path):
	# alignments made before 2.3.0 have no table, and their coordinates cannot be mapped to
	# the reference sequences; the run must stop and say how to recover
	rows = [("B.1", 15, "C", "T")]
	ref, target = pair(11, {5: ("C", "T")})
	with pytest.raises(FileNotFoundError) as error:
		run(tmp_path, rows, [(10, 20, ref, target)], table=False)
	assert "align again" in str(error.value)


def test_xmfa_snps_outside_the_table_are_left_out(tmp_path, caplog):
	# a SNP on a sequence the aligned reference does not have, or beyond the end of its
	# sequence, has no joined coordinate; it would otherwise be read from the next sequence
	rows = [("B.1", "chr1", 595, "C", "T"), ("X.1", "plasmid", 5, "C", "T"), ("X.2", "chr1", 605, "C", "T")]
	ref, target = pair(31, {595 - 590: ("C", "T"), 605 - 590: ("C", "T")})
	snps, called, parser = run(tmp_path, rows, [(590, 620, ref, target)], sequences=TWO_SEQUENCES)
	assert snps == {"B.1": 1}
	assert parser.joined == [595]
	assert "plasmid" in caplog.text and "X.2" in caplog.text


def test_mauve_prepare_writes_the_table(tmp_path):
	reference = tmp_path / "REF.fna"
	reference.write_text(">chr1 first\n" + "A" * 60 + "\n" + "C" * 40 + "\n>chr2\n" + "G" * 30 + "\n")
	query = tmp_path / "q.fna"
	query.write_text(">q\nACGT\n")
	output = str(tmp_path / "REF_q.xmfa")
	ProgressiveMauve().prepare(str(reference), str(query), output)
	assert os.path.isfile(output + SEQUENCES_SUFFIX)
	assert read_sequence_table(output + SEQUENCES_SUFFIX) == [("chr1", 100), ("chr2", 30)]


@pytest.mark.skipif(shutil.which("progressiveMauve") is None, reason="progressiveMauve is not installed")
def test_progressive_mauve_through_the_whole_alignment_step(tmp_path):
	# a two-sequence reference with SNPs on both sequences; the alignment step writes the
	# sequence table and the reader reports the SNPs per reference sequence
	rng = random.Random(11)
	chr1 = "".join(rng.choice("ACGT") for _ in range(6000))
	chr2 = "".join(rng.choice("ACGT") for _ in range(4000))
	joined = chr1 + chr2
	snps = {1500: ("B.1", "chr1", 1500), 4500: ("B.2", "chr1", 4500), 6000 + 2000: ("B.3", "chr2", 2000)}
	rows, query = [], list(joined)
	for pos, (snp_id, sequence, seq_pos) in snps.items():
		ref_base = joined[pos - 1]
		alt = "A" if ref_base != "A" else "C"
		rows.append((snp_id, sequence, seq_pos, ref_base, alt))
		if snp_id != "B.2":                     ## B.2 stays ancestral in the query
			query[pos - 1] = alt
	query = "".join(query)
	refdir = tmp_path / "refs"
	refdir.mkdir()
	(refdir / "REF.fna").write_text(">chr1\n{0}\n>chr2\n{1}\n".format(chr1, chr2))
	qfile = tmp_path / "sample.fna"
	qfile.write_text(">q1\n{0}\n>q2\n{1}\n".format(query[:5000], query[5000:]))
	db = str(tmp_path / "t.db")
	make_db(db, rows)
	obj = CanSNPer2([str(qfile)], workdir=str(tmp_path), refdir="refs", tmpdir="shared", logdir="logs", outdir="out",
		database=db, aligner="progressiveMauve")
	files = obj.align(str(qfile))
	assert files and all(os.path.isfile(f + SEQUENCES_SUFFIX) for f in files)
	states, info, called = obj.find_snps_multiproc(xmfa_obj=obj.aligner.reader(database=db), xmfa_files=files,
		references=obj.alignment_references)
	assert states == {"B.1": 1, "B.2": 2, "B.3": 1}
	assert {row[0]: (row[2], row[3]) for row in info} == {"B.1": ("chr1", "1500"), "B.2": ("chr1", "4500"), "B.3": ("chr2", "2000")}
