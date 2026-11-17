"""Regression tests for issue #39: SNPs added to the database as custom SNPs
were not called because the scan of an alignment block could stop early.

Each test feeds ParseXMFA a small database (SNP rows inserted the way
ModifyCanSNPer2Database.add_annotation does it) and a hand-made xmfa, and
checks which SNPs end up evaluated and called."""
import sqlite3

import pytest

from CanSNPer2.modules.ParseXMFA import ParseXMFA


def make_db(path, rows):
	'''rows: (snp_id, position, ancestral_base, derived_base) against reference REF'''
	conn = sqlite3.connect(path)
	conn.executescript("""
		CREATE TABLE nodes (id integer PRIMARY KEY, name text NOT NULL);
		CREATE TABLE tree (parent integer, child integer, rank_i integer);
		CREATE TABLE snp_annotation (node_id integer PRIMARY KEY, snp_id varchar(6), position integer,
			ancestral_base varchar(1), derived_base varchar(1), reference varchar(20),
			date datetime, genome_i integer);
		CREATE TABLE snp_references (id integer PRIMARY KEY, genome varchar(6));
		INSERT INTO snp_references (id, genome) VALUES (5, 'REF');
	""")
	for i, (snp_id, pos, anc, der) in enumerate(rows):
		conn.execute("INSERT INTO nodes (id, name) VALUES (?, ?)", (i + 1, snp_id))
		conn.execute("INSERT INTO snp_annotation (node_id, snp_id, position, ancestral_base, derived_base, reference, date, genome_i)"
			" VALUES (?, ?, ?, ?, ?, 'REF', '06/11/26', 5)", (i + 1, snp_id, pos, anc, der))
	conn.commit()
	conn.close()


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


def run(tmp_path, rows, blocks, name="t"):
	db, xmfa = tmp_path / (name + ".db"), tmp_path / ("REF_" + name + ".xmfa")
	make_db(db, rows)
	make_xmfa(xmfa, blocks)
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
	flags = {row[0]: row[8] for row in parser.get_snp_info()}
	assert flags == {"B.1": "edge", "M.1": "", "C.1": "edge"}, flags
	assert parser.edge_snps == {"B.1", "C.1"}


def test_query_contig_and_position_are_reported(tmp_path):
	# issue #23: every called position must say where it sits in the submitted assembly,
	# not only in the reference; the query count skips target gaps, and a "-" block maps
	# back through the reversed block
	rows = [("B.1", 15, "C", "T")]
	ref = list("C" * 11); target = list("C" * 11)
	ref[5] = "C"; target[5] = "T"
	target[3] = "-"                                  ## gap before the SNP: query is one base earlier
	ref_gap, target_gap = "".join(ref), "".join(target)
	_, _, parser = run(tmp_path, rows, [(10, 20, ref_gap, target_gap)], name="gap")
	info = {row[0]: row for row in parser.get_snp_info()}
	assert info["B.1"][3] == "query", "the contig name from the xmfa block header"
	assert info["B.1"][4] == "14", "start 10 + 5 columns - 1 skipped gap"

	rows = [("B.1", 15, "C", "T")]
	ref2, target2 = pair(11, {5: ("C", "T")})
	_, _, parser = run(tmp_path, rows, [(10, 20, ref2, target2, "-")], name="rev")
	info = {row[0]: row for row in parser.get_snp_info()}
	assert info["B.1"][4] == "15", "end 20 - 5 bases from the far side of the reversed block"


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
	assert info["B.1"][7] == "T" and info["B.2"][7] == "A"
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
	assert [row[7] for row in parser.get_snp_info()] == ["N"]
	assert "more than one alignment" in caplog.text
	# the same base twice is one observation
	snps, called, parser = run(tmp_path, rows, [(10, 20, ref1, target1), (10, 20, ref1, target1)], name="same")
	assert snps == {"B.1": 1} and len(parser.get_snp_info()) == 1
