"""Choosing the local reference: the representatives of the terminal nodes below the call of the
first round (genomes table), screened with a subsample of the reads mapped to all of them at
once (LocalReferences.candidates, screen)."""
import random
import shutil

import pytest

import helpers
from CanSNPer2.modules.ImportDatabase import import_from_sources
from CanSNPer2.modules.LocalReferences import candidates, screen

TREE = "R\nR\tA.1\nR\tA.1\tA.2\nR\tA.1\tA.3\nR\tB.1\n"
SNP_ROWS = [("R", 1, "P", "REF1", "chrI", 10, "A", "G", "", "yes"), ("A.1", 1, "P", "REF1", "chrI", 20, "T", "C", "", "yes"),
	("A.2", 1, "P", "REF1", "chrI", 30, "G", "A", "", "yes"), ("A.3", 1, "P", "REF1", "chrI", 35, "G", "A", "", "yes"),
	("B.1", 1, "P", "REF1", "chrI", 40, "A", "C", "", "yes")]
GENOMES = ("genome\tnode\trole\n" "GA2\tA.2\trepresentative\n" "GA2m\tA.2\tmember\n" "GA3\tA.3\trepresentative\n"
	"GB1\tB.1\trepresentative\n" "GOUT\tR\toutgroup\n")


def scheme(tmp_path):
	files = helpers.write_sources(tmp_path, SNP_ROWS, tree=TREE, name="src")
	(tmp_path / "g.txt").write_text(GENOMES)
	db = str(tmp_path / "s.db")
	import_from_sources(files["tree"], files["snp"], files["references"], db, genomes=str(tmp_path / "g.txt"))
	local = tmp_path / "refs" / "local"
	local.mkdir(parents=True)
	return db, str(tmp_path / "refs")


def test_candidates_below_a_node(tmp_path):
	db, refdir = scheme(tmp_path)
	for g in ("GA2", "GA3", "GB1"):
		(tmp_path / "refs" / "local" / (g + ".fna")).write_text(">c\nACGT\n")
	found, missing = candidates(db, "A.1", refdir)
	assert [c.genome for c in found] == ["GA2", "GA3"] and missing == []
	assert {c.node for c in found} == {"A.2", "A.3"}
	found, _ = candidates(db, "NA", refdir)
	assert [c.genome for c in found] == ["GA2", "GA3", "GB1"], "no call: every representative, never members or the outgroup"
	found, _ = candidates(db, "A.2", refdir)
	assert [c.genome for c in found] == ["GA2"]


def test_candidates_without_their_file_are_listed_as_missing(tmp_path):
	db, refdir = scheme(tmp_path)
	(tmp_path / "refs" / "local" / "GA2.fna").write_text(">c\nACGT\n")
	found, missing = candidates(db, "A.1", refdir)
	assert [c.genome for c in found] == ["GA2"] and missing == ["GA3"]


def test_candidates_are_capped_and_spread_over_the_subtrees(tmp_path):
	db, refdir = scheme(tmp_path)
	for g in ("GA2", "GA3", "GB1"):
		(tmp_path / "refs" / "local" / (g + ".fna")).write_text(">c\nACGT\n")
	found, _ = candidates(db, "R", refdir, limit=2)
	assert {c.genome for c in found} in ({"GA2", "GB1"}, {"GA3", "GB1"}), "one from each child subtree of R first"


@pytest.mark.skipif(shutil.which("minimap2") is None, reason="minimap2 is not installed")
def test_screen_picks_the_genome_the_reads_come_from(tmp_path):
	rng = random.Random(4)
	base = [rng.choice("ACGT") for _ in range(30000)]
	genomes = {}
	for name in ("G1", "G2", "G3"):
		g = [b if rng.random() > 0.08 else rng.choice("ACGT".replace(b, "")) for b in base]
		genomes[name] = "".join(g)
		(tmp_path / (name + ".fna")).write_text(">chr\n" + genomes[name] + "\n")
	r1, r2 = helpers.simulate_reads(genomes["G2"], 20, seed=3)
	(tmp_path / "r1.fq").write_text(r1)
	(tmp_path / "r2.fq").write_text(r2)
	ranked = screen([str(tmp_path / "r1.fq"), str(tmp_path / "r2.fq")],
		{name: str(tmp_path / (name + ".fna")) for name in genomes}, preset="sr", n_reads=2000, tmpdir=str(tmp_path))
	assert ranked[0][0] == "G2"
	assert ranked[0][1] > 0.8, "the share of the matched bases that went to G2"
