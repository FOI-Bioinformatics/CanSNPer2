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
	"GB1\tB.1\trepresentative\n" "GOUT\t\toutgroup\n")


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


# the descent: from the call down the tree, one representative per child subtree screened only
# where the counts of the first round cannot tell the subtrees apart

from CanSNPer2.modules.LocalReferences import descend


def rows(**states):
	return [{"Name": name, "state": state} for name, state in states.items()]


class StubScreen:
	'''a screen that answers with fixed winners and records the genomes it was given'''
	def __init__(self, *winners):
		self.winners = list(winners)
		self.calls = []

	def __call__(self, genomes):
		self.calls.append(sorted(genomes))
		win = self.winners.pop(0)
		return [(win, 0.7, 700)] + [(g, 0.3 / max(1, len(genomes) - 1), 100) for g in sorted(genomes) if g != win]


def with_files(tmp_path, *names):
	db, refdir = scheme(tmp_path)
	for g in names:
		(tmp_path / "refs" / "local" / (g + ".fna")).write_text(">c\nACGT\n")
	return db, refdir


def test_counts_that_point_to_one_subtree_need_no_screen(tmp_path):
	db, refdir = with_files(tmp_path, "GA2", "GA3", "GB1")
	stub = StubScreen()
	chosen, steps = descend(db, "NA", refdir, rows(**{"A.1": "derived", "A.2": "derived", "B.1": "no_reads"}), stub)
	assert chosen.genome == "GA2" and stub.calls == []
	assert [s["how"] for s in steps] == ["counts", "counts"]


def test_without_counts_the_representatives_of_each_subtree_are_screened(tmp_path):
	db, refdir = with_files(tmp_path, "GA2", "GA3", "GB1")
	stub = StubScreen("GB1")
	chosen, steps = descend(db, "NA", refdir, rows(), stub)
	assert chosen.genome == "GB1"
	assert stub.calls == [["GA2", "GA3", "GB1"]], "the representatives of both leaves of A, and of B"
	assert steps[0]["how"] == "screen" and steps[0]["share"] == 0.7


def test_the_descent_screens_again_at_each_branching(tmp_path):
	db, refdir = with_files(tmp_path, "GA2", "GA3", "GB1")
	stub = StubScreen("GA2", "GA3")
	chosen, steps = descend(db, "NA", refdir, rows(), stub)
	assert chosen.genome == "GA3"
	assert stub.calls[1] == ["GA2", "GA3"], "below A.1: the representatives of A.2 and A.3"
	assert [s["node"] for s in steps] == ["A.1", "A.3"]


def test_a_subtree_whose_first_canSNP_is_ancestral_is_left_out(tmp_path):
	db, refdir = with_files(tmp_path, "GA2", "GA3", "GB1")
	stub = StubScreen()
	chosen, _ = descend(db, "NA", refdir, rows(**{"A.1": "ancestral"}), stub)
	assert chosen.genome == "GB1" and stub.calls == []


def test_the_descent_starts_at_the_call(tmp_path):
	db, refdir = with_files(tmp_path, "GA2", "GA3", "GB1")
	stub = StubScreen("GA3")
	chosen, _ = descend(db, "A.1", refdir, rows(), stub)
	assert chosen.genome == "GA3" and stub.calls == [["GA2", "GA3"]]
	chosen, _ = descend(db, "A.2", refdir, rows(), StubScreen())
	assert chosen.genome == "GA2", "a terminal call: its own representative"


def test_subtrees_without_genome_files_are_passed_over(tmp_path):
	db, refdir = with_files(tmp_path, "GB1")
	chosen, _ = descend(db, "NA", refdir, rows(), StubScreen())
	assert chosen.genome == "GB1"
	db2, refdir2 = with_files(tmp_path / "none")
	assert descend(db2, "NA", refdir2, rows(), StubScreen()) == (None, [])


def test_reads_that_map_to_none_of_the_genomes_stop_the_descent(tmp_path):
	db, refdir = with_files(tmp_path, "GA2", "GA3", "GB1")
	chosen, steps = descend(db, "NA", refdir, rows(), lambda genomes: [(g, 0.0, 0) for g in sorted(genomes)])
	assert chosen is None


# when a local round can still change something

from CanSNPer2.modules.LocalReferences import open_canSNPs


def test_a_terminal_call_with_data_on_its_whole_path_needs_no_local_round(tmp_path):
	db, _ = scheme(tmp_path)
	assert open_canSNPs(db, "A.2", rows(**{"R": "derived", "A.1": "derived", "A.2": "derived"})) == []


def test_a_child_of_the_call_without_data_asks_for_a_local_round(tmp_path):
	db, _ = scheme(tmp_path)
	assert open_canSNPs(db, "A.1", rows(**{"R": "derived", "A.1": "derived", "A.2": "no_reads", "A.3": "ancestral"})) == ["A.2"]
	assert open_canSNPs(db, "A.1", rows(**{"R": "derived", "A.1": "derived", "A.2": "ancestral", "A.3": "mixed"})) == [], \
		"children with data (also a mixed one): the walk stopped for a reason another genome does not change"


def test_a_gap_on_the_path_asks_for_a_local_round(tmp_path):
	db, _ = scheme(tmp_path)
	assert open_canSNPs(db, "A.2", rows(**{"R": "low_depth", "A.1": "derived", "A.2": "derived"})) == ["R"]


def test_no_call_asks_for_a_local_round(tmp_path):
	db, _ = scheme(tmp_path)
	assert open_canSNPs(db, "NA", rows()) != []


def test_reads_that_hardly_map_to_the_best_genome_stop_the_descent(tmp_path):
	db, refdir = with_files(tmp_path, "GA2", "GA3", "GB1")
	def screen_fn(genomes):
		return [(g, 0.5, 100, 0.02) for g in sorted(genomes)]
	assert descend(db, "NA", refdir, rows(), screen_fn) == (None, [])


def test_up_to_three_representatives_of_each_subtree_are_screened(tmp_path):
	db, refdir = with_files(tmp_path, "GA2", "GA3", "GB1")
	stub = StubScreen("GA3", "GA3")
	chosen, steps = descend(db, "NA", refdir, rows(), stub)
	assert stub.calls[0] == ["GA2", "GA3", "GB1"], "both leaves of the A subtree stand for it"
	assert chosen.genome == "GA3" and [s["node"] for s in steps] == ["A.1", "A.3"]
