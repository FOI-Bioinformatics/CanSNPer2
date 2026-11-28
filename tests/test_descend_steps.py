"""The steps of the walk that chooses the local reference of a read sample (LocalReferences.descend
split into Descent, #103), each tested on its own: the tree and the representative genomes read
from the database, the node the walk starts at, one step down by the counts or by the screen, and
the genome it ends on. The walk as a whole is tested in test_local_selection.py."""

from CanSNPer2.modules.LocalReferences import Descent, tree_and_representatives
from test_local_selection import StubScreen, rows, scheme, with_files


def walk(db, refdir, states=None, screen=None, **kwargs):
    return Descent(db, refdir, rows(**(states or {})), screen or StubScreen(), **kwargs)


def test_tree_and_representatives_reads_the_tree_and_the_genomes_with_a_file(tmp_path):
    db, refdir = with_files(tmp_path, "GA2", "GB1")
    names, children, reps = tree_and_representatives(db, refdir)
    ids = {name: i for i, name in names.items()}
    assert set(names.values()) >= {"R", "A.1", "A.2", "A.3", "B.1"}
    assert [names[c] for c in children[ids["A.1"]]] == ["A.2", "A.3"]
    assert {c.genome for node in reps.values() for c in node} == {"GA2", "GB1"}, (
        "GA3 has no file and is not a candidate, nor are the member and the outgroup genome"
    )
    assert [c.genome for c in reps[ids["A.2"]]] == ["GA2"]


def test_tree_and_representatives_counts_a_genome_without_a_file_when_it_can_be_fetched(tmp_path):
    db, refdir = with_files(tmp_path, "GA2")
    _names, _children, reps = tree_and_representatives(db, refdir, can_fetch=True)
    assert {c.genome for node in reps.values() for c in node} == {"GA2", "GA3", "GB1"}


def test_the_walk_starts_at_the_call_or_at_the_root(tmp_path):
    db, refdir = with_files(tmp_path, "GA2", "GA3", "GB1")
    descent = walk(db, refdir)
    assert descent.names[descent.start_node("A.1")] == "A.1"
    assert descent.names[descent.start_node("NA")] == "R", "no call: the root of the tree"


def test_the_walk_has_no_start_without_a_tree(tmp_path):
    db, refdir = with_files(tmp_path / "empty")
    descent = walk(db, refdir)
    descent.children, descent.names = {}, {}
    assert descent.start_node("NA") is None


def test_the_subtree_of_a_node_holds_its_genomes(tmp_path):
    db, refdir = with_files(tmp_path, "GA2", "GA3", "GB1")
    descent = walk(db, refdir)
    a1 = descent.ids["A.1"]
    assert [descent.names[i] for i in descent.subtree(a1)] == ["A.1", "A.2", "A.3"]
    assert [c.genome for c in descent.genomes_in(a1)] == ["GA2", "GA3"]


def test_a_step_follows_the_counts_without_a_screen(tmp_path):
    db, refdir = with_files(tmp_path, "GA2", "GA3", "GB1")
    screen = StubScreen()
    descent = walk(db, refdir, {"A.1": "derived"}, screen)
    steps = []
    node, given_up = descent.step(descent.ids["R"], steps)
    assert descent.names[node] == "A.1" and not given_up
    assert screen.calls == []
    assert steps == [
        {"node": "A.1", "how": "counts", "share": None, "runner_up": None, "runner_up_node": None}
    ]


def test_a_step_screens_the_subtrees_the_counts_do_not_separate(tmp_path):
    db, refdir = with_files(tmp_path, "GA2", "GA3", "GB1")
    screen = StubScreen("GB1")
    descent = walk(db, refdir, {}, screen)
    steps = []
    node, given_up = descent.step(descent.ids["R"], steps)
    assert descent.names[node] == "B.1" and not given_up
    assert screen.calls == [["GA2", "GA3", "GB1"]]
    assert steps[0]["how"] == "screen" and steps[0]["share"] == 0.7
    assert steps[0]["runner_up_node"] == "A.1", "the best genome of the subtree it did not take"


def test_a_step_at_a_node_without_subtrees_ends_the_walk(tmp_path):
    db, refdir = with_files(tmp_path, "GA2")
    descent = walk(db, refdir)
    steps = []
    assert descent.step(descent.ids["A.2"], steps) == (None, False)
    assert steps == []


def test_a_first_screen_that_matches_nothing_gives_up_the_walk(tmp_path):
    db, refdir = with_files(tmp_path, "GA2", "GA3", "GB1")
    nothing = lambda genomes: [(g, 0.0, 0) for g in sorted(genomes)]
    descent = walk(db, refdir, {}, nothing)
    assert descent.step(descent.ids["R"], []) == (None, True)


def test_a_later_screen_that_matches_nothing_ends_the_walk_where_it_stands(tmp_path):
    db, refdir = with_files(tmp_path, "GA2", "GA3", "GB1")
    nothing = lambda genomes: [(g, 0.0, 0) for g in sorted(genomes)]
    descent = walk(db, refdir, {}, nothing)
    steps = [{"node": "A.1", "how": "counts"}]
    assert descent.step(descent.ids["A.1"], steps) == (None, False)


def test_the_walk_ends_on_the_representative_of_its_node(tmp_path):
    db, refdir = with_files(tmp_path, "GA2", "GA3", "GB1")
    descent = walk(db, refdir)
    assert descent.choose(descent.ids["A.2"]).genome == "GA2"
    assert descent.choose(descent.ids["A.1"]).genome in ("GA2", "GA3"), (
        "a node without a representative of its own: one from its subtree"
    )


def test_the_walk_ends_without_a_genome_when_none_can_be_fetched(tmp_path):
    db, refdir = scheme(tmp_path)
    descent = walk(db, refdir, fetch=lambda genomes: [])
    assert descent.choose(descent.ids["A.2"]) is None
