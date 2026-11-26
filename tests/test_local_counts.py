"""Round 2 of read typing with a local reference: the pileup on the local genome is read at the
lifted canSNP positions (counts of "-" sites complemented to the scheme reference strand) and
merged with the counts of round 1 (ParsePileup with lift, LocalReferences.merge_counts)."""

import pytest

from helpers import make_db
from CanSNPer2.modules.AlignmentReader import DERIVED, ANCESTRAL, NO_DATA
from CanSNPer2.modules.LiftOver import Lifted
from CanSNPer2.modules.LocalReferences import merge_counts, states_from_rows
from CanSNPer2.modules.ParsePileup import PileupReader, read_site_counts

ROWS = [("B.1", "chr", 100, "C", "G"), ("B.2", "chr", 200, "A", "T"), ("B.3", "chr", 300, "G", "A")]
LIFT = {"B.1": Lifted("c1", 50, "+", "G"), "B.2": Lifted("c2", 7, "-", "A")}  # B.3 not lifted


def read_lifted(tmp_path, text):
    db, path = tmp_path / "a.db", tmp_path / "16M_s.local.pileup"
    make_db(db, ROWS)
    path.write_text(text)
    reader = PileupReader(min_depth=5, min_allele_fraction=0.8, lift=LIFT, source="LOCAL1")
    snps = reader.run(str(path), reference="REF", database=str(db))
    return snps, read_site_counts(str(path) + ".counts")


def test_lifted_sites_are_read_and_minus_strand_counts_complemented(tmp_path):
    # c1:50 carries G (derived of B.1); c2:7 carries T on the local strand, A on the reference strand
    text = "c1\t50\tN\t6\tGGGggg\tIIIIII\nc2\t7\tN\t6\tTTTttt\tIIIIII\nc1\t51\tN\t6\tAAAAAA\tIIIIII\n"
    snps, rows = read_lifted(tmp_path, text)
    assert snps == {"B.1": DERIVED, "B.2": ANCESTRAL}
    by = {r["Name"]: r for r in rows}
    assert by["B.2"]["A"] == "6" and by["B.2"]["T"] == "0"
    assert by["B.2"]["lifted_position"] == "c2:7-" and by["B.1"]["lifted_position"] == "c1:50+"
    assert by["B.3"]["state"] == "not_lifted"
    assert {r["source"] for r in rows} == {"LOCAL1"}
    assert by["B.1"]["Sequence"] == "chr" and by["B.1"]["Pos"] == "100", "rows keep the scheme coordinates"


def row(name, state, depth, A=0, C=0, G=0, T=0, source="16M"):
    return {
        "Name": name,
        "Reference": "REF",
        "Sequence": "chr",
        "Pos": "1",
        "Ancestral base": "C",
        "Derived base": "G",
        "A": str(A),
        "C": str(C),
        "G": str(G),
        "T": str(T),
        "depth": str(depth),
        "fraction_derived": "0",
        "state": state,
        "source": source,
        "lifted_position": "",
    }


@pytest.mark.parametrize(
    "r1, r2, state, source",
    [
        (
            row("X", "no_reads", 0),
            row("X", "derived", 20, G=20, source="L"),
            "derived",
            "L",
        ),  # round 2 fills a gap
        (
            row("X", "derived", 30, G=30),
            row("X", "not_lifted", 0, source="L"),
            "derived",
            "16M",
        ),  # round 1 kept
        (
            row("X", "derived", 30, G=30),
            row("X", "low_depth", 2, G=2, source="L"),
            "derived",
            "16M",
        ),  # round 2 too thin
        (row("X", "low_depth", 3, G=3), row("X", "derived", 12, G=12, source="L"), "derived", "L"),
        (row("X", "ancestral", 30, C=30), row("X", "derived", 25, G=25, source="L"), "conflict", "16M,L"),
        (row("X", "mixed", 20, C=10, G=10), row("X", "derived", 25, G=25, source="L"), "derived", "L"),
    ],
)
def test_merge_counts(r1, r2, state, source):
    merged = merge_counts([r1], [r2], min_depth=5)
    assert [(m["state"], m["source"]) for m in merged] == [(state, source)]


def test_states_from_rows_give_the_tree_states_and_result_rows():
    rows = [
        row("A", "derived", 10, G=10),
        row("B", "ancestral", 8, C=8),
        row("C", "conflict", 9, C=5, G=4),
        row("D", "no_reads", 0),
        row("E", "mixed", 10, C=5, G=5),
    ]
    states, info, called = states_from_rows(rows)
    assert states == {"A": DERIVED, "B": ANCESTRAL, "C": NO_DATA, "E": NO_DATA}
    assert called == ["A"]
    assert {r[0]: r[8] for r in info} == {"A": "G", "B": "C", "C": "N", "E": "N"}
    assert info[0][4:6] == ["reads", "10"]
