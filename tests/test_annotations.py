"""Tree annotations (--annotations FILE).

Parts of the tree carry meaning the SNP calling itself does not know: a node may define a
type (A/M.1 = Type A, B.1 = Type B), a taxonomic rank (B.1, A.1 and B.16 can all be
subspecies while standing at different depths, because the database tree has uneven
resolution), and a branch may separate a feature such as an antibiotic resistance marker.
Annotations are read from a TSV (<snp name>\t<node|clade>\t<rank>\t<label>[\t<color>]),
travel in the embedded payload of both html views, and are never derived from tree depth."""

import json
import re
import sqlite3

import pytest

import helpers
from CanSNPer2.modules.NewickTree import NewickTree, load_annotations


def make_db(path):
    conn = sqlite3.connect(path)
    conn.executescript(
        helpers.SNP_ANNOTATION_TABLE
        + helpers.REFERENCES_TABLE
        + """
		CREATE TABLE nodes (id integer PRIMARY KEY, name text NOT NULL);
		CREATE TABLE tree (parent integer, child integer, rank_i integer, unique (parent, child, rank_i));
		INSERT INTO snp_references (id, genome, strain) VALUES (5, 'REF', 'ref strain');
		INSERT INTO nodes (id, name) VALUES (1, 'root'), (2, 'A.1'), (3, 'B.1'), (4, 'B.2'), (5, 'C.1');
		INSERT INTO tree (parent, child, rank_i) VALUES (1, 1, 1), (1, 2, 1), (2, 3, 1), (2, 4, 1), (4, 5, 1);
	"""
    )
    for node_id, snp_id, pos in ((2, "A.1", 10), (3, "B.1", 20), (4, "B.2", 30), (5, "C.1", 40)):
        conn.execute(
            "INSERT INTO snp_annotation (node_id, snp_id, sequence, sequence_position, ancestral_base, derived_base,"
            " reference, date, genome_i) VALUES (?, ?, 'chrI', ?, 'C', 'T', 'REF', '06/11/26', 5)",
            (node_id, snp_id, pos),
        )
    conn.commit()
    conn.close()


def tree(tmp_path):
    db = str(tmp_path / "t.db")
    make_db(db)
    return NewickTree(db, name="t", outdir=str(tmp_path))


def payload(path):
    page = open(path).read()
    match = re.search(r"var DATA = (\{.*\});", page)
    assert match, "the html must embed its json payload"
    return json.loads(match.group(1))


def test_load_annotations_parsing(tmp_path):
    path = str(tmp_path / "ann.tsv")
    open(path, "w").write(
        "# a comment\n"
        "snp\tkind\trank\tlabel\n"
        "B.1\tnode\ttype\tType B\t#d35400\n"
        "A.1\tnode\tsubspecies\tsubspecies A\n"
        "B.1\tclade\tmarker\tAMR feature\n"
        "\n"
    )
    ann = load_annotations(path)
    assert ann == [
        {"node": "B.1", "kind": "node", "rank": "type", "label": "Type B", "color": "#d35400"},
        {"node": "A.1", "kind": "node", "rank": "subspecies", "label": "subspecies A"},
        {"node": "B.1", "kind": "clade", "rank": "marker", "label": "AMR feature"},
    ]


def test_load_annotations_rejects_bad_rows(tmp_path):
    bad_kind = str(tmp_path / "bad1.tsv")
    open(bad_kind, "w").write("B.1\tbox\trank\tlabel\n")
    with pytest.raises(ValueError, match="kind must be node or clade"):
        load_annotations(bad_kind)
    too_few = str(tmp_path / "bad2.tsv")
    open(too_few, "w").write("B.1\tnode\tonly-three\n")
    with pytest.raises(ValueError, match="need <snp name>"):
        load_annotations(too_few)


def test_tree_viewer_carries_annotations(tmp_path):
    nt = tree(tmp_path)
    nt.draw_tree(
        {"A.1": 1},
        ["A.1"],
        save_tree=True,
        annotations=[
            {"node": "A.1", "kind": "node", "rank": "type", "label": "Type A"},
            {"node": "B.1", "kind": "clade", "rank": "marker", "label": "AMR feature"},
            {"node": "NOWHERE.9", "kind": "node", "rank": "type", "label": "wrong database"},
        ],
    )
    data = payload(nt.tree_file)
    assert [a["node"] for a in data["annotations"]] == ["A.1", "B.1"], "unknown names are dropped"
    assert "window.onerror" in open(nt.tree_file).read(), "the viewer must report its own script failures"
    assert "did not draw" in open(nt.tree_file).read(), "the viewer must notice an empty tree"


def test_dashboard_carries_annotations(tmp_path):
    nt = tree(tmp_path)
    out = nt.write_dashboard(
        [
            {
                "name": "s",
                "severity": "clean",
                "final": "A.1",
                "path": "A.1",
                "reasons": [],
                "states": {"A.1": 1},
                "called": ["A.1"],
                "edge": [],
                "no_data": [],
            }
        ],
        str(tmp_path / "summary.html"),
        [{"node": "B.1", "kind": "node", "rank": "subspecies", "label": "subspecies B"}],
    )
    assert payload(out)["annotations"][0]["rank"] == "subspecies"


def test_annotations_are_not_derived_from_depth(tmp_path):
    """the same rank may sit at different depths: A.1 is one level under the root, C.1 three
    levels below through B.1, yet a file can give them the same rank"""
    path = str(tmp_path / "ann.tsv")
    open(path, "w").write("A.1\tnode\tsubspecies\tsub A\nC.1\tnode\tsubspecies\tsub C\n")
    ann = load_annotations(path)
    nt = tree(tmp_path)
    entries = nt._tree_entries()
    depth = {e["name"]: e["depth"] for e in entries}
    assert depth["A.1"] != depth["C.1"], "the fixture tree must have the two at different depths"
    assert ann[0]["rank"] == ann[1]["rank"], "the file alone decides the rank"
