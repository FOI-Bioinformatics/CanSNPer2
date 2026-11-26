"""db modify and db remove on a small scheme, run in-process through the cansnper2 parsers.

Each test changes a database and compares its export (the tree links, the canSNPs and the
reference genomes) with the expected scheme, so that a change that leaves the database in a
state the export or a rebuild cannot read is found as well (#102)."""

import os
import sqlite3

import pytest

import helpers
import CanSNPer2.SNPDatabase as database_command
from CanSNPer2.cli import build_parsers
from CanSNPer2.modules.SchemeInfo import read_scheme_info

#      R
#     / \
#   A.1  B.1
#   / \
# A.2  A.4
#  |
# A.3
TREE = "R\nR\tA.1\nR\tA.1\tA.2\nR\tA.1\tA.2\tA.3\nR\tA.1\tA.4\nR\tB.1\n"
POSITIONS = {"R": 10, "A.1": 100, "A.2": 200, "A.3": 300, "A.4": 400, "B.1": 500, "C.1": 600, "C.2": 700}
SNP_HEADER = "\t".join(helpers.CANSNP_HEADER)


def snp_row(node, genome="REF1", sequence="CP000001.1"):
    return "{n}\t1\tPaper2020\t{g}\t{s}\t{p}\tT\tC\t\tyes".format(
        n=node, g=genome, s=sequence, p=POSITIONS[node]
    )


REFERENCES = (
    "genome\tstrain\tgenbank_id\trefseq_id\tassembly_name\n"
    "REF1\ts1\tGCA_000000001.1\tGCF_000000001.1\tASM1v1\n"
)


def cansnper2(*argv):
    """run cansnper2 in this process and return its exit status"""
    parser, _parsers = build_parsers()
    args = parser.parse_args([str(a) for a in argv])
    try:
        database_command.run_command(args)
    except SystemExit as e:
        return e.code or 0
    return 0


@pytest.fixture(autouse=True)
def in_tmp_path(tmp_path, monkeypatch):
    """the log folder of a call without --log-dir is ./logs"""
    monkeypatch.chdir(tmp_path)


@pytest.fixture
def scheme(tmp_path):
    files = helpers.write_sources(
        tmp_path / "src",
        [snp_row(n) for n in ("R", "A.1", "A.2", "A.3", "A.4", "B.1")],
        tree=TREE,
        references=REFERENCES,
        header=SNP_HEADER,
    )
    db = tmp_path / "s.db"
    status = cansnper2(
        "db", "create", "-d", db, "--tree", files["tree"], "--cansnps", files["snp"],
        "--references", files["references"], "--log-dir", tmp_path / "logs",
    )  # fmt: skip
    assert status == 0
    return db


def write(path, text):
    path.write_text(text)
    return path


def exported(db, tmp_path, name="exp"):
    """the links, the canSNP rows by node and the reference rows of an export of db"""
    out = tmp_path / name
    assert cansnper2("db", "export", "-d", db, "-o", out, "--log-dir", tmp_path / "logs") == 0
    stem = os.path.basename(str(db)).rsplit(".", 1)[0]
    links = set()
    for line in (out / (stem + "_tree.txt")).read_text().splitlines():
        path = [p for p in line.split("\t") if p]
        if len(path) > 1:
            links.add((path[-2], path[-1]))
    rows = (out / (stem + "_canSNPs.txt")).read_text().splitlines()
    snps = {row.split("\t")[0]: row for row in rows[1:]}
    references = (out / (stem + "_references.txt")).read_text()
    return links, snps, references, out


def node_names(db):
    conn = sqlite3.connect(str(db))
    names = {name for (name,) in conn.execute("SELECT name FROM nodes")}
    conn.close()
    return names


ALL_LINKS = {("R", "A.1"), ("A.1", "A.2"), ("A.2", "A.3"), ("A.1", "A.4"), ("R", "B.1")}


def test_the_export_of_the_scheme_is_the_scheme(scheme, tmp_path):
    links, snps, references, _ = exported(scheme, tmp_path)
    assert links == ALL_LINKS
    assert sorted(snps) == ["A.1", "A.2", "A.3", "A.4", "B.1", "R"]


def test_modify_adds_a_branch_of_two_levels_with_its_cansnps(scheme, tmp_path):
    mod = write(tmp_path / "mod.txt", "parent\tchild\nA.3\tC.1\nC.1\tC.2\n")
    new = write(tmp_path / "new.txt", SNP_HEADER + "\n" + snp_row("C.1") + "\n" + snp_row("C.2") + "\n")
    status = cansnper2(
        "db", "modify", "-d", scheme, "--mod-file", mod, "--parent", "A.3", "--cansnps", new,
        "--log-dir", tmp_path / "logs",
    )  # fmt: skip
    assert status == 0
    links, snps, _, _ = exported(scheme, tmp_path)
    assert links == ALL_LINKS | {("A.3", "C.1"), ("C.1", "C.2")}
    assert snps["C.1"].split("\t")[5] == "600" and snps["C.2"].split("\t")[5] == "700"
    info = read_scheme_info(str(scheme))
    assert info["status"] == "unknown" and info.get("modified")


def test_a_level_column_is_stored_with_the_link(scheme, tmp_path):
    mod = write(tmp_path / "mod.txt", "parent\tchild\tlevel\nB.1\tC.1\tsubtype\n")
    assert cansnper2("db", "modify", "-d", scheme, "--mod-file", mod, "--parent", "B.1") == 0
    conn = sqlite3.connect(str(scheme))
    rank = conn.execute(
        "SELECT rank.rank FROM tree JOIN nodes ON tree.child = nodes.id JOIN rank ON tree.rank_i = rank.rank_i"
        " WHERE nodes.name = 'C.1'"
    ).fetchall()
    conn.close()
    assert rank == [("subtype",)]


def test_the_same_modification_twice_changes_nothing_the_second_time(scheme, tmp_path):
    mod = write(tmp_path / "mod.txt", "parent\tchild\nA.3\tC.1\n")
    assert cansnper2("db", "modify", "-d", scheme, "--mod-file", mod, "--parent", "A.3") == 0
    first = exported(scheme, tmp_path, "first")[:3]
    assert cansnper2("db", "modify", "-d", scheme, "--mod-file", mod, "--parent", "A.3") == 0
    assert exported(scheme, tmp_path, "second")[:3] == first


def test_replace_removes_the_nodes_the_file_does_not_list(scheme, tmp_path):
    """below A.1 the file keeps A.2 and adds C.1; A.3 and A.4 go, with their canSNPs"""
    mod = write(tmp_path / "mod.txt", "parent\tchild\nA.1\tA.2\nA.1\tC.1\n")
    new = write(tmp_path / "new.txt", SNP_HEADER + "\n" + snp_row("C.1") + "\n")
    status = cansnper2(
        "db", "modify", "-d", scheme, "--mod-file", mod, "--parent", "A.1", "--replace", "--cansnps", new,
    )  # fmt: skip
    assert status == 0
    links, snps, _, _ = exported(scheme, tmp_path)
    assert links == {("R", "A.1"), ("A.1", "A.2"), ("A.1", "C.1"), ("R", "B.1")}
    assert sorted(snps) == ["A.1", "A.2", "B.1", "C.1", "R"]
    assert node_names(scheme) == {"R", "A.1", "A.2", "B.1", "C.1"}, "no node is left without a link"


def test_remove_takes_an_inner_node_with_its_descendants(scheme, tmp_path):
    assert cansnper2("db", "remove", "-d", scheme, "--parent", "A.2") == 0
    links, snps, _, _ = exported(scheme, tmp_path)
    assert links == {("R", "A.1"), ("A.1", "A.4"), ("R", "B.1")}
    assert sorted(snps) == ["A.1", "A.4", "B.1", "R"]
    assert node_names(scheme) == {"R", "A.1", "A.4", "B.1"}


def test_remove_takes_a_leaf(scheme, tmp_path):
    assert cansnper2("db", "remove", "-d", scheme, "--parent", "B.1") == 0
    links, snps, _, _ = exported(scheme, tmp_path)
    assert links == ALL_LINKS - {("R", "B.1")}
    assert "B.1" not in snps


@pytest.mark.parametrize("node", ["R", "NO.SUCH.NODE"])
def test_remove_of_the_root_or_an_unknown_node_leaves_the_database_unchanged(scheme, tmp_path, node):
    before = scheme.read_bytes()
    assert cansnper2("db", "remove", "-d", scheme, "--parent", node) == 1
    assert scheme.read_bytes() == before


def test_a_changed_database_is_rebuilt_from_its_export(scheme, tmp_path):
    """modify, then remove; the export of the result builds a database with the same export"""
    mod = write(tmp_path / "mod.txt", "parent\tchild\nA.4\tC.1\nC.1\tC.2\n")
    new = write(tmp_path / "new.txt", SNP_HEADER + "\n" + snp_row("C.1") + "\n" + snp_row("C.2") + "\n")
    assert (
        cansnper2("db", "modify", "-d", scheme, "--mod-file", mod, "--parent", "A.4", "--cansnps", new) == 0
    )
    assert cansnper2("db", "remove", "-d", scheme, "--parent", "A.2") == 0
    links, snps, references, out = exported(scheme, tmp_path)
    assert links == {("R", "A.1"), ("A.1", "A.4"), ("A.4", "C.1"), ("C.1", "C.2"), ("R", "B.1")}
    rebuilt = tmp_path / "rebuilt.db"
    status = cansnper2(
        "db", "create", "-d", rebuilt, "--tree", out / "s_tree.txt", "--cansnps", out / "s_canSNPs.txt",
        "--references", out / "s_references.txt",
    )  # fmt: skip
    assert status == 0
    again = exported(rebuilt, tmp_path, "again")
    assert again[0] == links and again[2] == references
    assert {k: v for k, v in again[1].items()} == snps
