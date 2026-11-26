"""Low-severity defects of building and modifying CanSNPer2 databases.

Positions above the range of an integer, empty or header-only source files, unsafe genome
names, repeated column names or accessions, a tree path of thousands of levels, Newick
names with special characters, trailing tabs, and the node ids of a modification."""

import re
import sqlite3
import sys

import pytest

import helpers
from CanSNPer2.modules.DatabaseConnection import SourceFileError
from CanSNPer2.modules.NewickTree import NewickTree
from CanSNPer2.modules.SourceCheck import check_sources
from test_create_input_checks import REFERENCES, SNP_HEADER, SNP_ROWS, TREE, create, refused, write_sources
from test_db_build_and_inputs import build_db, run_database_cli, tree_file, annotation_file

pytestmark = pytest.mark.skipif(sys.platform == "win32", reason="posix paths")


def problems(tmp_path, tree=TREE, snp=SNP_HEADER + SNP_ROWS, references=REFERENCES):
    files = write_sources(tmp_path, tree=tree, snp=snp, references=references)
    return check_sources(tree=files["tree"], annotation=files["snp"], references=files["references"])


def joined(found):
    return " | ".join(found)


def test_position_above_the_integer_range_is_refused(tmp_path):
    snp = SNP_HEADER + SNP_ROWS.replace("\t100\t", "\t99999999999999999999999\t")
    assert "2147483647" in joined(problems(tmp_path, snp=snp))
    assert problems(tmp_path, snp=SNP_HEADER + SNP_ROWS.replace("\t100\t", "\t2147483647\t")) == []


def test_replace_with_an_empty_modification_file_is_refused(tmp_path):
    db = build_db(tmp_path, "m.db", tree_file("\t"), annotation_file())
    mod = tmp_path / "mod.txt"
    mod.write_text("parent\tchild\n")
    before = open(db, "rb").read()
    result = run_database_cli(tmp_path, ["-db", db, "--mod_file", str(mod), "--parent", "R", "--replace"])
    assert result.returncode == 1
    assert "no rows" in result.stderr + result.stdout
    assert open(db, "rb").read() == before


def test_new_nodes_of_a_modification_continue_after_the_highest_id(tmp_path):
    db = build_db(tmp_path, "m.db", tree_file("\t"), annotation_file())
    mod = tmp_path / "mod.txt"
    mod.write_text("parent\tchild\nA.1\tC.1\nC.1\tC.2\n")
    top = sqlite3.connect(db).execute("SELECT MAX(id) FROM nodes").fetchone()[0]
    result = run_database_cli(tmp_path, ["-db", db, "--mod_file", str(mod), "--parent", "A.1"])
    assert result.returncode == 0, result.stderr
    ids = [
        row[0]
        for row in sqlite3.connect(db).execute(
            "SELECT id FROM nodes WHERE name IN ('C.1', 'C.2') ORDER BY id"
        )
    ]
    assert ids == [top + 1, top + 2]


def test_empty_snp_file_is_said_to_be_empty(tmp_path):
    found = joined(problems(tmp_path, snp=""))
    assert "empty" in found and "no column" not in found


def test_empty_references_file_is_said_to_be_empty(tmp_path):
    found = joined(problems(tmp_path, references=""))
    assert "empty" in found and "no column" not in found


def test_snp_file_with_only_a_header_is_refused(tmp_path):
    files = write_sources(tmp_path, snp=SNP_HEADER)
    message = refused(tmp_path, files)
    assert "no canSNP rows" in message


@pytest.mark.parametrize("name", ["a/b", "../x", ".hidden", "a b", "a\\b", "a:b"])
def test_unsafe_genome_names_are_refused(tmp_path, name):
    found = problems(
        tmp_path, snp=SNP_HEADER + SNP_ROWS.replace("REF1", name), references=REFERENCES.replace("REF1", name)
    )
    assert "genome name" in joined(found) and name in joined(found)


def test_duplicate_column_names_are_refused(tmp_path):
    header = SNP_HEADER.rstrip("\n") + "\tscore\n"
    rows = "".join(line + "\tx\n" for line in SNP_ROWS.splitlines())
    assert "'score'" in joined(problems(tmp_path, snp=header + rows))
    references = REFERENCES.replace("assembly_name", "strain")
    assert "'strain'" in joined(problems(tmp_path, references=references))


def test_two_genomes_with_one_refseq_id_are_refused(tmp_path):
    references = REFERENCES.replace("GCF_000000002.1", "GCF_000000001.1")
    found = joined(problems(tmp_path, references=references))
    assert "GCF_000000001.1" in found and "refseq_id" in found


def test_a_tree_path_of_3000_levels_can_be_created(tmp_path):
    names = ["N{0}".format(i) for i in range(3000)]
    tree = "\t".join(names) + "\n"
    snp_rows = ["{0}\t1\tP\tREF1\tCP000001.1\t{1}\tA\tG\t\tyes".format(names[-1], 10)]
    files = helpers.write_sources(
        tmp_path,
        snp_rows,
        tree=tree,
        references="genome\tstrain\tgenbank_id\trefseq_id\tassembly_name\nREF1\ts1\tGCA_000000001.1\tGCF_000000001.1\tasm\n",
    )
    db = create(tmp_path, files)
    assert sqlite3.connect(db).execute("SELECT COUNT(*) FROM nodes").fetchone()[0] == 3000
    newick = NewickTree(db).newickTree
    assert newick.count("(") == 3000 and newick.endswith("ROOT;")


def test_newick_of_a_tree_with_descending_ids_and_3000_levels(tmp_path):
    db = str(tmp_path / "deep.db")
    conn = sqlite3.connect(db)
    conn.executescript(
        """CREATE TABLE nodes (id integer PRIMARY KEY, name text NOT NULL);
		CREATE TABLE tree (parent integer, child integer, rank_i integer);
		"""
        + helpers.SNP_ANNOTATION_TABLE
    )
    ## the root has id 1 and the ids fall with the depth: a node is read before its parent
    n = 3000
    conn.execute("INSERT INTO nodes VALUES (1, 'N0')")
    conn.execute("INSERT INTO tree VALUES (1, 1, 1)")
    for level in range(1, n):
        node = n + 1 - level
        parent = 1 if level == 1 else node + 1
        conn.execute("INSERT INTO nodes VALUES (?, ?)", (node, "N{0}".format(level)))
        conn.execute("INSERT INTO tree VALUES (?, ?, 1)", (parent, node))
    conn.commit()
    conn.close()
    assert NewickTree(db).newickTree.endswith("ROOT;")


def newick_of(tmp_path, names):
    db = str(tmp_path / "q{0}.db".format(len(list(tmp_path.glob("q*.db")))))
    conn = sqlite3.connect(db)
    conn.executescript(
        """CREATE TABLE nodes (id integer PRIMARY KEY, name text NOT NULL);
		CREATE TABLE tree (parent integer, child integer, rank_i integer);
		"""
        + helpers.SNP_ANNOTATION_TABLE
    )
    for i, name in enumerate(names, start=1):
        conn.execute("INSERT INTO nodes VALUES (?, ?)", (i, name))
        conn.execute("INSERT INTO tree VALUES (?, ?, 1)", (1, i))
    conn.commit()
    conn.close()
    return NewickTree(db).newickTree


def test_newick_names_with_special_characters_are_quoted(tmp_path):
    tree = newick_of(tmp_path, ["R", "a(b)", "c,d", "e f", "g'h", "i:j;k"])
    ## a node without a SNP is written as -name-
    for quoted in ("'-a(b)-'", "'-c,d-'", "'-e f-'", "'-g''h-'", "'-i:j;k-'"):
        assert quoted in tree
    assert newick_of(tmp_path, ["R", "A.1", "B-2"]) == "((-A.1-,-B-2-)-R-)ROOT;"


def test_trailing_tab_on_header_and_rows_is_handled_alike(tmp_path):
    snp = "\n".join(line + "\t" for line in (SNP_HEADER + SNP_ROWS).splitlines()) + "\n"
    references = "\n".join(line + "\t" for line in REFERENCES.splitlines()) + "\n"
    files = write_sources(tmp_path, snp=snp, references=references)
    assert check_sources(tree=files["tree"], annotation=files["snp"], references=files["references"]) == []
    db = create(tmp_path, files)
    assert sqlite3.connect(db).execute("SELECT COUNT(*) FROM snp_annotation").fetchone()[0] == 4
