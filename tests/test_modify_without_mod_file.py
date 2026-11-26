"""--annotation and/or --references change an existing database without --mod_file.

The README documents --annotation and --references as ways to add SNPs and reference
genomes to a database, but ModifyTree demanded a --mod_file and refused the run with a
message about an empty modification source."""

import sqlite3

import pytest

import helpers
from test_db_build_and_inputs import annotation_file, build_db, run_database_cli, tree_file
from CanSNPer2.modules.SchemeInfo import read_scheme_info

SNP_LINE = "\t".join(helpers.SNP_HEADER) + "\n"
NEW_SNP = "B.1\tfrancisella\tBirdsell2014\tREF1\tCP000001.1\t300\tT\tC\n"
NEW_REFERENCE = (
    "genome\tstrain\tgenbank_id\trefseq_id\tassembly_name\n"
    "REF2\ttwo\tGCA_000000002.1\tGCF_000000002.1\tasm2\n"
)


def database(tmp_path):
    return build_db(tmp_path, "m.db", tree_file("\t"), annotation_file())


def count(db, table, where="1"):
    conn = sqlite3.connect(db)
    try:
        return conn.execute("SELECT COUNT(*) FROM {0} WHERE {1}".format(table, where)).fetchone()[0]
    finally:
        conn.close()


def test_annotation_alone_adds_a_snp(tmp_path):
    db = database(tmp_path)
    snp = tmp_path / "new_snp.txt"
    snp.write_text(SNP_LINE + NEW_SNP)
    result = run_database_cli(tmp_path, ["-db", db, "--annotation", str(snp)])
    assert result.returncode == 0, result.stderr
    assert count(db, "snp_annotation", "snp_id = 'B.1' AND sequence_position = 300") == 1
    assert count(db, "nodes") == 3, "the tree is not changed"
    info = read_scheme_info(db)
    assert info["status"] == "unknown" and "modified" in info["commit"]


def test_references_alone_adds_a_genome(tmp_path):
    db = database(tmp_path)
    references = tmp_path / "new_references.txt"
    references.write_text(NEW_REFERENCE)
    result = run_database_cli(tmp_path, ["-db", db, "--references", str(references)])
    assert result.returncode == 0, result.stderr
    assert count(db, "snp_references", "genome = 'REF2'") == 1
    assert "modified" in read_scheme_info(db)["commit"]


def test_references_and_annotation_without_mod_file(tmp_path):
    db = database(tmp_path)
    references = tmp_path / "new_references.txt"
    references.write_text(NEW_REFERENCE)
    snp = tmp_path / "new_snp.txt"
    snp.write_text(SNP_LINE + "B.1\tfrancisella\tBirdsell2014\tREF2\tCP000002.1\t300\tT\tC\n")
    result = run_database_cli(
        tmp_path, ["-db", db, "--references", str(references), "--annotation", str(snp)]
    )
    assert result.returncode == 0, result.stderr
    assert count(db, "snp_annotation", "snp_id = 'B.1'") == 1


def test_failed_annotation_run_leaves_the_database_as_it_was(tmp_path):
    db = database(tmp_path)
    snp = tmp_path / "new_snp.txt"
    snp.write_text(SNP_LINE + "Z\tfrancisella\tBirdsell2014\tREF1\tCP000001.1\t300\tT\tC\n")
    before = open(db, "rb").read()
    result = run_database_cli(tmp_path, ["-db", db, "--annotation", str(snp)])
    assert result.returncode == 1
    assert "'Z'" in result.stderr + result.stdout
    assert open(db, "rb").read() == before


def test_nothing_to_modify_names_the_options(tmp_path):
    db = database(tmp_path)
    before = open(db, "rb").read()
    result = run_database_cli(tmp_path, ["-db", db])
    assert result.returncode == 1
    text = result.stderr + result.stdout
    for option in ("--mod-file", "--cansnps", "--references"):
        assert option in text
    assert "both mod_database and mod_file" not in text
    assert open(db, "rb").read() == before
