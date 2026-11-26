import sqlite3

import pytest

import helpers
from CanSNPer2.modules.DatabaseConnection import LoadCanSNPAnnotation, SourceFileError


def make_db(path):
    """the nodes N1 and N2 without any annotation, for the loader to fill"""
    helpers.make_db(path, [])
    conn = sqlite3.connect(str(path))
    conn.execute("INSERT INTO nodes (id, name) VALUES (1, 'N1'), (2, 'N2')")
    conn.commit()
    conn.close()


def make_file(path, rows):
    with open(path, "w") as f:
        f.write(
            "#node_name\tsnp_id\tsequence\tsequence_position\tancestral_base\tderived_base\treference\tdate\tgenome_i\n"
        )
        for row in rows:
            f.write(row + "\n")


def test_load_annotations_loads_all_rows(tmp_path):
    db, ann = tmp_path / "test.db", tmp_path / "ann.tsv"
    make_db(db)
    make_file(
        ann, ["N1\t7\tchrI\t100\tA\tG\tref1\t2026-01-01\t5", "N2\t3\tchrII\t200\tC\tT\tref2\t2026-01-02\t6"]
    )
    ldr = LoadCanSNPAnnotation(str(db))
    ldr.load_annotations(str(ann))
    ldr.commit()
    rows = ldr.query(
        "SELECT node_id, snp_id, sequence, sequence_position, genome_i FROM snp_annotation ORDER BY node_id"
    ).fetchall()
    ldr.disconnect()
    assert rows == [(1, "7", "chrI", 100, 5), (2, "3", "chrII", 200, 6)]


def test_unknown_node_warns_and_skips(tmp_path, caplog):
    db, ann = tmp_path / "test.db", tmp_path / "ann.tsv"
    make_db(db)
    make_file(
        ann, ["NOPE\t7\tchrI\t100\tA\tG\tref1\t2026-01-01\t5", "N1\t8\tchrI\t150\tA\tT\tref1\t2026-01-01\t5"]
    )
    ldr = LoadCanSNPAnnotation(str(db))
    ldr.load_annotations(str(ann))
    ldr.commit()
    rows = ldr.query("SELECT node_id, snp_id FROM snp_annotation").fetchall()
    ldr.disconnect()
    assert rows == [(1, "8")]
    assert "NOPE was not found in the tree" in caplog.text


def test_short_row_warns_and_skips(tmp_path, caplog):
    db, ann = tmp_path / "test.db", tmp_path / "ann.tsv"
    make_db(db)
    make_file(ann, ["N1\t9\t150\tA", "N1\t8\tchrI\t150\tA\tT\tref1\t2026-01-01\t5"])
    ldr = LoadCanSNPAnnotation(str(db))
    ldr.load_annotations(str(ann))
    ldr.commit()
    rows = ldr.query("SELECT node_id, snp_id FROM snp_annotation").fetchall()
    ldr.disconnect()
    assert rows == [(1, "8")]
    assert "Skipped short row" in caplog.text


def test_missing_annotation_file_raises(tmp_path):
    db = tmp_path / "test.db"
    make_db(db)
    ldr = LoadCanSNPAnnotation(str(db))
    for missing in (None, str(tmp_path / "nope.tsv")):
        with pytest.raises(SourceFileError):
            ldr.load_annotations(missing)
    ldr.disconnect()
