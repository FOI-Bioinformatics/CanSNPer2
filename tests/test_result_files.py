"""The per-sample result files and snps_summary.txt.

<sample>_snps.txt and <sample>_not_called.txt start with the table (a header row and the rows);
the SNP path, the final SNP, the query file and the typing settings follow as # comment lines,
so that a table reader that skips comment lines (pandas.read_csv(comment="#")) reads the table.
Result files of earlier versions wrote these lines without the #; they are still read, so a
result of an earlier run is used again as before. snps_summary.txt starts with # lines on the
run (version, date, database, scheme, settings) and a header row."""

import csv
import os

import pytest

from test_typing_inputs import project, write, fake_alignment
from CanSNPer2.modules.CanSNPer2 import SNP_FILE_COLUMNS


def table_rows(path):
    """the rows of a result file as a reader that skips # lines sees them"""
    with open(path) as handle:
        return list(csv.reader((row for row in handle if not row.startswith("#")), delimiter="\t"))


def comments(path):
    return [row.rstrip("\n") for row in open(path) if row.startswith("#")]


def typed(tmp_path, called, states, **kwargs):
    q = write(tmp_path / "x.fasta")
    obj, db = project(tmp_path, [q], min_required_hits=2, **kwargs)
    fake_alignment(obj, called, states, q)
    obj.run(database=db)
    return obj, q


def test_the_snps_file_is_a_table_followed_by_comment_lines(tmp_path):
    obj, q = typed(tmp_path, ["A.1", "B.1"], {"A.1": 1, "B.1": 1})
    path = os.path.join(obj.results_dir, "x_snps.txt")
    rows = table_rows(path)
    assert rows[0] == SNP_FILE_COLUMNS
    assert [row[0] for row in rows[1:]] == ["A.1", "B.1"]
    assert all(len(row) == len(SNP_FILE_COLUMNS) for row in rows)
    lines = open(path).read().splitlines()
    first_comment = next(i for i, line in enumerate(lines) if line.startswith("#"))
    assert all(line.startswith("#") for line in lines[first_comment:]), "the comment lines follow the table"
    notes = comments(path)
    assert "# SNP path: A.1;B.1" in notes
    assert any(line.startswith("# Final SNP: B.1 found/depth: 2/") for line in notes)
    assert any(line.startswith("# Query file: " + os.path.abspath(q)) for line in notes)
    assert any(line.startswith("# Typing settings: database sha256:") for line in notes)


def test_the_snps_file_of_a_sample_without_a_call_keeps_its_header_row(tmp_path):
    obj, q = typed(tmp_path, [], {"A.1": 2})
    path = os.path.join(obj.results_dir, "x_snps.txt")
    assert table_rows(path) == [SNP_FILE_COLUMNS]
    assert "# SNP path: NA" in comments(path)
    assert "# Final SNP: NA" in comments(path)


def test_the_not_called_file_is_a_table_followed_by_comment_lines(tmp_path):
    obj, q = typed(tmp_path, ["A.1", "B.1"], {"A.1": 1, "B.1": 1, "C.1": 2})
    path = os.path.join(obj.results_dir, "x_not_called.txt")
    rows = table_rows(path)
    assert rows[0] == SNP_FILE_COLUMNS
    ## B.2 is in no alignment, so it is listed as a position without data (fake_alignment gives
    ## rows for the called SNPs only, so C.1 has no row of its own here)
    assert [row[0] for row in rows[1:]] == ["B.2"]
    assert rows[1][SNP_FILE_COLUMNS.index("Target base")] == "N/A"
    lines = open(path).read().splitlines()
    first_comment = next(i for i, line in enumerate(lines) if line.startswith("#"))
    assert all(line.startswith("#") for line in lines[first_comment:])
    assert any(line.startswith("# Query file: ") for line in comments(path))
    assert any(line.startswith("# Typing settings: ") for line in comments(path))


OLD = (
    "SNP path: A.1;B.1\nFinal SNP: B.1 found/depth: 2/3\nQuery file: /data/x.fasta sha256:{sha}\n{settings}\n"
)
NEW = "# SNP path: A.1;B.1\n# Final SNP: B.1 found/depth: 2/3\n# Query file: /data/x.fasta sha256:{sha}\n# {settings}\n"


@pytest.mark.parametrize("form", [OLD, NEW], ids=["without #", "with #"])
def test_both_forms_of_a_result_file_are_read(tmp_path, form):
    obj, db = project(tmp_path, [])
    path = str(tmp_path / "x_snps.txt")
    settings = obj.settings_line()
    with open(path, "w") as out:
        out.write("\t".join(SNP_FILE_COLUMNS) + "\n" + form.format(sha="ab" * 32, settings=settings))
    assert obj.read_result_file(path) == ("B.1", "A.1;B.1")
    assert obj.read_result_source(path) == ("/data/x.fasta", "ab" * 32)
    assert obj.read_result_settings(path) == settings
    assert obj.read_result_support(path) == "2/3"


@pytest.mark.parametrize("prefix", ["", "# "], ids=["without #", "with #"])
def test_a_result_of_an_earlier_run_is_used_again_in_both_forms(tmp_path, prefix):
    q = write(tmp_path / "in" / "x.fasta")
    obj, db = project(tmp_path, [q], keep_going=True)
    os.makedirs(obj.results_dir, exist_ok=True)
    with open(os.path.join(obj.results_dir, "x_snps.txt"), "w") as out:
        out.write("\t".join(SNP_FILE_COLUMNS) + "\n")
        for line in (
            "SNP path: A.1",
            "Final SNP: A.1 found/depth: 1/2",
            "Query file: /elsewhere/x.fasta sha256:" + obj.file_sha256(q),
            obj.settings_line(),
        ):
            out.write(prefix + line + "\n")
    obj.run(database=db)
    assert "x" in obj.reused
    report = obj.sample_reports["x"]
    assert report["final"] == "A.1" and report["support"] == "1/2"


def test_the_summary_starts_with_the_run_and_a_header_row(tmp_path):
    obj, q = typed(tmp_path, ["A.1", "B.1"], {"A.1": 1, "B.1": 1}, summary=True)
    path = os.path.join(obj.outdir, "snps_summary.txt")
    notes = comments(path)
    text = "\n".join(notes)
    assert notes and all(line.startswith("# ") for line in notes)
    assert "CanSNPer2 version" in text
    assert "Date" in text
    assert "Database: t.db sha256:" in text
    assert "Scheme" in text
    assert "min_required_hits=2" in text and "strictness=0.7" in text and "mask=0" in text
    lines = open(path).read().splitlines()
    assert lines[: len(notes)] == notes, "the comment lines come first"
    rows = table_rows(path)
    assert rows[0] == ["sample", "final_snp", "support", "snp_path"]
    assert rows[1] == ["x", "B.1", "2/2", "A.1;B.1"]


def test_the_summary_gives_na_support_for_a_sample_without_a_call(tmp_path):
    obj, q = typed(tmp_path, [], {"A.1": 2}, summary=True)
    rows = table_rows(os.path.join(obj.outdir, "snps_summary.txt"))
    assert rows[1] == ["x", "NA", "NA", "NA"]


# -- the per-sample tree page is opt-in now, and every per-sample file sits in
#    the samples subfolder of --outdir so that summary.html is the one file to find --


def test_the_tree_page_is_left_out_by_default(tmp_path):
    """the results overview shows the tree of every sample, so the standalone pages
    are back to being written on request (--save_tree)"""
    obj, q = typed(tmp_path, ["A.1", "B.1"], {"A.1": 1, "B.1": 1})
    assert not os.path.exists(os.path.join(obj.results_dir, "x_tree.html"))


def test_save_tree_writes_the_tree_page_on_request(tmp_path):
    obj, q = typed(tmp_path, ["A.1", "B.1"], {"A.1": 1, "B.1": 1}, save_tree=True)
    assert os.path.exists(os.path.join(obj.results_dir, "x_tree.html"))
    assert os.path.exists(os.path.join(obj.results_dir, "x_snps.txt")), "the tables are written either way"


def test_per_sample_files_sit_in_a_samples_subfolder(tmp_path):
    """summary.html greets you in --outdir, its links lead into samples/; a run with no
    sample output at all leaves no empty folder behind"""
    obj, q = typed(tmp_path, ["A.1", "B.1"], {"A.1": 1, "B.1": 1})
    assert os.path.exists(os.path.join(obj.outdir, "samples", "x_snps.txt"))
    assert os.path.exists(os.path.join(obj.outdir, "summary.html"))
    assert not os.path.exists(os.path.join(obj.outdir, "x_snps.txt")), "nothing per sample in --outdir"


def test_a_sample_whose_result_exists_prints_its_call_and_how_to_redo_it(tmp_path, capsys):
    """The skip was a debug line, so a second run into the same folder announced the
    alignments and then said nothing at all (issue #78)."""
    obj, q = typed(tmp_path, ["A.1", "B.1"], {"A.1": 1, "B.1": 1})
    first = capsys.readouterr().out
    assert "x: B.1" in first

    again, db = project(tmp_path, [q], min_required_hits=2)
    again.outdir = obj.outdir
    fake_alignment(again, ["A.1", "B.1"], {"A.1": 1, "B.1": 1}, q)
    again.run(database=db)
    printed = capsys.readouterr().out
    assert "x: B.1" in printed, "the call of the earlier run is reported"
    assert "--rerun" in printed, "and what to do to type it again"


# the "Final SNP:" line gives found/depth (kept for scripts that read it) and found/expected,
# the support; the root of the test tree has no canSNP, a label node, so the two differ (#99)


def test_the_final_snp_line_gives_the_depth_and_the_support(tmp_path):
    obj, q = typed(tmp_path, ["A.1", "B.1"], {"A.1": 1, "B.1": 1})
    path = os.path.join(obj.results_dir, "x_snps.txt")
    final = [line for line in comments(path) if line.startswith("# Final SNP:")]
    assert final == ["# Final SNP: B.1 found/depth: 2/3 found/expected: 2/2"]
    assert obj.read_result_support(path) == "2/2"
    assert obj.sample_supports["x"] == "2/2"


def test_a_reused_result_gives_the_same_support_in_the_summary(tmp_path):
    obj, q = typed(tmp_path, ["A.1", "B.1"], {"A.1": 1, "B.1": 1}, summary=True)
    first = table_rows(os.path.join(obj.outdir, "snps_summary.txt"))
    again, db = project(tmp_path, [q], min_required_hits=2, summary=True)
    fake_alignment(again, ["A.1", "B.1"], {"A.1": 1, "B.1": 1}, q)
    again.run(database=db)
    assert "x" in again.reused
    second = table_rows(os.path.join(again.outdir, "snps_summary.txt"))
    assert second == first
    assert second[1][2] == "2/2"


# competing paths (#91): B.1 and B.2 are both children of A.1, at the same depth


def test_competing_paths_are_written_to_the_result_file(tmp_path):
    obj, q = typed(tmp_path, ["A.1", "B.1", "B.2"], {"A.1": 1, "B.1": 1, "B.2": 1})
    assert obj.not_typed["x"].startswith("Two or more confirmed SNP paths compete")
    notes = comments(os.path.join(obj.results_dir, "x_snps.txt"))
    assert "# Final SNP: NA found/depth: 2/3" in notes
    assert "# Competing paths: B.1 found/expected: 2/2; B.2 found/expected: 2/2" in notes
    assert notes.index("# Competing paths: B.1 found/expected: 2/2; B.2 found/expected: 2/2") == 2


def test_a_called_sample_has_no_competing_paths_line(tmp_path):
    obj, q = typed(tmp_path, ["A.1", "B.1"], {"A.1": 1, "B.1": 1, "B.2": 2})
    notes = comments(os.path.join(obj.results_dir, "x_snps.txt"))
    assert not [line for line in notes if line.startswith("# Competing paths")]
