"""The per-sample result files.

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
	'''the rows of a result file as a reader that skips # lines sees them'''
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
	path = os.path.join(obj.outdir, "x_snps.txt")
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
	path = os.path.join(obj.outdir, "x_snps.txt")
	assert table_rows(path) == [SNP_FILE_COLUMNS]
	assert "# SNP path: NA" in comments(path)
	assert "# Final SNP: NA" in comments(path)


def test_the_not_called_file_is_a_table_followed_by_comment_lines(tmp_path):
	obj, q = typed(tmp_path, ["A.1", "B.1"], {"A.1": 1, "B.1": 1, "C.1": 2})
	path = os.path.join(obj.outdir, "x_not_called.txt")
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


OLD = "SNP path: A.1;B.1\nFinal SNP: B.1 found/depth: 2/3\nQuery file: /data/x.fasta sha256:{sha}\n{settings}\n"
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
	os.makedirs(obj.outdir, exist_ok=True)
	with open(os.path.join(obj.outdir, "x_snps.txt"), "w") as out:
		out.write("\t".join(SNP_FILE_COLUMNS) + "\n")
		for line in ("SNP path: A.1", "Final SNP: A.1 found/depth: 1/2",
				"Query file: /elsewhere/x.fasta sha256:" + obj.file_sha256(q), obj.settings_line()):
			out.write(prefix + line + "\n")
	obj.run(database=db)
	assert "x" in obj.reused
	assert obj.sample_reports["x"]["final"] == "A.1"
