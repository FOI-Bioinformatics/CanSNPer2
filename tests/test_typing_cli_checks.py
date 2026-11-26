"""Typing options and run checks found by the 2.3.0 stress test (items 23 to 29 of its list):
option values, run directories left behind, --workdir, --skip_alignment against a different
query, sample names, query messages and a warning repeated for every sample."""
import logging
import os
import subprocess
import sys

import pytest

from CanSNPer2.modules.CanSNPer2 import CanSNPer2, CanSNPer2Error
from test_aligners import ToyAligner, toy_case
from test_scheme_info import ROOT, TESTS, TYPE_SCRIPT, toy_database


def cli(tmp_path, *arguments, tmpdir=True):
	env = dict(os.environ)
	env["PYTHONPATH"] = os.pathsep.join([ROOT, TESTS, str(tmp_path), env.get("PYTHONPATH", "")])
	extra = ["--tmpdir", str(tmp_path / "tmp")] if tmpdir else []
	return subprocess.run([sys.executable, "-c", TYPE_SCRIPT, "--aligner", "toy", "--refdir", "refs"] + extra + list(arguments),
		capture_output=True, text=True, cwd=str(tmp_path), env=env)


def run_dirs(tmp_path):
	base = tmp_path / "tmp" / "CanSNPer2"
	return sorted(p for p in os.listdir(base) if p.startswith("run_")) if base.is_dir() else []


# 23: option values

@pytest.mark.parametrize("option,value", [("--alignment_timeout", "-1"), ("--alignment_timeout", "nan"),
	("--min_required_hits", "0"), ("--min_required_hits", "-3"), ("--mask", "-50")])
def test_option_values_out_of_range_are_refused(tmp_path, option, value):
	db = toy_database(tmp_path)
	result = cli(tmp_path, "--database", db, option, value, "sample.fna")
	assert result.returncode == 2, result.stderr
	assert option in result.stderr
	assert run_dirs(tmp_path) == []


def test_min_required_hits_help_describes_what_is_counted():
	from CanSNPer2 import CanSNPerTree
	source = open(CanSNPerTree.__file__).read()
	assert "Minimum sequential hits" not in source


# 24: run directories

def test_no_query_is_an_error_and_leaves_no_run_directory(tmp_path):
	db = toy_database(tmp_path)
	result = cli(tmp_path, "--database", db)
	assert result.returncode == 1
	assert "no query" in result.stderr.lower()
	assert run_dirs(tmp_path) == []


def test_an_unexpected_error_removes_the_run_directory(tmp_path, monkeypatch):
	obj, db, query = toy_case(tmp_path)
	monkeypatch.setattr(obj, "check_start", lambda: (_ for _ in ()).throw(RuntimeError("boom")))
	with pytest.raises(RuntimeError):
		obj.run(database=db)
	assert not os.path.exists(obj.tmpdir)


def test_the_kept_message_is_printed_only_when_files_are_kept(tmp_path):
	db = toy_database(tmp_path)
	result = cli(tmp_path, "--database", db, "--min_required_hits", "1", "-o", "out", "sample.fna")
	assert result.returncode == 0, result.stderr
	assert "are kept in" not in result.stdout
	assert run_dirs(tmp_path) == []
	result = cli(tmp_path, "--database", db, "--min_required_hits", "1", "-o", "out2", "--keep_temp", "sample.fna")
	assert "kept in" in result.stdout + result.stderr


# 25: --workdir

def test_workdir_holds_the_log_and_tmp_folders(tmp_path):
	db = toy_database(tmp_path)
	(tmp_path / "wd").mkdir()
	os.rename(tmp_path / "refs", tmp_path / "wd" / "refs")
	result = cli(tmp_path, "--database", db, "--min_required_hits", "1", "--workdir", "wd", "-o", "res",
		"--logdir", "lg", "--tmpdir", "tmpA", str(tmp_path / "sample.fna"), tmpdir=False)
	assert result.returncode == 0, result.stderr
	assert not (tmp_path / "lg").exists() and not (tmp_path / "tmpA").exists()
	assert any(name.startswith("CanSNPer2") for name in os.listdir(tmp_path / "wd" / "lg"))
	assert (tmp_path / "wd" / "res" / "samples" / "sample_snps.txt").exists()


# 26: --skip_alignment

def test_skip_alignment_refuses_alignments_of_another_query(tmp_path):
	obj, db, query = toy_case(tmp_path, keep_temp=True)
	obj.run(database=db)
	kept = obj.tmpdir
	query.write_text(">q\nACGTACGA\n")  # same name, other content
	obj = CanSNPer2([str(query)], workdir=str(tmp_path), refdir="refs", tmpdir=kept, logdir="logs", outdir="out2",
		database=db, aligner=ToyAligner(), skip_alignment=True, keep_going=True)
	obj.run(database=db)
	assert "sample" not in obj.typed
	assert "different query" in obj.not_typed["sample"]


def test_skip_alignment_with_the_same_query_is_typed(tmp_path):
	obj, db, query = toy_case(tmp_path, keep_temp=True)
	obj.run(database=db)
	obj = CanSNPer2([str(query)], workdir=str(tmp_path), refdir="refs", tmpdir=obj.tmpdir, logdir="logs", outdir="out2",
		database=db, aligner=ToyAligner(), skip_alignment=True, min_required_hits=1)
	obj.run(database=db)
	assert "sample" in obj.typed


def test_missing_alignments_stop_the_run_without_keep_going(tmp_path, caplog):
	obj, db, query = toy_case(tmp_path)
	empty = tmp_path / "empty_run"
	empty.mkdir()
	obj = CanSNPer2([str(query)], workdir=str(tmp_path), refdir="refs", tmpdir=str(empty), logdir="logs", outdir="out2",
		database=db, aligner=ToyAligner(), skip_alignment=True)
	with pytest.raises(CanSNPer2Error):
		obj.run(database=db)
	assert "xmfa" not in caplog.text


def test_skip_alignment_with_a_missing_tmpdir_does_not_create_it(tmp_path):
	db = toy_database(tmp_path)
	result = cli(tmp_path, "--database", db, "--skip_alignment", "--tmpdir", "nowhere", "sample.fna", tmpdir=False)
	assert result.returncode == 1
	assert not (tmp_path / "nowhere").exists()
	assert "nowhere" in result.stderr


# 27: sample names

def test_an_empty_sample_name_is_refused(tmp_path):
	obj, db, query = toy_case(tmp_path)
	hidden = tmp_path / ".fna"
	hidden.write_text(">q\nACTTACGT\n")
	with pytest.raises(CanSNPer2Error, match="empty sample name"):
		obj.check_sample_names([str(hidden)])


def test_a_sample_named_summary_is_refused_with_summary_and_save_tree(tmp_path):
	obj, db, query = toy_case(tmp_path, summary=True, save_tree=True)
	named = tmp_path / "summary.fna"
	named.write_text(">q\nACTTACGT\n")
	with pytest.raises(CanSNPer2Error, match="summary"):
		obj.check_sample_names([str(named)])


# 28: query messages

def test_a_query_name_starting_with_a_dash_is_not_passed_as_an_option(tmp_path, monkeypatch):
	obj, db, query = toy_case(tmp_path)
	monkeypatch.chdir(tmp_path)
	(tmp_path / "-x.fna").write_text(">q\nACTTACGT\n")
	obj.query_name = "-x"
	commands, logs, names = obj.create_alignment_commands("-x.fna")
	assert not any(arg.startswith("-x") for command in commands for arg in command)


def test_tmpdir_root_is_not_called_empty(tmp_path):
	db = toy_database(tmp_path)
	result = cli(tmp_path, "--database", db, "--tmpdir", "/", "sample.fna", tmpdir=False)
	assert "cannot be empty" not in result.stderr


def test_a_directory_as_query_is_named(tmp_path):
	obj, db, query = toy_case(tmp_path, keep_going=True)
	(tmp_path / "q").mkdir()
	obj.query = [str(tmp_path / "q")]
	obj.run(database=db)
	assert "is a directory" in obj.not_typed["q"]


def test_a_missing_query_stops_the_run_before_any_alignment(tmp_path):
	obj, db, query = toy_case(tmp_path)
	obj.query = [str(query), str(tmp_path / "missing.fna")]
	with pytest.raises(CanSNPer2Error, match="missing.fna"):
		obj.run(database=db)
	assert not [f for f in os.listdir(obj.tmpdir) if f.endswith(".tsv")] if os.path.isdir(obj.tmpdir) else True


# 29: a warning repeated for every sample

def test_the_root_without_a_snp_is_warned_about_once(tmp_path, caplog):
	obj, db, query = toy_case(tmp_path)
	second = tmp_path / "second.fna"
	second.write_text(">q\nACTTACGT\n")
	obj.query = [str(query), str(second)]
	caplog.set_level(logging.WARNING)
	obj.run(database=db)
	assert caplog.text.count("has no snp_annotation root") <= 1


# moved from tests/release_test.sh: --keep-going must not hide a run that typed nothing
# (issues #43 and #191)

def test_keep_going_with_every_sample_failed_exits_1_and_names_the_sample(tmp_path):
	db = toy_database(tmp_path)
	(tmp_path / "q").mkdir()  # a directory as query fails the sample, --keep_going goes on
	result = cli(tmp_path, "--database", db, "--keep_going", "-o", "out", "q")
	assert result.returncode == 1, result.stderr
	assert "No sample was typed" in result.stderr
	assert "q: " in result.stderr.split("No sample was typed", 1)[1]
	assert sorted(os.listdir(tmp_path / "out")) == ["summary.html"]
	assert '"q"' in (tmp_path / "out" / "summary.html").read_text()
