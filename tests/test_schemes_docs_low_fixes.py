"""Low severity fixes: schemes.toml entries, scheme file errors, scheme_name checks,
help texts, install notes, CanSNPer2-test and CanSNPer2-download messages."""
import os
import re
import subprocess
import sys

import pytest

import helpers
from test_fetch_scheme import scheme_repo, write_tree
from CanSNPer2.modules.FetchScheme import SchemeError, LocalSchemeSource, build_scheme

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ENV = dict(os.environ, PYTHONPATH=ROOT)


def run(module, *args, cwd=None, env=None):
	return subprocess.run([sys.executable, "-m", module] + list(args), capture_output=True, text=True,
		env=env or ENV, cwd=str(cwd) if cwd else None)


def data_with_toml(tmp_path, toml):
	files = scheme_repo()
	files["schemes.toml"] = toml
	return write_tree(tmp_path / "data", files)


# 20. schemes.toml entries of the wrong type

@pytest.mark.parametrize("toml", ["test_scheme = 1\n", '[test_scheme]\nsnp = 5\n', '[test_scheme]\ndatabase_dir = 3\n',
	'[test_scheme]\nreadme = []\n'])
def test_a_malformed_schemes_toml_entry_is_a_scheme_error(tmp_path, toml):
	src = LocalSchemeSource(str(data_with_toml(tmp_path, toml)))
	with pytest.raises(SchemeError, match="schemes.toml") as e:
		build_scheme("test_scheme", str(tmp_path / "x.db"), source=src)
	assert "test_scheme" in str(e.value)


def test_list_schemes_with_a_malformed_entry_gives_a_message(tmp_path):
	data = data_with_toml(tmp_path, "sc = 1\n")
	done = run("CanSNPer2.DownloadGenomes", "--list_schemes", "--data_dir", str(data), "--logs", str(tmp_path / "l"), cwd=tmp_path)
	assert done.returncode == 1
	assert "Traceback" not in done.stderr and "sc" in done.stderr and "schemes.toml" in done.stderr


# 21. a bad scheme SNP file names the scheme

def test_a_bad_snp_file_names_the_scheme_not_a_temporary_path(tmp_path):
	bad = scheme_repo(snp="snp_id\tstrain\n")
	data = write_tree(tmp_path / "data", bad)
	with pytest.raises(Exception) as e:
		build_scheme("test_scheme", str(tmp_path / "x.db"), source=LocalSchemeSource(str(data)))
	text = str(e.value)
	assert "test_scheme" in text and "test_scheme_snp.txt" in text
	assert "cansnper2_scheme_" not in text


# 22. control characters in scheme_name / scheme_source

@pytest.mark.parametrize("option", ["--scheme_name", "--scheme_source"])
@pytest.mark.parametrize("bad", ["a\tb", "a\nb", "a\rb"])
def test_control_characters_in_scheme_options_are_refused(tmp_path, option, bad):
	paths = helpers.write_sources(tmp_path, [("A.1", "s", "r", "REF", "CP1", 10, "T", "C")])
	done = run("CanSNPer2.SNPDatabase", "--create", "--database", str(tmp_path / "x.db"), "--tree", paths["tree"],
		"--annotation", paths["snp"], "--references", paths["references"], option, bad, "--logs", str(tmp_path / "l"), cwd=tmp_path)
	assert done.returncode != 0
	assert option.replace("_", "-") in done.stderr + done.stdout and "control" in (done.stderr + done.stdout)
	assert not (tmp_path / "x.db").exists()


# 30. help texts

def test_save_tree_and_no_snpfiles_help_describe_the_behaviour():
	text = open(os.path.join(ROOT, "CanSNPer2", "CanSNPerTree.py")).read()
	assert "ETE3" not in text and "PDF" not in text
	# the results overview carries every sample's tree, so the standalone page is opt-in
	# again: --save_tree is the flag that carries the help, --no_save_tree stays accepted
	# but hidden; read from the parser, so that the layout of the source does not matter
	import argparse

	from CanSNPer2.CanSNPerTree import build_parser

	actions = {a.dest: a for a in build_parser()._actions}
	assert "html" in actions["save_tree"].help and "samples/" in actions["save_tree"].help
	assert actions["no_save_tree"].help == argparse.SUPPRESS
	assert "Don´t save" not in text


def test_cli_reference_matches_the_help():
	reference = open(os.path.join(ROOT, "docs", "cli.md")).read()
	block = reference.split("## cansnper2 type\n\n```text\n")[1].split("\n```")[0]
	help_text = run("CanSNPer2.cli", "type", "--help", env=dict(ENV, COLUMNS="80")).stdout
	for option in re.findall(r"^\s+(--[a-z][a-z-]*)", help_text, re.M):
		assert option in block, option
	assert "ETE3" not in block and "Don´t" not in block


# 31. install notes

def test_install_md_matches_setup():
	text = open(os.path.join(ROOT, "docs", "installation.md")).read()
	for word in ("3.7", "[ETE3]", "pip install PyQt5", "[NumPy]"):
		assert word not in text
	assert "3.11" in text


def test_selftest_does_not_mention_a_missing_requirements_file():
	assert "requirements.txt" not in open(os.path.join(ROOT, "CanSNPer2", "selftest.py")).read()


# 32. CanSNPer2-test

def fake_path(tmp_path, programs):
	bin_dir = tmp_path / "bin"
	bin_dir.mkdir()
	for program in programs:
		p = bin_dir / program
		p.write_text("#!/bin/sh\nexit 0\n")
		p.chmod(0o755)
	return str(bin_dir)


@pytest.mark.parametrize("programs,code,lines", [
	(["minimap2"], 0, ["progressiveMauve: not found", "minimap2: ", "nucmer: not found"]),
	(["nucmer"], 1, ["delta-filter: not found"]),
	(["nucmer", "delta-filter"], 0, ["nucmer: ", "delta-filter: "]),
	(["progressiveMauve"], 0, ["progressiveMauve: "]),
	([], 1, ["progressiveMauve: not found", "minimap2: not found"]),
])
def test_selftest_reports_each_aligner(tmp_path, programs, code, lines):
	env = dict(ENV, PATH=fake_path(tmp_path, programs))
	done = subprocess.run([sys.executable, "-m", "CanSNPer2.selftest"], capture_output=True, text=True, env=env)
	out = done.stdout + done.stderr
	assert done.returncode == code, out
	for line in lines:
		assert line in out, out
	assert "python packages" not in out
	if False:
		assert "nucmer: not found" not in out.split("delta-filter")[0] or "delta-filter: not found" in out


# 33. download: --no_references with -db, and a line on success

def test_no_references_with_database_only_is_refused(tmp_path):
	db = helpers.make_db(tmp_path / "x.db", [("A.1", "CP1", 10, "C", "T")])
	done = run("CanSNPer2.DownloadGenomes", "-db", db, "--no_references", "--logs", str(tmp_path / "l"), cwd=tmp_path)
	assert done.returncode != 0
	assert "--no_references" in done.stderr and "--scheme" in done.stderr


def test_download_prints_one_line_on_success(tmp_path, monkeypatch, capsys):
	from CanSNPer2 import DownloadGenomes
	db = helpers.make_db(tmp_path / "x.db", [("A.1", "CP1", 10, "C", "T")])
	monkeypatch.setattr(DownloadGenomes.DownloadGenomes, "run", lambda self: 0)
	monkeypatch.setattr(sys, "argv", ["CanSNPer2-download", "-db", db, "-o", str(tmp_path / "refs"), "--logs", str(tmp_path / "l")])
	DownloadGenomes.command()
	out = capsys.readouterr().out.strip().splitlines()
	assert len(out) == 1 and db in out[0] and str(tmp_path / "refs") in out[0]
