"""The scheme_info table: which scheme a database holds, its status, where it came from and
which CanSNPer2 version built it, and the notice CanSNPer2 logs when it opens a database."""
import datetime
import logging
import os
import sqlite3
import subprocess
import sys

import pytest

import helpers
from test_fetch_scheme import COMMIT, FakeRepo, scheme_repo
from CanSNPer2 import __version__
from CanSNPer2.modules.FetchScheme import DATA_REPO, GitHubSchemeSource, build_scheme, provenance_path
from CanSNPer2.modules.SchemeInfo import (MISSING_TEXT, describe, read_scheme_info, scheme_notice,
	write_scheme_info)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TESTS = os.path.dirname(os.path.abspath(__file__))


def run(tmp_path, module, arguments, extra_path=()):
	env = dict(os.environ)
	env["PYTHONPATH"] = os.pathsep.join([ROOT, TESTS] + list(extra_path) + [env.get("PYTHONPATH", "")])
	return subprocess.run([sys.executable, "-m", module] + arguments, capture_output=True, text=True,
		cwd=str(tmp_path), env=env)


def create(tmp_path, db, *extra):
	paths = helpers.write_sources(tmp_path / "src", ["A.1\tfrancisella\tBirdsell2014\tREF1\tCP000001.1\t100\tT\tC",
		"B.1\tfrancisella\tBirdsell2014\tREF1\tCP000001.1\t200\tG\tA"])
	return run(tmp_path, "CanSNPer2.SNPDatabase", ["--create", "--database", db, "--tree", paths["tree"],
		"--annotation", paths["snp"], "--references", paths["references"], "--logs", "logs"] + list(extra)), paths


def test_create_writes_format_and_version(tmp_path):
	db = str(tmp_path / "s.db")
	result, _paths = create(tmp_path, db, "--scheme_name", "my_scheme", "--scheme_status", "draft",
		"--scheme_source", "my-folder")
	assert result.returncode == 0, result.stderr
	assert read_scheme_info(db) == {"format": "2.3", "built_with": __version__,
		"built": datetime.date.today().isoformat(), "scheme": "my_scheme", "status": "draft", "source": "my-folder"}


def test_create_without_scheme_options_writes_the_build_keys_only(tmp_path):
	db = str(tmp_path / "s.db")
	result, _paths = create(tmp_path, db)
	assert result.returncode == 0, result.stderr
	assert sorted(read_scheme_info(db)) == ["built", "built_with", "format"]


def test_import_writes_format_and_version_and_scheme_options(tmp_path):
	db = str(tmp_path / "s.db")
	result, paths = create(tmp_path, db)
	assert result.returncode == 0, result.stderr
	imported = str(tmp_path / "i.db")
	result = run(tmp_path, "CanSNPer2.SNPDatabase", ["--import", "--database", imported, "--tree", paths["tree"],
		"--annotation", paths["snp"], "--references", paths["references"], "--scheme_name", "my_scheme",
		"--scheme_status", "verified", "--logs", "logs"])
	assert result.returncode == 0, result.stderr
	info = read_scheme_info(imported)
	assert info["format"] == "2.3" and info["built_with"] == __version__
	assert info["scheme"] == "my_scheme" and info["status"] == "verified"


def test_an_unknown_scheme_status_is_refused(tmp_path):
	db = str(tmp_path / "s.db")
	result, _paths = create(tmp_path, db, "--scheme_status", "tested")
	assert result.returncode == 1
	assert "--scheme_status" in result.stderr and "verified, draft or unknown" in result.stderr
	assert "Traceback" not in result.stderr
	assert not os.path.exists(db)


def test_scheme_options_are_refused_without_create_or_import(tmp_path):
	db = helpers.make_db(tmp_path / "s.db", [("A.1", "chr", 3, "G", "T")])
	result = run(tmp_path, "CanSNPer2.SNPDatabase", ["--database", db, "--export", "-o", "out",
		"--scheme_name", "x", "--logs", "logs"])
	assert result.returncode == 1
	assert "--create" in result.stderr and "Traceback" not in result.stderr


def test_info_prints_the_table(tmp_path):
	db = str(tmp_path / "s.db")
	result, _paths = create(tmp_path, db, "--scheme_name", "my_scheme", "--scheme_status", "verified")
	assert result.returncode == 0, result.stderr
	result = run(tmp_path, "CanSNPer2.SNPDatabase", ["--database", db, "--info", "--logs", "logs"])
	assert result.returncode == 0, result.stderr
	lines = result.stdout.splitlines()
	assert "scheme\tmy_scheme" in lines and "status\tverified" in lines
	assert "format\t2.3" in lines and "built_with\t{v}".format(v=__version__) in lines
	assert any(line.startswith("built\t") for line in lines)


def test_info_of_a_database_without_the_table(tmp_path):
	db = helpers.make_db(tmp_path / "old.db", [("A.1", "chr", 3, "G", "T")])
	result = run(tmp_path, "CanSNPer2.SNPDatabase", ["--database", db, "--info", "--logs", "logs"])
	assert result.returncode == 0, result.stderr
	assert MISSING_TEXT in result.stdout


def test_scheme_download_records_status_commit_and_known_issues(tmp_path):
	files = scheme_repo(status="draft")
	files["schemes.toml"] += '[test_scheme]\nknown_issues = ["B.1", "A.1"]\n'
	db = str(tmp_path / "test_scheme.db")
	build_scheme("test_scheme", db, source=GitHubSchemeSource(ref="dev", fetcher=FakeRepo(files)))
	info = read_scheme_info(db)
	assert info["scheme"] == "test_scheme" and info["status"] == "draft"
	assert info["source"] == DATA_REPO and info["ref"] == "dev" and info["commit"] == COMMIT
	assert info["known_issues"] == "B.1,A.1"
	assert info["format"] == "2.3" and info["built_with"] == __version__
	assert os.path.exists(provenance_path(db))
	level, text = scheme_notice(db)
	assert level == logging.WARNING
	assert "draft" in text and "B.1" in text


def test_a_scheme_without_known_issues_records_an_empty_value(tmp_path):
	db = str(tmp_path / "test_scheme.db")
	build_scheme("test_scheme", db, source=GitHubSchemeSource(ref="master", fetcher=FakeRepo(scheme_repo())))
	info = read_scheme_info(db)
	assert info["known_issues"] == "" and info["status"] == "verified"
	assert scheme_notice(db)[0] == logging.INFO


def test_unknown_commit_is_recorded_as_unknown(tmp_path):
	from CanSNPer2.modules.FetchScheme import NotFound, RAW_URL, SchemeError
	files = scheme_repo()
	def fetcher(url, accept=None, timeout=None):
		if url.startswith("https://api.github.com/"):
			raise SchemeError("rate limit exceeded")
		prefix = RAW_URL.format(repo=DATA_REPO, ref="dev", path="")
		if url.startswith(prefix) and url[len(prefix):] in files:
			return files[url[len(prefix):]].encode()
		raise NotFound(url)
	db = str(tmp_path / "t.db")
	build_scheme("test_scheme", db, source=GitHubSchemeSource(ref="dev", fetcher=fetcher))
	info = read_scheme_info(db)
	assert info["commit"] == "unknown" and info["ref"] == "dev"
	assert "commit unknown" in describe(info)


def test_describe_gives_one_line():
	info = {"scheme": "brucella_suis", "status": "draft", "source": "FOI-Bioinformatics/CanSNPer2-data",
		"ref": "dev", "commit": "abc1234" + "0" * 33, "built": "2026-11-21", "built_with": "2.3.0", "format": "2.3"}
	assert describe(info) == ("brucella_suis (draft), FOI-Bioinformatics/CanSNPer2-data dev@abc1234, "
		"built 2026-11-21 with CanSNPer2 2.3.0")


def test_notice_levels(tmp_path):
	db = helpers.make_db(tmp_path / "n.db", [("A.1", "chr", 3, "G", "T")])
	assert read_scheme_info(db) == {}
	assert scheme_notice(db) == (logging.WARNING, MISSING_TEXT)
	conn = sqlite3.connect(db)
	write_scheme_info(conn, {"scheme": "s", "status": "verified", "format": "2.3"})
	conn.close()
	assert scheme_notice(db)[0] == logging.INFO
	for change in ({"status": "unknown"}, {"status": "verified", "known_issues": "A.1"}):
		conn = sqlite3.connect(db)
		write_scheme_info(conn, change)
		conn.close()
		assert scheme_notice(db)[0] == logging.WARNING
	conn = sqlite3.connect(db)
	conn.execute("DELETE FROM scheme_info WHERE key = 'status'")
	conn.commit()
	conn.close()
	assert scheme_notice(db)[0] == logging.WARNING, "a database that records no status is unverified"


## typing through the command line, with the toy aligner of test_aligners

TOY_MODULE = '''from CanSNPer2.modules import Aligners
from test_aligners import ToyAligner

class CLIToy(ToyAligner):
	def __init__(self, path=""):
		super(CLIToy, self).__init__(0)

Aligners.ALIGNERS["toy"] = CLIToy
'''

TYPE_SCRIPT = '''import sys
import toy_cli
from CanSNPer2.CanSNPerTree import main
sys.argv[0] = "CanSNPer2"
main()
'''


def toy_database(tmp_path, status=None):
	db = str(tmp_path / "toy.db")
	helpers.make_db(db, [("A.1", "chr", 3, "G", "T"), ("A.2", "chr", 6, "C", "A")], genome="REF_A")
	conn = sqlite3.connect(db)
	conn.execute("INSERT INTO nodes (id, name) VALUES (0, 'root')")
	conn.executemany("INSERT INTO tree (parent, child, rank_i) VALUES (?, ?, 1)", [(0, 0), (0, 1), (1, 2)])
	conn.commit()
	if status:
		write_scheme_info(conn, {"scheme": "toy", "status": status, "format": "2.3"})
	conn.close()
	refdir = tmp_path / "refs"
	refdir.mkdir()
	(refdir / "REF_A.fna").write_text(">chr\nACGTACGT\n")
	(tmp_path / "sample.fna").write_text(">q\nACTTACGT\n")
	(tmp_path / "toy_cli.py").write_text(TOY_MODULE)
	return db


def type_toy(tmp_path, db):
	env = dict(os.environ)
	env["PYTHONPATH"] = os.pathsep.join([ROOT, TESTS, str(tmp_path), env.get("PYTHONPATH", "")])
	return subprocess.run([sys.executable, "-c", TYPE_SCRIPT, "--database", db, "--aligner", "toy",
		"--refdir", "refs", "--tmpdir", str(tmp_path / "tmp"), "--logdir", "logs", "-o", "out",
		"--min_required_hits", "1", "sample.fna"], capture_output=True, text=True, cwd=str(tmp_path), env=env)


def test_typing_warns_once_for_a_draft_database(tmp_path):
	db = toy_database(tmp_path, status="draft")
	result = type_toy(tmp_path, db)
	assert result.returncode == 0, result.stdout + result.stderr
	assert os.path.exists(tmp_path / "out" / "sample_snps.txt")
	notices = [line for line in result.stderr.splitlines() if "toy (draft)" in line]
	assert len(notices) == 1, result.stderr
	assert "WARN" in notices[0]


def test_a_database_without_scheme_info_is_typed_with_one_warning(tmp_path):
	db = toy_database(tmp_path)
	result = type_toy(tmp_path, db)
	assert result.returncode == 0, result.stdout + result.stderr
	assert os.path.exists(tmp_path / "out" / "sample_snps.txt")
	assert result.stderr.count(MISSING_TEXT) == 1, result.stderr


def test_export_logs_the_notice_and_writes_no_table_rows(tmp_path):
	db = toy_database(tmp_path, status="draft")
	result = run(tmp_path, "CanSNPer2.SNPDatabase", ["--database", db, "--export", "-o", "out", "--logs", "logs"])
	assert result.returncode == 0, result.stderr
	assert result.stderr.count("toy (draft)") == 1, result.stderr
	for name in os.listdir(tmp_path / "out"):
		text = (tmp_path / "out" / name).read_text()
		assert "built_with" not in text and "draft" not in text


def test_download_logs_the_notice(tmp_path):
	db = toy_database(tmp_path, status="draft")
	## without reference genomes in the database nothing is fetched from NCBI
	conn = sqlite3.connect(db)
	conn.execute("DELETE FROM snp_references")
	conn.commit()
	conn.close()
	result = run(tmp_path, "CanSNPer2.DownloadGenomes", ["--database", db, "-o", "refs", "--logs", "logs"])
	assert result.returncode == 0, result.stderr
	assert result.stderr.count("toy (draft)") == 1, result.stderr


def test_a_path_with_uri_characters_is_read(tmp_path):
	folder = tmp_path / "a ?b#c%"
	folder.mkdir()
	db = helpers.make_db(folder / "n.db", [("A.1", "chr", 3, "G", "T")])
	conn = sqlite3.connect(db)
	write_scheme_info(conn, {"scheme": "s"})
	conn.close()
	assert read_scheme_info(db) == {"scheme": "s"}
