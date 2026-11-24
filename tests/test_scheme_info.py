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
	paths = helpers.write_sources(tmp_path / "src", ["A.1\t1\tBirdsell2014\tREF1\tCP000001.1\t100\tT\tC\t\tyes",
		"A.1\t2\tBirdsell2014\tREF1\tCP000001.1\t150\tA\tG\t480\tno",
		"B.1\t1\tBirdsell2014\tREF1\tCP000001.1\t200\tG\tA\t\tyes"])
	return run(tmp_path, "CanSNPer2.SNPDatabase", ["--create", "--database", db, "--tree", paths["tree"],
		"--annotation", paths["snp"], "--references", paths["references"], "--logs", "logs"] + list(extra)), paths


def test_create_writes_format_and_version(tmp_path):
	db = str(tmp_path / "s.db")
	result, _paths = create(tmp_path, db, "--scheme_name", "my_scheme", "--scheme_status", "draft",
		"--scheme_source", "my-folder")
	assert result.returncode == 0, result.stderr
	assert read_scheme_info(db) == {"format": "2.4", "built_with": __version__,
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
	assert info["format"] == "2.4" and info["built_with"] == __version__
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
	assert "format\t2.4" in lines and "built_with\t{v}".format(v=__version__) in lines
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
	assert info["format"] == "2.4" and info["built_with"] == __version__
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


def type_toy(tmp_path, db, *extra):
	env = dict(os.environ)
	env["PYTHONPATH"] = os.pathsep.join([ROOT, TESTS, str(tmp_path), env.get("PYTHONPATH", "")])
	return subprocess.run([sys.executable, "-c", TYPE_SCRIPT, "--database", db, "--aligner", "toy",
		"--refdir", "refs", "--tmpdir", str(tmp_path / "tmp"), "--logdir", "logs", "-o", "out",
		"--min_required_hits", "1", "sample.fna"] + list(extra), capture_output=True, text=True, cwd=str(tmp_path), env=env)


def test_typing_warns_once_for_a_draft_database(tmp_path):
	db = toy_database(tmp_path, status="draft")
	result = type_toy(tmp_path, db)
	assert result.returncode == 0, result.stdout + result.stderr
	assert os.path.exists(tmp_path / "out" / "samples" / "sample_snps.txt")
	notices = [line for line in result.stderr.splitlines() if "toy (draft)" in line]
	assert len(notices) == 1, result.stderr
	assert "WARN" in notices[0]


def test_a_database_without_scheme_info_is_typed_with_one_warning(tmp_path):
	db = toy_database(tmp_path)
	result = type_toy(tmp_path, db)
	assert result.returncode == 0, result.stdout + result.stderr
	assert os.path.exists(tmp_path / "out" / "samples" / "sample_snps.txt")
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
	## without reference genomes in the database nothing is fetched from NCBI, and the
	## command ends with exit 1 after the notice
	conn = sqlite3.connect(db)
	conn.execute("DELETE FROM snp_references")
	conn.commit()
	conn.close()
	result = run(tmp_path, "CanSNPer2.DownloadGenomes", ["--database", db, "-o", "refs", "--logs", "logs"])
	assert result.returncode == 1, result.stderr
	assert "nothing to download" in result.stderr
	assert result.stderr.count("toy (draft)") == 1, result.stderr


def test_a_path_with_uri_characters_is_read(tmp_path):
	folder = tmp_path / "a ?b#c%"
	folder.mkdir()
	db = helpers.make_db(folder / "n.db", [("A.1", "chr", 3, "G", "T")])
	conn = sqlite3.connect(db)
	write_scheme_info(conn, {"scheme": "s"})
	conn.close()
	assert read_scheme_info(db) == {"scheme": "s"}


def test_info_of_a_file_that_is_not_a_database(tmp_path):
	text = tmp_path / "notes.db"
	text.write_text("this is not a database\n" * 20)
	result = run(tmp_path, "CanSNPer2.SNPDatabase", ["--database", str(text), "--info", "--logs", "logs"])
	assert result.returncode == 1
	assert "is not an sqlite database" in result.stderr
	assert "Traceback" not in result.stderr


def test_known_issues_that_are_not_a_list_are_not_recorded(tmp_path, caplog):
	files = scheme_repo()
	files["schemes.toml"] += '[test_scheme]\nknown_issues = 5\n'
	db = str(tmp_path / "t.db")
	build_scheme("test_scheme", db, source=GitHubSchemeSource(ref="master", fetcher=FakeRepo(files)))
	assert read_scheme_info(db)["known_issues"] == ""
	assert "is not a list of snp_ids" in caplog.text


def test_a_local_checkout_with_uncommitted_scheme_files_is_marked(tmp_path):
	import shutil
	from test_fetch_scheme import write_tree
	from CanSNPer2.modules.FetchScheme import LocalSchemeSource
	if not shutil.which("git"):
		pytest.skip("git is not installed")
	data = write_tree(tmp_path / "data", scheme_repo())
	git = lambda *a: subprocess.run(["git", "-C", str(data)] + list(a), check=True, capture_output=True, text=True)
	git("init", "-q")
	git("add", ".")
	git("-c", "user.name=t", "-c", "user.email=t@t", "commit", "-q", "-m", "schemes")
	head = git("rev-parse", "HEAD").stdout.strip()
	build_scheme("test_scheme", str(tmp_path / "a.db"), source=LocalSchemeSource(str(data)))
	assert read_scheme_info(str(tmp_path / "a.db"))["commit"] == head
	(data / "schemes/test_scheme/test_scheme_tree.txt").write_text("R\nR\tA.1\nR\tA.1\tB.1\nR\tA.1\tC.1\n")
	build_scheme("test_scheme", str(tmp_path / "b.db"), source=LocalSchemeSource(str(data)))
	info = read_scheme_info(str(tmp_path / "b.db"))
	assert info["commit"] == head + "+uncommitted"
	assert "@{sha}+uncommitted".format(sha=head[:7]) in describe(info)


DOWNLOAD_SCRIPT = '''import sys
from CanSNPer2 import DownloadGenomes
DownloadGenomes.DownloadGenomes.run = lambda self: 0  ## no reference is fetched from NCBI
sys.argv[0] = "CanSNPer2-download"
DownloadGenomes.main()
'''


def test_a_draft_scheme_download_reports_the_draft_once(tmp_path):
	from test_fetch_scheme import write_tree
	data = write_tree(tmp_path / "data", scheme_repo(status="draft"))
	env = dict(os.environ)
	env["PYTHONPATH"] = os.pathsep.join([ROOT, TESTS, env.get("PYTHONPATH", "")])
	result = subprocess.run([sys.executable, "-c", DOWNLOAD_SCRIPT, "--scheme", "test_scheme", "--data_dir", str(data),
		"--logs", "logs"], capture_output=True, text=True, cwd=str(tmp_path), env=env)
	assert result.returncode == 0, result.stderr
	assert result.stderr.count("is a draft") + result.stderr.count("(draft)") == 1, result.stderr


## a database changed after it was built ------------------------------------------------

def verified_toy(tmp_path, commit="5b24118aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"):
	'''the toy database recorded as a verified scheme of a repository'''
	db = toy_database(tmp_path)
	conn = sqlite3.connect(db)
	## the modification commands read the rank table, which make_db does not write
	conn.execute("CREATE TABLE rank (rank_i integer PRIMARY KEY, rank VARCHAR(15))")
	conn.execute("INSERT INTO rank (rank_i, rank) VALUES (1, 'family')")
	conn.commit()
	info = {"scheme": "toy", "status": "verified", "source": "owner/data", "format": "2.3",
		"built_with": "2.2.9", "built": "2026-01-01"}
	if commit is not None:
		info["commit"] = commit
	write_scheme_info(conn, info)
	conn.close()
	return db


def test_a_removal_marks_a_verified_database_as_modified(tmp_path):
	db = verified_toy(tmp_path)
	result = run(tmp_path, "CanSNPer2.SNPDatabase", ["--database", db, "--remove", "--parent", "A.2", "--logs", "logs"])
	assert result.returncode == 0, result.stderr
	result = run(tmp_path, "CanSNPer2.SNPDatabase", ["--database", db, "--info", "--logs", "logs"])
	assert result.returncode == 0, result.stderr
	lines = result.stdout.splitlines()
	assert "status\tunknown" in lines, lines
	assert "commit\t5b24118aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa+modified" in lines, lines
	assert "modified\t{d}".format(d=datetime.date.today().isoformat()) in lines, lines
	assert "built_with\t{v}".format(v=__version__) in lines, lines
	assert "source\towner/data" in lines and "built\t2026-01-01" in lines, lines
	result = type_toy(tmp_path, db)
	assert result.returncode == 0, result.stdout + result.stderr
	notices = [line for line in result.stderr.splitlines() if "toy (unknown)" in line]
	assert len(notices) == 1 and "WARN" in notices[0], result.stderr
	assert "@5b24118+modified" in notices[0] and "modified {d}".format(d=datetime.date.today().isoformat()) in notices[0]


def test_a_mod_file_marks_a_database_without_a_commit_as_modified(tmp_path):
	db = verified_toy(tmp_path, commit=None)
	(tmp_path / "mod.txt").write_text("child\tparent\nC.1\tA.1\n")
	result = run(tmp_path, "CanSNPer2.SNPDatabase", ["--database", db, "--mod_file", "mod.txt", "--parent", "A.1",
		"--logs", "logs"])
	assert result.returncode == 0, result.stderr
	info = read_scheme_info(db)
	assert info["status"] == "unknown" and info["commit"] == "modified", info
	assert info["modified"] == datetime.date.today().isoformat() and info["built_with"] == __version__
	assert "toy (unknown), owner/data@modified" in describe(info), describe(info)


def test_a_failed_modification_leaves_the_scheme_info(tmp_path):
	db = verified_toy(tmp_path)
	before = read_scheme_info(db)
	result = run(tmp_path, "CanSNPer2.SNPDatabase", ["--database", db, "--remove", "--parent", "NOT_A_NODE", "--logs", "logs"])
	assert result.returncode == 1
	assert read_scheme_info(db) == before


def test_info_prints_none_for_an_empty_value(tmp_path):
	db = verified_toy(tmp_path)
	conn = sqlite3.connect(db)
	write_scheme_info(conn, {"known_issues": ""})
	conn.close()
	result = run(tmp_path, "CanSNPer2.SNPDatabase", ["--database", db, "--info", "--logs", "logs"])
	assert result.returncode == 0, result.stderr
	assert "known_issues\tnone" in result.stdout.splitlines(), result.stdout


## a database genome without its reference file -----------------------------------------

def two_genome_toy(tmp_path):
	'''the toy database with A.2 moved to a second genome REF_B, which has no file in refs'''
	db = toy_database(tmp_path, status="verified")
	conn = sqlite3.connect(db)
	conn.execute("INSERT INTO snp_references (id, genome, strain) VALUES (?, 'REF_B', 'REF_B')", (helpers.GENOME_I + 1,))
	conn.execute("UPDATE snp_annotation SET genome_i=? WHERE snp_id='A.2'", (helpers.GENOME_I + 1,))
	conn.commit()
	conn.close()
	return db


def test_a_missing_reference_file_stops_typing(tmp_path):
	db = two_genome_toy(tmp_path)
	result = type_toy(tmp_path, db)
	assert result.returncode == 1, result.stdout + result.stderr
	assert "REF_B" in result.stderr and "--allow_missing_references" in result.stderr, result.stderr
	assert "Traceback" not in result.stderr


def test_allow_missing_references_types_with_one_warning(tmp_path):
	db = two_genome_toy(tmp_path)
	result = type_toy(tmp_path, db, "--allow_missing_references")
	assert result.returncode == 0, result.stdout + result.stderr
	assert os.path.exists(tmp_path / "out" / "samples" / "sample_snps.txt")
	warnings = [line for line in result.stderr.splitlines() if "REF_B" in line]
	assert len(warnings) == 1 and "WARN" in warnings[0], result.stderr
	assert "not compared" in warnings[0], warnings


def test_allow_missing_references_still_stops_when_no_genome_has_a_file(tmp_path):
	db = two_genome_toy(tmp_path)
	os.rename(str(tmp_path / "refs" / "REF_A.fna"), str(tmp_path / "refs" / "OTHER.fna"))
	result = type_toy(tmp_path, db, "--allow_missing_references")
	assert result.returncode == 1, result.stdout + result.stderr
	assert "REF_A, REF_B" in result.stderr and "cansnper2 ref get" in result.stderr, result.stderr
	assert "Traceback" not in result.stderr
