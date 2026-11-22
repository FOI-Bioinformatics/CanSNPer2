"""CanSNPer2-database: mode flags, messages, export, convert and log folder (low severity fixes)."""
import logging
import os
import sqlite3
import stat
import subprocess
import sys

import pytest

import helpers
from CanSNPer2.modules.InitiateDatabase import CanSNPDatabase
from CanSNPer2.modules import ExportDatabase
from CanSNPer2.modules.ExportDatabase import ExportError, export_source
from CanSNPer2.modules.ConvertPositions import ConversionError, convert_database, convert_table
from CanSNPer2.modules.LogSetup import setup_logging

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TREE = "R\nR\tA.1\nR\tA.1\tA.2\nR\tB.1\n"
ROWS = ["R\ts1\tP\tREF1\tCP000001.1\t10\tA\tG",
	"A.1\ts1\tP\tREF1\tCP000001.1\t20\tT\tC",
	"A.2\ts1\tP\tREF1\tCP000001.1\t30\tG\tA",
	"B.1\ts1\tP\tREF1\tCP000001.1\t40\tC\tA"]


def run(tmp_path, *args, logs=True, cwd=None):
	env = dict(os.environ, PYTHONPATH=ROOT)
	extra = ["--logs", str(tmp_path / "logs")] if logs else []
	return subprocess.run([sys.executable, "-m", "CanSNPer2.SNPDatabase"] + list(args) + extra,
		capture_output=True, text=True, cwd=str(cwd or tmp_path), env=env)


def out(result):
	return result.stdout + result.stderr


@pytest.fixture
def scheme(tmp_path):
	files = helpers.write_sources(tmp_path / "src", ROWS, tree=TREE)
	db = str(tmp_path / "scheme.db")
	CanSNPDatabase(database=db, create=True, tree=files["tree"], annotation=files["snp"],
		references=files["references"], source_type="CanSNPer", verbose=False)
	return db, files


# 10
def test_create_into_missing_folder(tmp_path, scheme):
	_, files = scheme
	r = run(tmp_path, "--create", "-db", str(tmp_path / "nodir" / "x.db"), "--tree", files["tree"],
		"--annotation", files["snp"], "--references", files["references"])
	assert r.returncode == 1
	assert "folder" in out(r) and "does not exist" in out(r)
	assert "connet" not in out(r) and "Traceback" not in out(r)


def test_no_typo_left():
	for dirpath, _, names in os.walk(os.path.join(ROOT, "CanSNPer2")):
		for name in names:
			if name.endswith(".py"):
				with open(os.path.join(dirpath, name)) as f:
					assert "connet" not in f.read()


# 12
@pytest.mark.parametrize("flags,named", [
	(["--export", "--remove", "--parent", "A.1"], ["--export", "--remove"]),
	(["--export", "--create"], ["--export", "--create"]),
	(["--info", "--export"], ["--info", "--export"]),
	(["--export", "--convert"], ["--export", "--convert"]),
	(["--import", "--create"], ["--import", "--create"]),
	(["--mod_file", "x.txt", "--export"], ["--mod_file", "--export"]),
])
def test_conflicting_modes_refused(tmp_path, scheme, flags, named):
	db, _ = scheme
	before = open(db, "rb").read()
	r = run(tmp_path, "-db", db, *flags)
	assert r.returncode == 1
	for flag in named:
		assert flag in out(r)
	assert "cannot be combined" in out(r)
	assert open(db, "rb").read() == before
	assert not os.path.exists(str(tmp_path / "scheme_tree.txt"))


def test_single_modes_still_work(tmp_path, scheme):
	db, _ = scheme
	assert run(tmp_path, "-db", db, "--info").returncode == 0
	r = run(tmp_path, "-db", db, "--export", "--export_format", "newick", "-o", str(tmp_path / "o"))
	assert r.returncode == 0, out(r)


# 13
def test_message_for_database_alone(tmp_path, scheme):
	db, _ = scheme
	r = run(tmp_path, "-db", db)
	assert r.returncode == 1
	assert "nothing to do" in out(r).lower() and "--export" in out(r)


def test_message_for_replace_alone(tmp_path, scheme):
	db, _ = scheme
	r = run(tmp_path, "-db", db, "--replace")
	assert r.returncode == 1
	assert "--replace" in out(r) and "--mod_file" in out(r) and "are empty" not in out(r)


def test_message_for_export_format_without_export(tmp_path, scheme):
	db, _ = scheme
	r = run(tmp_path, "-db", db, "--export_format", "newick")
	assert r.returncode == 1
	assert "--export_format" in out(r) and "--export" in out(r)


def test_nonpositive_joined_position(tmp_path):
	ref = tmp_path / "REF.fna"
	ref.write_text(">c\nACGTACGT\n")
	src = tmp_path / "old.txt"
	src.write_text("snp_id\tstrain\treference\tgenome\tposition\tderived_base\tancestral_base\n"
		"A\tREF\tt\tREF\t0\tT\tA\n")
	with pytest.raises(ConversionError) as e:
		convert_table(str(src), str(tmp_path / "o.txt"), {"REF": str(ref)})
	assert "not a positive" in str(e.value) and "beyond" not in str(e.value)


# 14
def test_database_is_a_directory(tmp_path):
	(tmp_path / "d.db").mkdir()
	r = run(tmp_path, "-db", str(tmp_path / "d.db"), "--info")
	assert r.returncode == 1
	assert "directory" in out(r) and "Traceback" not in out(r)


# 15
def test_export_refuses_to_overwrite(tmp_path, scheme):
	db, _ = scheme
	o = tmp_path / "o"
	assert run(tmp_path, "-db", db, "--export", "-o", str(o)).returncode == 0
	tree = o / "scheme_tree.txt"
	tree.write_text("mine\n")
	r = run(tmp_path, "-db", db, "--export", "-o", str(o))
	assert r.returncode == 1 and "already exists" in out(r)
	assert tree.read_text() == "mine\n"
	r = run(tmp_path, "-db", db, "--export", "--export_format", "newick", "-o", str(o))
	assert r.returncode == 0
	r = run(tmp_path, "-db", db, "--export", "--export_format", "newick", "-o", str(o))
	assert r.returncode == 1 and "already exists" in out(r)


def test_export_leaves_nothing_when_a_write_fails(tmp_path, scheme, monkeypatch):
	db, _ = scheme
	o = tmp_path / "o"
	o.mkdir()
	real = ExportDatabase.os.replace
	calls = []
	def failing(a, b):
		calls.append(b)
		raise OSError("disk full")
	monkeypatch.setattr(ExportDatabase.os, "replace", failing)
	with pytest.raises(ExportError):
		export_source(db, str(o))
	assert os.listdir(str(o)) == []


def test_export_snp_failure_leaves_no_tree(tmp_path, scheme, monkeypatch):
	db, _ = scheme
	o = tmp_path / "o"
	real_print = print
	monkeypatch.setattr(ExportDatabase, "SNP_COLUMNS", None)  ## printing *None fails after the tree is written
	with pytest.raises(Exception):
		export_source(db, str(o))
	assert not os.path.exists(str(o)) or os.listdir(str(o)) == []


@pytest.mark.skipif(os.geteuid() == 0, reason="root can write anywhere")
def test_export_to_unwritable_folder(tmp_path, scheme):
	db, _ = scheme
	o = tmp_path / "ro"
	o.mkdir()
	o.chmod(stat.S_IRUSR | stat.S_IXUSR)
	try:
		r = run(tmp_path, "-db", db, "--export", "-o", str(o))
		assert r.returncode == 1 and "Traceback" not in out(r)
		r = run(tmp_path, "-db", db, "--export", "--export_format", "newick", "-o", str(o))
		assert r.returncode == 1 and "Traceback" not in out(r)
	finally:
		o.chmod(stat.S_IRWXU)


# 17
def test_convert_current_format_stops_before_download(tmp_path, scheme):
	db, files = scheme
	r = run(tmp_path, "--convert", "--database", db, "--refdir", str(tmp_path / "refs"), "-o", str(tmp_path / "o"))
	assert r.returncode == 1 and "current format" in out(r)
	assert not os.path.exists(str(tmp_path / "refs"))
	r = run(tmp_path, "--convert", "--annotation", files["snp"], "--refdir", str(tmp_path / "refs"),
		"-o", str(tmp_path / "o"))
	assert r.returncode == 1 and "current format" in out(r)
	assert not os.path.exists(str(tmp_path / "refs"))


# 18
def test_convert_writes_upper_case_bases(tmp_path):
	ref = tmp_path / "REF.fna"
	ref.write_text(">c\nACGTACGT\n")
	src = tmp_path / "old.txt"
	src.write_text("snp_id\tstrain\treference\tgenome\tposition\tderived_base\tancestral_base\n"
		"A\tREF\tt\tREF\t2\tg\tc\n")
	dst = tmp_path / "o.txt"
	convert_table(str(src), str(dst), {"REF": str(ref)})
	assert dst.read_text().splitlines()[1].split("\t")[-2:] == ["G", "C"]
	db = helpers.make_db(tmp_path / "old.db", [("A", None, 2, "c", "g")], genome="REF", joined=True)
	target = tmp_path / "new.db"
	convert_database(db, str(target), {"REF": str(ref)})
	conn = sqlite3.connect(str(target))
	assert conn.execute("SELECT ancestral_base, derived_base FROM snp_annotation").fetchall() == [("C", "G")]


# 19
def test_log_folder_not_writable(tmp_path):
	blocked = tmp_path / "blocked"
	blocked.write_text("a file, so the log folder cannot be made")
	setup_logging(str(blocked / "logs"), "tag")
	logging.getLogger("x").warning("still logged")


@pytest.mark.skipif(os.geteuid() == 0, reason="root can write anywhere")
def test_info_in_read_only_directory(tmp_path, scheme):
	db, _ = scheme
	ro = tmp_path / "ro"
	ro.mkdir()
	ro.chmod(stat.S_IRUSR | stat.S_IXUSR)
	try:
		r = run(tmp_path, "-db", db, "--info", logs=False, cwd=ro)
		assert r.returncode == 0, out(r)
		assert out(r).count("log") >= 1 and "Traceback" not in out(r)
	finally:
		ro.chmod(stat.S_IRWXU)
