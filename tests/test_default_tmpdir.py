"""The default --tmpdir. /tmp (or TMPDIR) can be read-only, full, or hold a CanSNPer2 folder
owned by another user; the run then stopped before typing anything. Without --tmpdir the next
usable location is taken. A --tmpdir that was given is used as it is, and a bad one is an error."""
import os
import subprocess
import sys

import pytest

from CanSNPer2.CanSNPerTree import default_tmpdir_candidates
from test_scheme_info import ROOT, TESTS, TYPE_SCRIPT, toy_database

pytestmark = pytest.mark.skipif(hasattr(os, "geteuid") and os.geteuid() == 0, reason="root can write anywhere")


def type_with(tmp_path, db, system_tmp, *extra):
	env = dict(os.environ)
	env["PYTHONPATH"] = os.pathsep.join([ROOT, TESTS, str(tmp_path), env.get("PYTHONPATH", "")])
	env["TMPDIR"] = str(system_tmp)
	return subprocess.run([sys.executable, "-c", TYPE_SCRIPT, "--database", db, "--aligner", "toy",
		"--refdir", "refs", "--logdir", "logs", "-o", "out", "--min_required_hits", "1", "sample.fna"] + list(extra),
		capture_output=True, text=True, cwd=str(tmp_path), env=env)


def test_candidates_follow_tmpdir_then_the_output_folder(tmp_path, monkeypatch):
	(tmp_path / "sys").mkdir()
	monkeypatch.setenv("TMPDIR", str(tmp_path / "sys"))
	import tempfile
	tempfile.tempdir = None
	assert default_tmpdir_candidates("out") == [str(tmp_path / "sys" / "CanSNPer2"), os.path.join("out", ".cansnper2_tmp")]


def test_an_unwritable_default_folder_falls_back_to_the_output_folder(tmp_path):
	"""the system folder is writable but its CanSNPer2 folder is not, as when another user made it"""
	db = toy_database(tmp_path)
	system = tmp_path / "sys"
	locked = system / "CanSNPer2"
	locked.mkdir(parents=True)
	locked.chmod(0o500)
	try:
		result = type_with(tmp_path, db, system)
	finally:
		locked.chmod(0o700)
	assert result.returncode == 0, result.stdout + result.stderr
	assert os.path.exists(tmp_path / "out" / "samples" / "sample_snps.txt")
	assert "trying another location" in result.stderr


def test_a_given_tmpdir_that_is_not_writable_is_still_an_error(tmp_path):
	db = toy_database(tmp_path)
	locked = tmp_path / "locked"
	locked.mkdir()
	locked.chmod(0o500)
	try:
		result = type_with(tmp_path, db, tmp_path, "--tmpdir", str(locked))
	finally:
		locked.chmod(0o700)
	assert result.returncode != 0
	assert "--tmp-dir" in result.stderr and "cannot be created" in result.stderr
