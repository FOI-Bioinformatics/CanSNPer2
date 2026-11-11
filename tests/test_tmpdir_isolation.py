"""Regression tests for the shared temporary directory: every CanSNPer2 run used to work
directly in <tmpdir>/CanSNPer2 (by default /tmp/CanSNPer2) and cleanup() removed every
file in it, so two runs on the same machine deleted each other's alignments and staged
genomes, and a run whose tmpdir held a subdirectory crashed at the end.

Each aligning run now works in a private run_* directory below the shared one and only
removes that directory; --skip_mauve keeps reading the directory it is given."""
import os

from CanSNPer2.modules.CanSNPer2 import CanSNPer2, RUN_MARKER
from CanSNPer2.CanSNPerTree import resolve_tmpdir


def make_run(tmp_path, **kwargs):
	return CanSNPer2(["query.fna"], workdir=str(tmp_path), tmpdir="shared", logdir="logs", outdir="out", **kwargs)


def test_each_run_gets_its_own_directory(tmp_path):
	first = make_run(tmp_path)
	second = make_run(tmp_path)
	shared = os.path.join(str(tmp_path), "shared")
	assert first.tmpdir != second.tmpdir
	for run in (first, second):
		assert os.path.dirname(run.tmpdir) == shared
		assert os.path.basename(run.tmpdir).startswith("run_")
		assert os.path.isfile(os.path.join(run.tmpdir, RUN_MARKER))


def test_cleanup_leaves_other_runs_and_foreign_files(tmp_path):
	first = make_run(tmp_path)
	second = make_run(tmp_path)
	shared = os.path.join(str(tmp_path), "shared")
	foreign = os.path.join(shared, "not_from_cansnper2.txt")
	with open(foreign, "w") as fh:
		fh.write("keep me")
	other_alignment = os.path.join(second.tmpdir, "REF_query.xmfa")
	with open(other_alignment, "w") as fh:
		fh.write("#FormatVersion Mauve1\n")
	'''a subdirectory in the run directory made the old cleanup fail on os.remove'''
	os.makedirs(os.path.join(first.tmpdir, "subdir"))
	with open(os.path.join(first.tmpdir, "REF_query.xmfa"), "w") as fh:
		fh.write("#FormatVersion Mauve1\n")

	first.cleanup()

	assert not os.path.exists(first.tmpdir)
	assert os.path.isfile(foreign)
	assert os.path.isfile(other_alignment)
	assert os.path.isdir(shared)


def test_skip_mauve_reads_the_given_directory_and_never_removes_it(tmp_path):
	shared = os.path.join(str(tmp_path), "shared")
	os.makedirs(shared)
	alignment = os.path.join(shared, "REF_query.xmfa")
	with open(alignment, "w") as fh:
		fh.write("#FormatVersion Mauve1\n")
	run = make_run(tmp_path, skip_mauve=True)
	assert run.tmpdir == shared
	run.cleanup()
	assert os.path.isfile(alignment)


def test_resolve_tmpdir(tmp_path):
	base = str(tmp_path / "scratch")
	assert resolve_tmpdir(base) == base + "/CanSNPer2"
	assert resolve_tmpdir(base + "/CanSNPer2/") == base + "/CanSNPer2"
	'''a kept run directory is used as it is, so its alignments can be reused with --skip_mauve'''
	kept = tmp_path / "scratch" / "CanSNPer2" / "run_abc123"
	kept.mkdir(parents=True)
	(kept / RUN_MARKER).write_text("")
	assert resolve_tmpdir(str(kept)) == str(kept)
	'''a directory that only looks like a run directory is not trusted without the marker'''
	lookalike = tmp_path / "run_xyz"
	lookalike.mkdir()
	assert resolve_tmpdir(str(lookalike)) == str(lookalike) + "/CanSNPer2"
