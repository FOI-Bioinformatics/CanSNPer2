"""The steps a typing run is set up in (CanSNPerTree.run_command split into functions, #103),
each tested on its own: the log level, the folders of the run read against --workdir, the
default and the given --tmp-dir, and whether a run typed anything. The whole command is
exercised through `cansnper2 type` in test_default_tmpdir.py and test_typing_cli_checks.py."""

import logging
import os
import tempfile
from types import SimpleNamespace

import pytest

from CanSNPer2.CanSNPerTree import (
    in_workdir,
    log_level,
    resolve_temp_folder,
    resolve_tmpdir,
    typed_nothing,
)


def options(**changes):
    fields = dict(
        debug=None,
        verbose=None,
        supress=logging.WARNING,
        workdir=".",
        outdir="results",
        tmpdir=None,
        skip_alignment=False,
    )
    fields.update(changes)
    return SimpleNamespace(**fields)


def test_log_level_follows_debug_then_verbose_then_supress():
    assert log_level(options()) == logging.WARNING
    assert log_level(options(supress=logging.ERROR)) == logging.ERROR
    assert log_level(options(verbose=logging.INFO)) == logging.INFO
    assert log_level(options(verbose=logging.INFO, debug=logging.DEBUG)) == logging.DEBUG


def test_a_folder_of_the_run_is_read_against_the_workdir(tmp_path):
    args = options(workdir=str(tmp_path))
    assert in_workdir(args, "logs/") == os.path.join(str(tmp_path), "logs")
    assert in_workdir(args, "/elsewhere/logs") == "/elsewhere/logs", "an absolute path stands"


def test_the_default_temp_folder_is_the_first_that_can_be_written_to(tmp_path, monkeypatch):
    (tmp_path / "systmp").mkdir()
    monkeypatch.setattr(tempfile, "tempdir", str(tmp_path / "systmp"))
    args = options(workdir=str(tmp_path))
    resolve_temp_folder(args, logging.getLogger(__name__))
    assert args.tmpdir == os.path.join(str(tmp_path / "systmp"), "CanSNPer2")
    assert os.path.isdir(args.tmpdir)


def test_the_output_folder_holds_the_temp_folder_when_the_system_one_cannot_be_used(
    tmp_path, monkeypatch, caplog
):
    blocked = tmp_path / "blocked"
    blocked.mkdir(mode=0o500)
    monkeypatch.setattr(tempfile, "tempdir", str(blocked))
    args = options(workdir=str(tmp_path))
    with caplog.at_level(logging.WARNING):
        resolve_temp_folder(args, logging.getLogger(__name__))
    assert args.tmpdir == os.path.join(str(tmp_path), "results", ".cansnper2_tmp", "CanSNPer2")
    assert "trying another location" in caplog.text


def test_a_given_temp_folder_is_read_against_the_workdir_and_gets_the_run_folder(tmp_path):
    args = options(workdir=str(tmp_path), tmpdir="scratch")
    resolve_temp_folder(args, logging.getLogger(__name__))
    assert args.tmpdir == os.path.join(str(tmp_path), "scratch", "CanSNPer2")
    assert os.path.isdir(args.tmpdir)
    assert resolve_tmpdir(args.tmpdir) == args.tmpdir, "the run folder is not added twice"


def test_skip_alignment_does_not_create_the_temp_folder_it_reads(tmp_path):
    args = options(workdir=str(tmp_path), tmpdir="gone", skip_alignment=True)
    with pytest.raises(SystemExit) as stopped:
        resolve_temp_folder(args, logging.getLogger(__name__))
    assert stopped.value.code == 1
    assert not os.path.exists(os.path.join(str(tmp_path), "gone"))


def test_a_run_that_typed_nothing_is_a_failure():
    assert typed_nothing(SimpleNamespace(typed={}, reused=set(), not_typed={"x": "NA"}, skipped={}))
    assert typed_nothing(SimpleNamespace(typed={}, reused=set(), not_typed={}, skipped={"x": "no"}))
    assert not typed_nothing(SimpleNamespace(typed={"x": "B.1"}, reused=set(), not_typed={}, skipped={}))
    assert not typed_nothing(SimpleNamespace(typed={}, reused={"x"}, not_typed={}, skipped={})), (
        "a sample whose earlier result was used again counts as typed"
    )
    assert not typed_nothing(SimpleNamespace(typed={}, reused=set(), not_typed={}, skipped={})), (
        "a run without queries is not a failed run"
    )
