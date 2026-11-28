"""The steps of one round of alignments (CanSNPer2.run_alignments split into methods, #103),
each tested on its own: starting an alignment, collecting the ones that finished, stopping the
ones that run too long, and what the round returns. The round as a whole is tested in
test_parallel_samples.py (--threads) and test_download_inputs_timeout.py (--alignment-timeout)."""

import sys
import time

import pytest

from CanSNPer2.modules.CanSNPer2 import AlignmentError
from test_typing_inputs import project

EXIT = "import sys; sys.exit(int(sys.argv[1]))"
SLEEP = "import time; time.sleep(30)"


def runner(tmp_path, **kwargs):
    obj, _db = project(tmp_path, [], **kwargs)
    return obj


def started(obj, tmp_path, code="0", name="a"):
    command = [sys.executable, "-c", EXIT, code]
    return obj.start_alignment(command, str(tmp_path / (name + ".log")), ("REF", name))


def finished(obj, running):
    """collect the alignments once they have all ended"""
    for job in running:
        job.process.wait()
    return obj.collect_alignments(running)


def test_an_alignment_writes_its_output_to_its_log_file(tmp_path):
    obj = runner(tmp_path, aligner="progressiveMauve")
    job = obj.start_alignment(
        [sys.executable, "-c", "import sys; print('to the log'); sys.stderr.write('and this\\n')"],
        str(tmp_path / "a.log"),
        ("REF", "a"),
    )
    assert job.names == ("REF", "a") and job.logpath == str(tmp_path / "a.log")
    assert job.started <= time.time()
    code, failed, had_finished = finished(obj, [job])
    assert (code, failed, had_finished) == (0, 0, True)
    assert sorted(open(str(tmp_path / "a.log")).read().splitlines()) == ["and this", "to the log"], (
        "stderr goes to the log file too"
    )


def test_collecting_takes_the_alignments_that_ended_out_of_the_round(tmp_path):
    obj = runner(tmp_path, aligner="progressiveMauve")
    running = [
        started(obj, tmp_path, name="a"),
        obj.start_alignment([sys.executable, "-c", SLEEP], str(tmp_path / "b.log"), ("REF", "b")),
    ]
    running[0].process.wait()
    code, failed, had_finished = obj.collect_alignments(running)
    assert (code, failed, had_finished) == (0, 0, True)
    assert [job.names[1] for job in running] == ["b"], "the one still running stays"
    obj.stop_alignments(running)
    assert running[0].process.poll() is not None


def test_collecting_says_nothing_happened_while_an_alignment_runs(tmp_path):
    obj = runner(tmp_path, aligner="progressiveMauve")
    running = [obj.start_alignment([sys.executable, "-c", SLEEP], str(tmp_path / "b.log"), ("REF", "b"))]
    assert obj.collect_alignments(running) == (0, 0, False)
    assert len(running) == 1
    obj.stop_alignments(running)


def test_a_failed_alignment_is_counted(tmp_path, caplog):
    obj = runner(tmp_path, aligner="progressiveMauve")
    code, failed, _ = finished(obj, [started(obj, tmp_path, code="4")])
    assert (code, failed) == (0, 1)
    assert "non zero exitcode: 4" in caplog.text


def test_an_exit_code_the_aligner_retries_is_returned(tmp_path, caplog):
    obj = runner(tmp_path, aligner="progressiveMauve")
    code, failed, _ = finished(obj, [started(obj, tmp_path, code="11")])
    assert (code, failed) == (11, 0), "progressiveMauve exit code 11 asks for one retry"
    assert "exitcode: 11" in caplog.text


def test_an_alignment_that_runs_too_long_is_stopped(tmp_path, caplog):
    obj = runner(tmp_path, aligner="progressiveMauve", alignment_timeout=0.02)  ## 1.2 seconds
    running = [obj.start_alignment([sys.executable, "-c", SLEEP], str(tmp_path / "b.log"), ("REF", "b"))]
    assert not obj.timed_out(running)
    running[0].started = time.time() - 2
    assert obj.timed_out(running)
    obj.stop_timed_out_alignments(running)
    assert running[0].process.poll() is not None
    assert "did not finish within 0.02 minutes" in caplog.text


def test_no_timeout_with_the_option_set_to_zero(tmp_path):
    obj = runner(tmp_path, aligner="progressiveMauve", alignment_timeout=0)
    job = started(obj, tmp_path)
    job.started = time.time() - 10**6
    assert not obj.timed_out([job])
    finished(obj, [job])


def test_the_round_returns_the_retry_code_before_anything_else(tmp_path):
    obj = runner(tmp_path, aligner="progressiveMauve")
    assert obj.alignment_round_result(11, 3) == 11


def test_the_round_raises_when_an_alignment_failed(tmp_path):
    obj = runner(tmp_path, aligner="progressiveMauve")
    assert obj.alignment_round_result(0, 0) == 0
    with pytest.raises(AlignmentError, match="2 progressiveMauve processes did not run correctly"):
        obj.alignment_round_result(0, 2)


def test_keep_going_turns_a_failure_into_one_without_stopping_the_run(tmp_path):
    obj = runner(tmp_path, aligner="progressiveMauve", keep_going=True)
    assert obj.alignment_round_result(0, 2) == 1
    assert obj.alignment_round_result(0, 0) == 0, (
        "a round where every process succeeded still reports success (issue #126)"
    )
