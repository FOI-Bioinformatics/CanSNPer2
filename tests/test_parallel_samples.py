"""--jobs (#100): samples typed at the same time give the same results, in the same order, as
samples typed one after the other; and --threads limits the alignments of an assembly that run
at the same time."""

import json
import os
import sys
import time

import pytest

from CanSNPer2.modules.CanSNPer2 import CanSNPer2Error
from test_result_files import comments
from test_typing_inputs import fake_alignment, project, write

STATES = {
    "a": (["A.1", "B.1"], {"A.1": 1, "B.1": 1}),
    "b": (["A.1", "B.2"], {"A.1": 1, "B.2": 1}),
    "c": (["C.1"], {"C.1": 1, "A.1": 2}),
    "d": ([], {"A.1": 2}),
}


def typed_samples(tmp_path, jobs, names="abcd", **kwargs):
    """samples with given states (no aligner); the states of each sample are looked up by the
    sample name, so that they follow it into the worker that types it"""
    queries = [write(tmp_path / "in" / (name + ".fasta")) for name in names]
    obj, db = project(tmp_path, queries, min_required_hits=1, jobs=jobs, summary=True, **kwargs)
    obj.align = lambda q, references=(): ["REF_x.paf"]

    def states(**_kwargs):
        called, sites = STATES[obj.query_name]
        rows = [[name, "REF", "chr1", "1", "chr1", "1", "C", "T", "T", ""] for name in called]
        return sites, rows, list(called)

    obj.find_snps_multiproc = states
    obj.run(database=db)
    return obj


def text_files(outdir, root):
    """the text of the output files, with the folder of the run and the dates left out"""
    found = {}
    for folder, _dirs, names in os.walk(outdir):
        for name in names:
            path = os.path.join(folder, name)
            text = open(path).read().replace(str(root), "<root>")
            found[os.path.relpath(path, outdir)] = "\n".join(
                line
                for line in text.splitlines()
                if "Date" not in line and '"date"' not in line and "Sample." not in line
            )
    return found


def test_jobs_give_the_same_results_as_one_sample_at_a_time(tmp_path):
    one = typed_samples(tmp_path / "one", jobs=1)
    three = typed_samples(tmp_path / "three", jobs=3)
    assert three.typed == one.typed and three.not_typed == one.not_typed
    assert list(three.sample_reports) == list(one.sample_reports) == ["a", "b", "c", "d"]
    a, b = text_files(one.outdir, tmp_path / "one"), text_files(three.outdir, tmp_path / "three")
    assert sorted(a) == sorted(b)
    for name in a:
        if name.endswith(".html"):
            continue  ## the pages embed the run folder and time
        assert a[name] == b[name], name


def test_jobs_keep_the_order_of_the_samples_in_the_summary_and_run_json(tmp_path):
    obj = typed_samples(tmp_path, jobs=4)
    rows = [
        line.split("\t")[0]
        for line in open(os.path.join(obj.outdir, "snps_summary.txt"))
        if not line.startswith("#")
    ]
    assert rows == ["sample", "a", "b", "c", "d"]
    run = json.load(open(os.path.join(obj.outdir, "run.json")))
    assert [s["sample"] for s in run["samples"]] == ["a", "b", "c", "d"]
    assert run["samples"][0]["result"] == "samples/a.json"
    assert comments(os.path.join(obj.results_dir, "b_snps.txt"))[1].startswith("# Final SNP: B.2")


def test_a_failed_sample_with_jobs_and_keep_going_is_reported(tmp_path):
    (tmp_path / "in").mkdir(parents=True)
    (tmp_path / "in" / "x").mkdir()
    queries = [
        write(tmp_path / "in" / "a.fasta"),
        str(tmp_path / "in" / "x"),
        write(tmp_path / "in" / "b.fasta"),
    ]
    obj, db = project(tmp_path, queries, min_required_hits=1, jobs=2, keep_going=True)
    fake_alignment(obj, ["A.1", "B.1"], {"A.1": 1, "B.1": 1}, queries[0])
    obj.run(database=db)
    assert list(obj.sample_reports) == ["a", "x", "b"]
    assert obj.sample_reports["x"]["severity"] == "failed"
    assert "is a directory" in obj.not_typed["x"]
    assert set(obj.typed) == {"a", "b"}


def test_a_failed_sample_with_jobs_stops_the_run_without_keep_going(tmp_path):
    """the states of sample e are not known, so reading them fails while the sample is typed"""
    with pytest.raises(CanSNPer2Error, match="use --keep-going"):
        typed_samples(tmp_path, jobs=2, names="abed")


def test_a_sample_that_fails_while_it_is_typed_is_reported_with_keep_going(tmp_path):
    obj = typed_samples(tmp_path, jobs=2, names="abed", keep_going=True)
    assert list(obj.sample_reports) == ["a", "b", "e", "d"]
    assert obj.sample_reports["e"]["severity"] == "failed" and obj.not_typed["e"].startswith("failed:")


def test_jobs_reuse_earlier_results(tmp_path):
    typed_samples(tmp_path, jobs=2)
    again = typed_samples(tmp_path, jobs=2)
    assert again.reused == {"a", "b", "c", "d"}
    assert again.sample_reports["a"]["severity"] == "reused"


# --threads limits the alignments of one assembly that run at the same time


SLEEPER = "import sys, time; open(sys.argv[1], 'a').write('%f start\\n' % time.time()); time.sleep(0.4); open(sys.argv[1], 'a').write('%f end\\n' % time.time())"


def most_at_once(log):
    events = sorted((float(t), kind) for t, kind in (line.split() for line in open(log)))
    running = most = 0
    for _t, kind in events:
        running += 1 if kind == "start" else -1
        most = max(most, running)
    return most


@pytest.mark.parametrize("threads,expected", [(2, 2), (1, 1), (8, 5)])
def test_threads_limit_the_alignments_that_run_at_the_same_time(tmp_path, threads, expected):
    obj, _db = project(tmp_path, [], threads=threads)
    log = tmp_path / "events.txt"
    commands = [[sys.executable, "-c", SLEEPER, str(log)] for _ in range(5)]
    logs = [str(tmp_path / "align{0}.log".format(i)) for i in range(5)]
    names = [("REF{0}".format(i), "x") for i in range(5)]
    started = time.time()
    assert obj.run_alignments(commands, logs, names) == 0
    assert len(open(log).read().splitlines()) == 10, "every alignment ran"
    assert most_at_once(log) == expected
    assert time.time() - started < 10
