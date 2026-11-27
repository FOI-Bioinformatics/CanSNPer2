"""A sample whose result from an earlier run in the same output folder is used again (#117) is
reported as that run reported it: its status, reasons, support and site states, so that it is
counted under its status and drawn in the tree of summary.html. That it comes from an earlier run
is a mark on the sample, not its status."""

import json
import os
import re

import pytest

from test_result_json import validate
from test_typing_inputs import project, write

## called canSNPs and site states per sample (A.1 > B.1, B.2; C.1 elsewhere in the tree)
STATES = {
    "a": (["A.1", "B.1"], {"A.1": 1, "B.1": 1, "B.2": 2, "C.1": 2}),
    "b": (["A.1", "B.2"], {"A.1": 1, "B.2": 1, "B.1": 2}),
    "r": (["A.1", "B.1", "C.1"], {"A.1": 1, "B.1": 1, "C.1": 1, "B.2": 2}),
    "n": ([], {"A.1": 2, "B.1": 2}),
    "c": (["C.1"], {"C.1": 1, "A.1": 2}),
}


def run(tmp_path, names, jobs=1, **kwargs):
    """type the samples of names in one output folder; the states of each sample are looked up
    by its name, so that they follow it into a worker with --jobs"""
    queries = [write(tmp_path / "in" / (name + ".fasta"), ">%s\nACGT\n" % name) for name in names]
    obj, db = project(tmp_path, queries, min_required_hits=1, jobs=jobs, summary=True, **kwargs)
    obj.align = lambda q, references=(): ["REF_x.paf"]

    def states(**_kwargs):
        called, sites = STATES[obj.query_name]
        rows = [[name, "REF", "chr1", "1", "chr1", "1", "C", "T", "T", ""] for name in called]
        return dict(sites), rows, list(called)

    obj.find_snps_multiproc = states
    obj.run(database=db)
    return obj


def without_mark(report):
    return {key: value for key, value in report.items() if key not in ("reused", "earlier_run")}


@pytest.fixture
def earlier(tmp_path):
    """run 1 types a, b, r and n; run 2 adds c to the same folder"""
    first = run(tmp_path, "abrn")
    again = run(tmp_path, "abrnc")
    return first, again


def test_the_first_run_covers_every_status(earlier):
    first, _again = earlier
    assert {first.sample_reports[name]["severity"] for name in "abrn"} >= {"clean", "review", "not called"}


def test_a_reused_sample_is_reported_as_the_earlier_run_reported_it(earlier):
    first, again = earlier
    assert again.reused == {"a", "b", "r", "n"}
    for name in "abrn":
        assert without_mark(again.sample_reports[name]) == without_mark(first.sample_reports[name]), name
        assert again.sample_reports[name]["reused"] is True
    assert again.sample_reports["c"].get("reused") is False


def test_a_reused_sample_has_its_states_for_the_tree(earlier):
    _first, again = earlier
    report = again.sample_reports["a"]
    assert report["states"] == {"A.1": 1, "B.1": 1, "B.2": 2, "C.1": 2}
    assert report["called"] == ["A.1", "B.1"]
    assert report["reasons"] == []
    assert not any("earlier run" in reason for reason in again.sample_reports["r"]["reasons"])


def test_the_earlier_run_is_dated(earlier):
    first, again = earlier
    date = json.load(open(os.path.join(first.results_dir, "a.json")))["date"]
    assert again.sample_reports["a"]["earlier_run"] == date


def test_run_json_counts_reused_samples_apart_from_their_status(earlier):
    first, again = earlier
    document = json.load(open(os.path.join(again.outdir, "run.json")))
    statuses = {name: first.sample_reports[name]["severity"] for name in "abrn"}
    statuses["c"] = "clean"
    counts = document["counts"]
    for status in ("clean", "review", "not called", "failed"):
        assert counts[status] == sum(1 for s in statuses.values() if s == status), status
    assert counts["reused"] == 4 and counts["samples"] == 5
    entries = {s["sample"]: s for s in document["samples"]}
    assert {name: entries[name]["status"] for name in entries} == statuses
    assert [entries[name]["reused"] for name in "abrnc"] == [True, True, True, True, False]
    assert entries["a"]["result"] == "samples/a.json"


def test_run_json_with_reused_samples_matches_the_schema(earlier):
    _first, again = earlier
    document = json.load(open(os.path.join(again.outdir, "run.json")))
    validate(document, "cansnper2-run")
    broken = json.loads(json.dumps(document))
    del broken["samples"][0]["reused"]
    with pytest.raises(Exception, match="reused"):
        validate(broken, "cansnper2-run")


def test_the_console_marks_the_earlier_run_and_gives_the_hint_once(tmp_path, capsys):
    run(tmp_path, "ab")
    capsys.readouterr()
    run(tmp_path, "abc")
    out = capsys.readouterr().out.splitlines()
    assert "a: B.1 (earlier run)" in out and "b: B.2 (earlier run)" in out and "c: C.1" in out
    hints = [line for line in out if "--rerun" in line]
    assert hints == [
        "2 samples were taken from an earlier run in this folder; use --rerun to type them again"
    ]


def test_the_alignment_line_counts_only_the_samples_that_are_aligned(tmp_path, capsys):
    run(tmp_path, "ab")
    first = [line for line in capsys.readouterr().out.splitlines() if line.startswith("Run alignments")]
    assert first[0].startswith("Run alignments for 2 queries against") and "earlier run" not in first[0]
    run(tmp_path, "abc")
    again = [line for line in capsys.readouterr().out.splitlines() if line.startswith("Run alignments")]
    assert again[0].startswith("Run alignments for 1 of 3 queries against")
    assert again[0].endswith("; 2 have a result from an earlier run")


def test_a_run_without_reused_samples_gives_no_hint(tmp_path, capsys):
    run(tmp_path, "ab")
    assert not any("--rerun" in line for line in capsys.readouterr().out.splitlines())


def test_reused_samples_are_not_counted_as_skipped(earlier, caplog):
    _first, again = earlier
    assert again.skipped == {}


def test_jobs_restore_reused_samples_too(tmp_path):
    first = run(tmp_path, "abrn", jobs=2)
    again = run(tmp_path, "abrnc", jobs=3)
    assert list(again.sample_reports) == ["a", "b", "r", "n", "c"]
    for name in "abrn":
        assert without_mark(again.sample_reports[name]) == without_mark(first.sample_reports[name]), name
        assert again.sample_reports[name]["reused"] is True


# without the sample JSON of the earlier run there is nothing to restore


def test_without_the_sample_json_the_sample_is_shown_as_reused(tmp_path):
    first = run(tmp_path, "ab")
    os.remove(os.path.join(first.results_dir, "a.json"))
    again = run(tmp_path, "ab")
    report = again.sample_reports["a"]
    assert report["severity"] == "reused" and report["reused"] is True
    assert report["final"] == "B.1" and report["states"] == {}
    assert again.sample_reports["b"]["severity"] == first.sample_reports["b"]["severity"]


def test_a_sample_json_of_another_input_is_not_restored(tmp_path):
    first = run(tmp_path, "ab")
    path = os.path.join(first.results_dir, "a.json")
    document = json.load(open(path))
    document["input"]["sha256"] = "0" * 64
    json.dump(document, open(path, "w"))
    again = run(tmp_path, "ab")
    assert again.sample_reports["a"]["severity"] == "reused"


def test_a_sample_json_with_another_call_is_not_restored(tmp_path):
    first = run(tmp_path, "ab")
    path = os.path.join(first.results_dir, "a.json")
    document = json.load(open(path))
    document["final_snp"] = "C.1"
    json.dump(document, open(path, "w"))
    again = run(tmp_path, "ab")
    assert again.sample_reports["a"]["severity"] == "reused"
    assert again.sample_reports["a"]["final"] == "B.1"


def test_an_unreadable_sample_json_is_not_restored(tmp_path):
    first = run(tmp_path, "ab")
    with open(os.path.join(first.results_dir, "a.json"), "w") as out:
        out.write("{ not json")
    again = run(tmp_path, "ab")
    assert again.sample_reports["a"]["severity"] == "reused"


# summary.html


def test_the_dashboard_marks_samples_from_an_earlier_run(earlier):
    _first, again = earlier
    page = open(os.path.join(again.outdir, "summary.html")).read()
    payload = json.loads(page.split("var DATA = ", 1)[1].split(";\n", 1)[0])
    rows = {s["name"]: s for s in payload["samples"]}
    assert rows["a"]["severity"] == "clean" and rows["a"]["reused"] is True and rows["a"]["states"]
    assert rows["c"]["reused"] is False
    key = page.split('id="statusKey"', 1)[1].split("</details>", 1)[0]
    assert "earlier run" in key


@pytest.fixture
def page(tmp_path):
    """summary.html of a run with a sample typed in this run, one taken from an earlier run and one
    taken from an earlier result without its sample JSON"""
    from test_report_pages import dashboard, reports

    rows = reports()
    fresh, earlier_sample = dict(rows[0], name="fresh", reused=False), dict(rows[0], name="earlier")
    earlier_sample.update(reused=True, earlier_run="2026-11-20 08:15")
    old = dict(rows[0], name="old", severity="reused", reused=True, states={})
    return dashboard(tmp_path, samples=[fresh, earlier_sample, old])


def needs_node():
    from test_report_pages import NODE

    if NODE is None:
        pytest.skip("node is not installed")


def test_the_page_counts_samples_from_an_earlier_run_apart(page):
    needs_node()
    from test_report_pages import run_page, text_of

    out = run_page(page)
    assert text_of(out["html"]["counts"]) == (
        "3 samples: 2 clean, 0 review, 0 not called, 0 failed, 1 reused; 2 from an earlier run"
    )


def test_the_page_marks_the_row_and_the_detail_of_an_earlier_run(page):
    needs_node()
    from test_report_pages import run_page, text_of

    out = run_page(page, "paint(1);")
    rows = {text_of(re.findall(r"<td[^>]*>(.*?)</td>", r["html"])[0]).strip(): r["html"] for r in out["rows"]}
    assert "earlier run" in rows["earlier"] and "2026-11-20 08:15" in rows["earlier"]
    assert "earlier run" not in text_of(rows["fresh"])
    assert "earlier run" not in text_of(rows["old"]), "the reused badge already says it"
    detail = text_of(out["html"]["detail"])
    assert "Typed by an earlier run in this folder (2026-11-20 08:15); --rerun types it again." in detail
