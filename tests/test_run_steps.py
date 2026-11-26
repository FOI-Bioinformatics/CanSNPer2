"""The steps of a typing run (CanSNPer2._run split into methods, #103), each tested on its own
without an aligner: the call of a sample from given site states, the review reasons of the
results overview and the Reads line of a read sample."""

import logging
from types import SimpleNamespace

from CanSNPer2.modules.CanSNPer2 import SampleCall
from test_typing_inputs import project, write

# rows in the columns of SNP_FILE_COLUMNS, as the alignment readers give them
ROWS = [[name, "REF", "chr1", "1", "chr1", "1", "C", "T", "T", ""] for name in ("A.1", "B.1")]


def sample(tmp_path, **kwargs):
    q = write(tmp_path / "x.fasta")
    obj, db = project(tmp_path, [q], min_required_hits=2, **kwargs)
    obj.query_name = "x"
    obj.read_database()
    return obj, q


def a_call(snp="B.1", message="", path=("A.1", "B.1"), off_chain=(), support="2/2"):
    return SampleCall(
        snp=snp,
        message=message,
        final_snp=None,
        no_data_snps=[],
        path_names=list(path),
        off_chain=list(off_chain),
        support=support,
        final_line="",
        path_text="",
        called_rows=[],
    )


def test_call_sample_gives_the_final_snp_path_and_support(tmp_path):
    obj, _q = sample(tmp_path)
    call = obj.call_sample({"A.1": 1, "B.1": 1}, ROWS, ["A.1", "B.1"])
    assert call.snp == "B.1"
    assert call.path_names == ["A.1", "B.1"] and call.off_chain == []
    assert call.support == "2/2"
    ## the root of the test tree has no canSNP: depth 3, two canSNPs on the path
    assert call.final_line == "Final SNP: B.1 found/depth: 2/3 found/expected: 2/2"
    assert [row[0] for row in call.called_rows] == ["A.1", "B.1"]
    assert obj.result_paths["x"] == "A.1;B.1" and obj.sample_supports["x"] == "2/2"


def test_call_sample_without_a_call_is_na_with_the_reason(tmp_path):
    obj, _q = sample(tmp_path)
    call = obj.call_sample({"A.1": 2, "B.1": 2}, ROWS, [])
    assert call.snp == "NA" and call.support == "NA"
    assert call.final_line == "Final SNP: NA"
    assert call.message
    assert obj.sample_supports["x"] == "NA"


def test_a_clean_call_needs_no_review(tmp_path):
    obj, _q = sample(tmp_path)
    assert obj.review_reasons(a_call(), set(), [], []) == ("clean", [])


def test_a_sample_without_a_call_is_not_called_with_its_reason(tmp_path):
    obj, _q = sample(tmp_path)
    severity, reasons = obj.review_reasons(a_call(snp="NA", message="no path", path=()), set(), [], [])
    assert severity == "not called" and reasons == ["no path"]


def test_review_reasons_name_each_cause(tmp_path, caplog):
    obj, _q = sample(tmp_path, mixed_sites=2)
    caplog.set_level(logging.WARNING)
    severity, reasons = obj.review_reasons(
        a_call(off_chain=["C.1"]), edge_flags={"B.1", "C.1"}, path_missing=["A.1"], mixed=["S1", "S2"]
    )
    assert severity == "review"
    assert [reason.split(" ", 2)[:2] for reason in reasons] == [
        ["1", "SNP"],
        ["SNPs", "flagged"],
        ["1", "derived"],
        ["possibly", "mixed:"],
    ]
    assert "B.1" in reasons[1] and "C.1" not in reasons[1], "only edge flags on the path"
    assert "x is possibly mixed" in caplog.text


def test_fewer_mixed_sites_than_the_threshold_need_no_review(tmp_path):
    obj, _q = sample(tmp_path, mixed_sites=3)
    assert obj.review_reasons(a_call(), set(), [], ["S1", "S2"]) == ("clean", [])


def site(depth, state="derived"):
    return {"depth": str(depth), "state": state}


def test_the_reads_line_gives_the_median_depth_over_sites_with_reads(tmp_path):
    obj, q = sample(tmp_path)
    obj.read_samples = {q: SimpleNamespace(platform="ont", files=[q], name="x")}
    rows = [site(0, "no_data"), site(10), site(20), site(31, "mixed")]
    line = obj.reads_line(q, rows, ["S4"], None)
    assert line == "Reads: platform=ont files=1 sites_with_reads=3/4 median_depth=20 mixed_sites=1"


def test_the_reads_line_names_the_local_reference(tmp_path):
    obj, q = sample(tmp_path)
    obj.read_samples = {q: SimpleNamespace(platform="illumina", files=[q, q], name="x")}
    obj.local_references = True
    rows = [site(4), site(6, "conflict")]
    assert obj.reads_line(q, rows, [], None).endswith("median_depth=5 mixed_sites=0 local_reference=none")
    done = {"genome": "GCA_1.1", "node": "B.1", "share": 0.5, "screens": 2, "rounds": 1}
    assert obj.reads_line(q, rows, [], done).endswith(
        " local_reference=GCA_1.1 local_node=B.1 screen_share=0.50 screens=2 local_rounds=1 conflicts=1"
    )
