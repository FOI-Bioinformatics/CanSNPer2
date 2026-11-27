"""samples/<sample>.json, run.json and cansnper2_mqc.tsv (#97): written from the same values as
the text result files, and valid against the JSON Schemas in docs/schemas."""

import json
import os
import pathlib

import pytest

from CanSNPer2.modules import ResultJSON
from test_result_files import comments, table_rows
from test_typing_inputs import fake_alignment, project, write

ROOT = pathlib.Path(__file__).resolve().parent.parent
SCHEMAS = ROOT / "docs" / "schemas"


def validate(document, name):
    """check a document against docs/schemas/<name>.schema.json"""
    jsonschema = pytest.importorskip("jsonschema")
    referencing = pytest.importorskip("referencing")
    schemas = {path.name: json.loads(path.read_text()) for path in SCHEMAS.glob("*.schema.json")}
    registry = referencing.Registry().with_resources(
        (key, referencing.Resource.from_contents(value)) for key, value in schemas.items()
    )
    schema = schemas[name + ".schema.json"]
    jsonschema.Draft202012Validator(schema, registry=registry).validate(document)


def typed(tmp_path, called, states, queries=None, **kwargs):
    q = write(tmp_path / "in" / "x.fasta")
    obj, db = project(tmp_path, queries or [q], min_required_hits=2, **kwargs)
    fake_alignment(obj, called, states, q)
    obj.run(database=db)
    return obj, q


def load(path):
    return json.loads(pathlib.Path(path).read_text())


def test_the_sample_json_gives_what_the_result_file_gives(tmp_path):
    obj, q = typed(tmp_path, ["A.1", "B.1"], {"A.1": 1, "B.1": 1, "B.2": 2})
    document = load(os.path.join(obj.results_dir, "x.json"))
    validate(document, "cansnper2-sample")
    notes = comments(os.path.join(obj.results_dir, "x_snps.txt"))
    assert document["sample"] == "x" and document["status"] == "clean"
    assert document["final_snp"] == "B.1"
    assert document["support"] == {"found": 2, "expected": 2}
    assert "# SNP path: " + ";".join(document["snp_path"]) in notes
    assert document["input"] == {
        "kind": "assembly",
        "files": [os.path.abspath(q)],
        "sha256": obj.file_sha256(q),
    }
    settings_line = [line[2:] for line in notes if line.startswith("# Typing settings:")][0]
    sha, settings = ResultJSON.settings_fields(settings_line)
    assert document["settings"] == settings and document["database"]["sha256"] == sha
    assert document["site_states"]["B.2"] == "ancestral" and document["site_states"]["A.1"] == "derived"
    assert document["derived_cansnps"] == ["A.1", "B.1"]
    assert document["result_files"] == {
        "snps": "samples/x_snps.txt",
        "not_called": "samples/x_not_called.txt",
    }
    assert document["reads"] is None


def test_a_sample_without_a_call_has_no_final_snp(tmp_path):
    obj, q = typed(tmp_path, [], {"A.1": 2})
    document = load(os.path.join(obj.results_dir, "x.json"))
    validate(document, "cansnper2-sample")
    assert document["status"] == "not called"
    assert document["final_snp"] is None and document["support"] is None and document["snp_path"] == []
    assert document["reasons"]


def test_run_json_lists_every_sample_with_the_failed_one(tmp_path):
    q = write(tmp_path / "in" / "x.fasta")
    (tmp_path / "in" / "d").mkdir()
    obj, q = typed(
        tmp_path,
        ["A.1", "B.1"],
        {"A.1": 1, "B.1": 1},
        queries=[q, str(tmp_path / "in" / "d")],
        keep_going=True,
    )
    document = load(os.path.join(obj.outdir, "run.json"))
    validate(document, "cansnper2-run")
    assert [(s["sample"], s["status"], s["final_snp"], s["result"]) for s in document["samples"]] == [
        ("x", "clean", "B.1", "samples/x.json"),
        ("d", "failed", None, None),
    ]
    assert document["counts"] == {
        "samples": 2,
        "clean": 1,
        "review": 0,
        "not called": 0,
        "failed": 1,
        "reused": 0,
    }
    assert [s["reused"] for s in document["samples"]] == [False, False]
    assert document["settings"] == {
        "min_required_hits": 2,
        "strictness": 0.7,
        "mask": 0,
        "aligner": obj.aligner.name,
    }
    assert document["read_settings"] == []
    assert not os.path.exists(os.path.join(obj.results_dir, "d.json"))


def test_a_reused_sample_keeps_its_json(tmp_path):
    obj, q = typed(tmp_path, ["A.1", "B.1"], {"A.1": 1, "B.1": 1})
    path = os.path.join(obj.results_dir, "x.json")
    first = pathlib.Path(path).read_text()
    again, db = project(tmp_path, [q], min_required_hits=2)
    fake_alignment(again, ["A.1", "B.1"], {"A.1": 1, "B.1": 1}, q)
    again.run(database=db)
    assert "x" in again.reused
    assert pathlib.Path(path).read_text() == first
    document = load(os.path.join(again.outdir, "run.json"))
    validate(document, "cansnper2-run")
    """the sample keeps the status of the earlier run and is marked as reused (#117)"""
    assert document["samples"][0]["status"] == "clean" and document["samples"][0]["reused"] is True
    assert document["samples"][0]["result"] == "samples/x.json"
    assert document["samples"][0]["support"] == {"found": 2, "expected": 2}


def test_no_snp_files_writes_no_sample_json(tmp_path):
    obj, q = typed(tmp_path, ["A.1", "B.1"], {"A.1": 1, "B.1": 1}, export=False)
    assert not os.path.exists(os.path.join(obj.results_dir, "x.json"))
    document = load(os.path.join(obj.outdir, "run.json"))
    assert document["samples"][0]["result"] is None and document["samples"][0]["final_snp"] == "B.1"


def test_the_multiqc_table_is_written_with_summary(tmp_path):
    obj, q = typed(tmp_path, ["A.1", "B.1"], {"A.1": 1, "B.1": 1}, summary=True)
    text = (pathlib.Path(obj.outdir) / "cansnper2_mqc.tsv").read_text()
    assert text.startswith("# id: 'cansnper2'\n")
    assert "# plot_type: 'table'" in text
    assert table_rows(os.path.join(obj.outdir, "cansnper2_mqc.tsv")) == [
        ["Sample", "Final SNP", "Support", "Status"],
        ["x", "B.1", "2/2", "clean"],
    ]


def test_without_summary_there_is_no_multiqc_table(tmp_path):
    obj, q = typed(tmp_path, ["A.1", "B.1"], {"A.1": 1, "B.1": 1})
    assert not os.path.exists(os.path.join(obj.outdir, "cansnper2_mqc.tsv"))


@pytest.mark.parametrize(
    "text,expected",
    [
        (
            "Reads: platform=ont files=1 sites_with_reads=3/4 median_depth=20.5 mixed_sites=0",
            {
                "platform": "ont",
                "files": 1,
                "sites_with_reads": "3/4",
                "median_depth": 20.5,
                "mixed_sites": 0,
            },
        ),
        (
            "min_required_hits=3 strictness=0.7 mask=0 aligner=minimap2",
            {"min_required_hits": 3, "strictness": 0.7, "mask": 0, "aligner": "minimap2"},
        ),
    ],
)
def test_key_value_lines_are_read_with_their_numbers(text, expected):
    assert ResultJSON.fields(text.split(":", 1)[1] if text.startswith("Reads:") else text) == expected


@pytest.mark.parametrize(
    "support,expected", [("2/3", {"found": 2, "expected": 3}), ("NA", None), (None, None)]
)
def test_support_is_split_into_found_and_expected(support, expected):
    assert ResultJSON.support_fields(support) == expected


def test_the_schema_check_rejects_a_document_that_does_not_match(tmp_path):
    jsonschema = pytest.importorskip("jsonschema")
    obj, q = typed(tmp_path, ["A.1", "B.1"], {"A.1": 1, "B.1": 1})
    document = load(os.path.join(obj.results_dir, "x.json"))
    del document["support"]
    with pytest.raises(jsonschema.ValidationError):
        validate(document, "cansnper2-sample")


def test_the_sample_json_lists_competing_paths(tmp_path):
    obj, q = typed(tmp_path, ["A.1", "B.1", "B.2"], {"A.1": 1, "B.1": 1, "B.2": 1})
    document = load(os.path.join(obj.results_dir, "x.json"))
    validate(document, "cansnper2-sample")
    assert document["status"] == "not called" and document["final_snp"] is None
    assert document["competing_paths"] == [
        {"node": "B.1", "support": {"found": 2, "expected": 2}, "path_no_data": []},
        {"node": "B.2", "support": {"found": 2, "expected": 2}, "path_no_data": []},
    ]
    called, _q = typed(tmp_path / "other", ["A.1", "B.1"], {"A.1": 1, "B.1": 1})
    assert load(os.path.join(called.results_dir, "x.json"))["competing_paths"] == []
