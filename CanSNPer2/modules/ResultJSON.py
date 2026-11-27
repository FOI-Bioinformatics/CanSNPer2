"""Machine-readable results (#97): samples/<sample>.json for each typed sample, run.json for the
run, and cansnper2_mqc.tsv, a MultiQC custom content table.

The documents are built from the same values as the text results: the sample report that
summary.html shows (CanSNPer2.add_sample_report), the run information of snps_summary.txt
(CanSNPer2.run_info) and the "Typing settings:" and "Reads:" lines of the result files, so that
the JSON and the text files cannot disagree. Their layout is described by the JSON Schemas in
docs/schemas and in docs/output-formats.md; FORMAT_VERSION changes when a field changes meaning
or is removed, not when one is added."""

import json
import os

from CanSNPer2.modules.NewickTree import SNP_KINDS

FORMAT_VERSION = "1.0"

## the site states as the alignment readers give them: 1 derived, 2 ancestral, 3 another base, 0 no data
STATE_NAMES = SNP_KINDS


def number(text):
    """an int or float when the text is one, the text otherwise"""
    for kind in (int, float):
        try:
            return kind(text)
        except ValueError:
            pass
    return text


def fields(text):
    """the key=value words of a line as a dict, numbers converted"""
    values = {}
    for word in text.split():
        key, sep, value = word.partition("=")
        if sep:
            values[key] = number(value)
    return values


def settings_fields(line):
    """the settings of a "Typing settings: database sha256:<sha> key=value ..." line; returns the
    sha256 of the database and the settings"""
    text = line.split(":", 1)[1] if line.startswith("Typing settings:") else line
    sha = None
    for word in text.split():
        if word.startswith("sha256:"):
            sha = word.split(":", 1)[1]
    return sha, fields(text)


def support_fields(support):
    """found/expected as {"found": n, "expected": m}, None for NA"""
    found, sep, expected = str(support).partition("/")
    if not sep:
        return None
    try:
        return {"found": int(found), "expected": int(expected)}
    except ValueError:
        return None


def database_fields(info, sha256=None):
    """the database, scheme and reference genomes of the run information"""
    return {
        "file": info.get("database"),
        "sha256": sha256 or info.get("database_sha256"),
        "scheme": info.get("scheme") or None,
        "scheme_status": info.get("status") or None,
        "scheme_line": info.get("scheme_line") or None,
        "reference_genomes": list(info.get("references") or []),
        "notice": info.get("notice") or None,
    }


def sample_document(
    report, info, settings_line, query_files, query_sha256, kind, reads_line, result_files, competing=()
):
    """samples/<sample>.json of one sample.
    report: the sample report of summary.html; info: the run information; settings_line and
    reads_line: the lines of the result file (reads_line None for an assembly); result_files:
    the paths of the text result files, relative to the output folder; competing: the
    confirmed paths at the same depth of a sample left NA by them (#91)"""
    sha, settings = settings_fields(settings_line)
    final = report.get("final")
    states = report.get("states") or {}
    return {
        "format": "cansnper2-sample",
        "format_version": FORMAT_VERSION,
        "cansnper2_version": info.get("version"),
        "date": info.get("date"),
        "sample": report["name"],
        "input": {
            "kind": kind,
            "files": [os.path.abspath(f) for f in query_files],
            "sha256": query_sha256,
        },
        "status": report.get("severity"),
        "final_snp": final if final and final != "NA" else None,
        "support": support_fields(report.get("support")),
        "snp_path": [] if report.get("path") in (None, "", "NA") else report["path"].split(";"),
        "path_snps": list(report.get("path_snps") or []),
        "reasons": list(report.get("reasons") or []),
        "path_no_data": list(report.get("path_no_data") or []),
        "edge_flagged": list(report.get("edge") or []),
        "derived_cansnps": list(report.get("called") or []),
        "no_data_cansnps": list(report.get("no_data") or []),
        "competing_paths": [
            {
                "node": path["node"],
                "support": {"found": path["found"], "expected": path["expected"]},
                "path_no_data": list(path.get("no_data") or []),
            }
            for path in competing
        ],
        "site_states": {snp: STATE_NAMES.get(state, "no_data") for snp, state in sorted(states.items())},
        "reads": fields(reads_line.split(":", 1)[1]) if reads_line else None,
        "settings": settings,
        "database": database_fields(info, sha),
        "result_files": dict(result_files),
    }


def run_document(info, reports, sample_files, counts):
    """run.json: the run information and one entry per sample, in the order of the run.
    sample_files: the samples/<sample>.json path of each sample that has one, relative to
    the output folder"""
    settings = {
        "min_required_hits": info.get("min_required_hits"),
        "strictness": info.get("strictness"),
        "mask": info.get("mask"),
        "aligner": info.get("aligner"),
    }
    read_settings = [fields(line) for line in info.get("read_settings") or []]
    return {
        "format": "cansnper2-run",
        "format_version": FORMAT_VERSION,
        "cansnper2_version": info.get("version"),
        "date": info.get("date"),
        "database": database_fields(info),
        "settings": settings,
        "read_settings": read_settings,
        "counts": dict(counts),
        "samples": [
            {
                "sample": report["name"],
                "status": report.get("severity"),
                "final_snp": report["final"] if report.get("final") not in (None, "", "NA") else None,
                "support": support_fields(report.get("support")),
                "snp_path": [] if report.get("path") in (None, "", "NA") else report["path"].split(";"),
                "reasons": list(report.get("reasons") or []),
                "result": sample_files.get(report["name"]),
            }
            for report in reports
        ],
    }


MQC_HEADER = """# id: 'cansnper2'
# section_name: 'CanSNPer2'
# description: 'canSNP typing: the final SNP of each sample, its support (derived canSNPs on the SNP path out of the canSNPs on it) and its status in summary.html.'
# format: 'tsv'
# plot_type: 'table'
"""


def mqc_table(reports):
    """cansnper2_mqc.tsv: a MultiQC custom content table, one row per sample"""
    lines = [MQC_HEADER.rstrip("\n"), "\t".join(["Sample", "Final SNP", "Support", "Status"])]
    for report in reports:
        lines.append(
            "\t".join(
                [
                    report["name"],
                    report.get("final") or "NA",
                    str(report.get("support") or "NA"),
                    report.get("severity") or "",
                ]
            )
        )
    return "\n".join(lines) + "\n"


def write_json(path, document):
    """write a document, through a temporary file so that a reader never sees half of it"""
    part = path + ".part"
    with open(part, "w") as out:
        json.dump(document, out, indent=2)
        out.write("\n")
    os.replace(part, path)
