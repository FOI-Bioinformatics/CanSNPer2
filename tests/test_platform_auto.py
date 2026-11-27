"""The platform of a read sample in a run (#119): the warning when the reads do not look like
the platform in effect, --platform auto, the clr preset, and what the result files record."""

import json
import logging
import shutil
import subprocess

import pytest

import helpers
from test_read_typing import comments, make_scheme, run, write_pair

pytestmark = pytest.mark.skipif(
    not (shutil.which("minimap2") and shutil.which("samtools")), reason="minimap2 and samtools are needed"
)


def long_reads(tmp_path, genome, name, quality=None, length=2000, depth=30, seed=7):
    """single-end long reads of a genome; quality replaces the quality of every base"""
    text, _ = helpers.simulate_reads(genome, depth, seed=seed, length=length, paired=False)
    if quality:
        lines = text.split("\n")
        lines[3::4] = [quality * len(line) for line in lines[1::4]]
        text = "\n".join(lines)
    path = tmp_path / (name + ".fastq")
    path.write_text(text)
    return str(path)


def reads_line(tmp_path, name):
    lines = comments(tmp_path / "out" / "samples" / (name + "_snps.txt"))
    return [line for line in lines if line.startswith("Reads: ")][0]


def settings_line(tmp_path, name):
    lines = comments(tmp_path / "out" / "samples" / (name + "_snps.txt"))
    return [line for line in lines if line.startswith("Typing settings:")][0]


def warnings_of(caplog):
    return [r.getMessage() for r in caplog.records if r.levelno >= logging.WARNING]


def test_reads_that_do_not_look_like_the_platform_are_warned_about(tmp_path, caplog):
    db, _, b = make_scheme(tmp_path)
    reads = long_reads(tmp_path, b, "lr")
    with caplog.at_level(logging.WARNING):
        obj = run(tmp_path, [reads], db=db)
    assert obj.typed == {"lr": "B.3"}, "the call is not changed by the warning"
    said = [w for w in warnings_of(caplog) if "look like" in w]
    assert len(said) == 1
    assert "hifi" in said[0] and "illumina" in said[0] and "--platform auto" in said[0]
    assert "platform=illumina" in settings_line(tmp_path, "lr"), "--platform is not overridden"


def test_reads_that_match_the_platform_are_not_warned_about(tmp_path, caplog):
    db, a, _ = make_scheme(tmp_path)
    files = write_pair(tmp_path, "sa", [(a, 30, 1)])
    with caplog.at_level(logging.WARNING):
        run(tmp_path, files, db=db)
    assert [w for w in warnings_of(caplog) if "look like" in w] == []


def test_auto_uses_the_platform_of_the_sample(tmp_path, caplog):
    db, _, b = make_scheme(tmp_path)
    reads = long_reads(tmp_path, b, "lr")
    with caplog.at_level(logging.WARNING):
        obj = run(tmp_path, [reads], db=db, platform="auto")
    assert obj.typed == {"lr": "B.3"}
    assert "platform=hifi" in settings_line(tmp_path, "lr")
    assert [w for w in warnings_of(caplog) if "look like" in w] == []


def test_auto_types_short_and_long_reads_in_one_run(tmp_path):
    db, a, b = make_scheme(tmp_path)
    short = write_pair(tmp_path, "sa", [(a, 30, 1)])
    long_read = long_reads(tmp_path, b, "lr")
    obj = run(tmp_path, short + [long_read], db=db, platform="auto")
    assert obj.typed == {"sa": "A.3", "lr": "B.3"}
    assert "platform=illumina" in settings_line(tmp_path, "sa")
    assert "platform=hifi" in settings_line(tmp_path, "lr")


def test_the_reads_line_and_the_sample_json_record_the_detection(tmp_path):
    db, _, b = make_scheme(tmp_path)
    reads = long_reads(tmp_path, b, "lr")
    run(tmp_path, [reads], db=db, platform="hifi")
    line = reads_line(tmp_path, "lr")
    assert "detected=hifi" in line and "detected_from=reads" in line
    assert "detected_median_length=2000" in line and "detected_mean_baseq=40.0" in line
    document = json.loads((tmp_path / "out" / "samples" / "lr.json").read_text())
    assert document["reads"]["detected"] == "hifi"
    assert document["reads"]["detected_median_length"] == 2000
    assert document["reads"]["detected_mean_baseq"] == 40.0


def test_an_assembly_run_has_no_detection(tmp_path):
    db, a, _ = make_scheme(tmp_path)
    query = tmp_path / "asm.fna"
    query.write_text(">chr\n" + a + "\n")
    obj = run(tmp_path, [str(query)], db=db)
    assert obj.typed == {"asm": "A.3"}
    document = json.loads((tmp_path / "out" / "samples" / "asm.json").read_text())
    assert document["reads"] is None


def test_reads_without_base_qualities_are_typed_with_the_clr_preset(tmp_path, caplog):
    db, _, b = make_scheme(tmp_path)
    reads = long_reads(tmp_path, b, "clr", quality="!", length=5000)
    with caplog.at_level(logging.WARNING):
        obj = run(tmp_path, [reads], db=db, platform="clr")
    assert obj.typed == {"clr": "B.3"}
    assert "platform=clr min_depth=10 min_allele_fraction=0.7" in settings_line(tmp_path, "clr")
    assert "min_baseq=0" in settings_line(tmp_path, "clr")
    assert warnings_of(caplog) == []
    assert "detected=clr" in reads_line(tmp_path, "clr")


def test_reads_without_base_qualities_point_at_min_baseq(tmp_path, caplog):
    """the hifi preset filters at base quality 20, which discards every base of such a file"""
    db, _, b = make_scheme(tmp_path)
    reads = long_reads(tmp_path, b, "flat", quality="!", length=5000)
    with caplog.at_level(logging.WARNING):
        run(tmp_path, [reads], db=db, platform="hifi")
    said = warnings_of(caplog)
    assert [w for w in said if "no base qualities" in w and "--min-baseq 0" in w]
    assert [w for w in said if "look like" in w and "clr" in w]


def test_auto_takes_the_platform_of_a_bam_from_its_read_group(tmp_path, caplog):
    db, a, _ = make_scheme(tmp_path)
    files = write_pair(tmp_path, "sa", [(a, 30, 1)])
    ref = tmp_path / "refs" / "REF.fna"
    sam = subprocess.run(
        ## minimap2 takes the read group line with escaped tabs
        ["minimap2", "-a", "-x", "sr", "-R", r"@RG\tID:1\tSM:sa\tPL:ONT", str(ref)] + files,
        capture_output=True,
        check=True,
    ).stdout
    bam = str(tmp_path / "rg.bam")
    subprocess.run(["samtools", "sort", "-o", bam, "-"], input=sam, check=True, capture_output=True)
    with caplog.at_level(logging.WARNING):
        run(tmp_path, [bam], db=db)
    said = [w for w in warnings_of(caplog) if "look like" in w]
    assert len(said) == 1 and "ont" in said[0]
    assert "detected=ont" in reads_line(tmp_path, "rg")
    assert "detected_from=read_group" in reads_line(tmp_path, "rg")
