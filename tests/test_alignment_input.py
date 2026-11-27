"""Read typing from BAM and CRAM files (#96): reads already mapped to a reference genome of the
scheme are counted at the canSNP positions without mapping them again."""

import logging
import os
import shutil
import subprocess

import pytest

from CanSNPer2.modules.CanSNPer2 import CanSNPer2Error
from CanSNPer2.modules.ReadInput import input_kind
from test_read_typing import comments, make_scheme, run, write_pair

pytestmark = pytest.mark.skipif(
    not (shutil.which("minimap2") and shutil.which("samtools")), reason="minimap2 and samtools are needed"
)


def mapped(tmp_path, name, reads, sort=True, cram=False, rename=None):
    """the reads mapped to the reference of the toy scheme, as BAM (or CRAM)"""
    ref = tmp_path / "refs" / "REF.fna"
    if rename:
        ref_text = ref.read_text().replace(">chr", ">" + rename, 1)
        ref = tmp_path / "renamed.fna"
        ref.write_text(ref_text)
    sam = subprocess.run(
        ["minimap2", "-a", "-x", "sr", str(ref)] + [str(r) for r in reads],
        capture_output=True,
        check=True,
    ).stdout
    out = tmp_path / (name + (".cram" if cram else ".bam"))
    if sort:
        command = ["samtools", "sort", "-o", str(out)]
        if cram:
            command += ["-O", "cram", "--reference", str(ref)]
        command.append("-")
    else:
        command = ["samtools", "view", "-b", "-o", str(out), "-"]
    subprocess.run(command, input=sam, check=True, capture_output=True)
    return str(out)


def site_counts(tmp_path, name):
    return (tmp_path / "out" / "samples" / (name + "_site_counts.tsv")).read_text()


@pytest.fixture
def scheme(tmp_path):
    db, a, _ = make_scheme(tmp_path)
    reads = write_pair(tmp_path, "sa", [(a, 30, 1)])
    return db, reads


def test_bam_and_cram_files_are_told_by_their_content(tmp_path, scheme):
    db, reads = scheme
    bam = mapped(tmp_path, "sb", reads)
    cram = mapped(tmp_path, "sc", reads, cram=True)
    assert input_kind(bam) == "alignment" and input_kind(cram) == "alignment"
    assert input_kind(reads[0]) == "reads"
    os.rename(bam, tmp_path / "no_ending")
    assert input_kind(str(tmp_path / "no_ending")) == "alignment"


@pytest.mark.parametrize(
    "sort,cram", [(True, False), (False, False), (True, True)], ids=["bam", "unsorted bam", "cram"]
)
def test_a_mapped_file_is_typed_like_the_reads(tmp_path, scheme, sort, cram):
    db, reads = scheme
    fastq = run(tmp_path, reads, db=db)
    expected = site_counts(tmp_path, "sa")
    alignment = mapped(tmp_path, "sa_mapped", reads, sort=sort, cram=cram)
    obj = run(tmp_path, [alignment], db=db)
    assert obj.typed == {"sa_mapped": fastq.typed["sa"]} == {"sa_mapped": "A.3"}
    assert site_counts(tmp_path, "sa_mapped") == expected
    notes = comments(tmp_path / "out" / "samples" / "sa_mapped_snps.txt")
    settings = [line for line in notes if line.startswith("Typing settings:")][0]
    assert "max_reads=0" in settings, "every read of the file is used"
    reads_line = [line for line in notes if line.startswith("Reads: ")][0]
    assert "files=1" in reads_line


def test_a_bam_and_read_files_are_typed_in_one_run(tmp_path, scheme):
    db, reads = scheme
    alignment = mapped(tmp_path, "sb", reads)
    obj = run(tmp_path, reads + [alignment], db=db)
    assert obj.typed == {"sa": "A.3", "sb": "A.3"}


def test_a_bam_of_another_reference_is_refused(tmp_path, scheme):
    db, reads = scheme
    alignment = mapped(tmp_path, "sb", reads, rename="contig_7")
    with pytest.raises(CanSNPer2Error, match="use --keep-going"):
        run(tmp_path, [alignment], db=db)
    obj = run(tmp_path, [alignment], db=db, keep_going=True)
    reason = obj.not_typed["sb"]
    assert "contig_7" in reason and "REF" in reason and "chr" in reason


def test_local_references_need_the_reads(tmp_path, scheme):
    db, reads = scheme
    alignment = mapped(tmp_path, "sb", reads)
    with pytest.raises(CanSNPer2Error, match="--local-references.*BAM or CRAM"):
        run(tmp_path, [alignment], db=db, local_references=True)


def test_max_reads_is_not_applied_to_a_mapped_file(tmp_path, scheme, caplog):
    db, reads = scheme
    alignment = mapped(tmp_path, "sb", reads)
    caplog.set_level(logging.WARNING)
    run(tmp_path, [alignment], db=db, max_reads=100)
    assert "--max-reads is not applied to BAM or CRAM" in caplog.text
    notes = comments(tmp_path / "out" / "samples" / "sb_snps.txt")
    assert any("max_reads=0" in line for line in notes if line.startswith("Typing settings:"))


def test_an_assembly_and_a_bam_cannot_be_typed_in_one_run(tmp_path, scheme):
    db, reads = scheme
    alignment = mapped(tmp_path, "sb", reads)
    genome = tmp_path / "ga.fna"
    genome.write_text((tmp_path / "refs" / "REF.fna").read_text())
    with pytest.raises(CanSNPer2Error, match="reads and assemblies cannot be typed in one run"):
        run(tmp_path, [str(genome), alignment], db=db)


@pytest.mark.parametrize("name,sample", [("x.bam", "x"), ("x.sorted.bam", "x.sorted"), ("x.CRAM", "x")])
def test_the_sample_is_named_after_the_file(tmp_path, scheme, name, sample):
    db, reads = scheme
    alignment = mapped(tmp_path, "tmp_name", reads, cram=name.lower().endswith(".cram"))
    os.rename(alignment, tmp_path / name)
    obj = run(tmp_path, [str(tmp_path / name)], db=db)
    assert list(obj.typed) == [sample]
