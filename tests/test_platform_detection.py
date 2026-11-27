"""The platform of read input (#119): the clr preset, reading the platform from the head of a
read file, and --platform auto, which uses the platform of each sample."""

import gzip

import pytest

from CanSNPer2.modules.ParsePileup import AUTO, DEFAULT_PLATFORM, PLATFORM_CHOICES, PLATFORMS
from CanSNPer2.modules.ReadInput import (
    ReadInputError,
    pair_read_files,
    read_samplesheet,
)
from CanSNPer2.modules.ReadPlatform import (
    HIFI_BASEQ,
    LONG_READ,
    NO_BASEQ,
    RECORDS,
    describe,
    detect_platform,
    reads_line_fields,
)

"""Phred+33 quality characters: Q40, Q20, Q10 and Q0"""
Q40, Q20, Q10, Q0 = "I", "5", "+", "!"


def records(n, length, quality, name="r{i}"):
    """n FASTQ records of one length, every base of the same quality"""
    return "".join(
        "@{name}\n{seq}\n+\n{q}\n".format(name=name.format(i=i), seq="A" * length, q=quality * length)
        for i in range(n)
    )


def fastq(tmp_path, name, text, gz=False):
    path = tmp_path / name
    if gz:
        with gzip.open(path, "wt") as out:
            out.write(text)
    else:
        path.write_text(text)
    return str(path)


def test_clr_is_a_platform_of_its_own():
    """PacBio continuous long reads carry no base qualities, so the hifi preset discards every
    base; they need their own preset (#119)"""
    clr = PLATFORMS["clr"]
    assert clr["preset"] == "map-pb"
    assert clr["min_baseq"] == 0, "the qualities of CLR reads carry no information"
    assert clr["min_depth"] == PLATFORMS["hifi"]["min_depth"], (
        "measured on eight CLR read sets (#93): the depth of the other platforms gives the same "
        "calls from 20x and adds calls at 10x, with no call on another branch"
    )
    assert PLATFORMS["hifi"] == {
        "preset": "map-hifi",
        "min_mapq": 20,
        "min_baseq": 20,
        "min_depth": 5,
        "min_allele_fraction": 0.8,
        "max_reads": 300000,
    }, "hifi keeps its values: no command changes"


def test_auto_is_a_choice_but_not_a_preset():
    assert AUTO not in PLATFORMS, "every lookup of a preset must see a platform"
    assert PLATFORM_CHOICES == sorted(PLATFORMS) + [AUTO]
    assert "clr" in PLATFORM_CHOICES


@pytest.mark.parametrize(
    "n, length, quality, platform",
    [
        (200, 150, Q40, "illumina"),
        (200, 75, Q20, "illumina"),
        (200, 5000, Q10, "ont"),
        (200, 1200, Q20, "ont"),
        (200, 15000, Q40, "hifi"),
        (200, 10000, Q0, "clr"),
    ],
)
def test_the_platform_is_read_from_the_lengths_and_the_qualities(tmp_path, n, length, quality, platform):
    found = detect_platform(fastq(tmp_path, "r.fastq", records(n, length, quality)), "reads")
    assert found.platform == platform
    assert found.source == "reads"
    assert found.median_length == length
    assert found.reads == n


def test_the_thresholds_are_the_ones_the_detection_documents():
    assert LONG_READ == 500 and HIFI_BASEQ == 25 and NO_BASEQ == 3


@pytest.mark.parametrize(
    "name, length, quality, platform",
    [
        ("m64011_190830_220126/{i}/ccs", 3000, Q10, "hifi"),
        ("m64011_190830_220126/{i}/0_3000", 3000, Q40, "clr"),
        (
            "01234567-89ab-cdef-0123-456789abcdef runid=abc sampleid=s{i}",
            3000,
            Q40,
            "ont",
        ),
        ("M01234:12:000000000-ABCDE:1:1101:15589:{i} 1:N:0:1", 150, Q40, "illumina"),
    ],
)
def test_read_names_that_name_the_instrument_decide(tmp_path, name, length, quality, platform):
    """the names say which instrument wrote the file; modern ONT reads reach the base quality of
    HiFi reads, and PacBio CLR files sometimes carry qualities"""
    found = detect_platform(fastq(tmp_path, "r.fastq", records(50, length, quality, name=name)), "reads")
    assert found.platform == platform
    assert found.source == "read_names"


def test_a_name_form_that_says_illumina_is_not_trusted_on_long_reads(tmp_path):
    """something has rewritten the names of a file of long reads"""
    name = "M01234:12:000000000-ABCDE:1:1101:15589:{i}"
    found = detect_platform(fastq(tmp_path, "r.fastq", records(50, 9000, Q40, name=name)), "reads")
    assert (found.platform, found.source) == ("hifi", "reads")


def test_short_nanopore_reads_are_told_by_their_names(tmp_path):
    """a nanopore run of short fragments is still ONT: the names decide, not the lengths. The
    name is the one dorado writes, the read identifier followed by tab separated tags"""
    name = "01234567-89ab-cdef-0123-456789abcde{i}\tqs:f:12.8\tch:i:181\tRG:Z:x_dna_r10.4.1_e8.2"
    found = detect_platform(fastq(tmp_path, "r.fastq", records(10, 300, Q20, name=name)), "reads")
    assert (found.platform, found.source) == ("ont", "read_names")


def test_a_gzip_compressed_file_is_read(tmp_path):
    path = fastq(tmp_path, "r.fastq.gz", records(100, 8000, Q10), gz=True)
    assert detect_platform(path, "reads").platform == "ont"


def test_only_the_head_of_a_file_is_read(tmp_path):
    """the records after the first RECORDS do not change the result"""
    text = records(RECORDS + 200, 8000, Q10) + records(5000, 150, Q40)
    found = detect_platform(fastq(tmp_path, "big.fastq", text), "reads")
    assert found.platform == "ont"
    assert found.reads == RECORDS


@pytest.mark.parametrize("text", ["", ">c\nACGT\n", "not a read file\n"])
def test_a_file_that_holds_no_reads_gives_no_platform(tmp_path, text):
    found = detect_platform(fastq(tmp_path, "odd.txt", text), "reads")
    assert found.platform is None


def test_a_missing_file_gives_no_platform(tmp_path):
    assert detect_platform(str(tmp_path / "gone.fastq"), "reads").platform is None


def test_the_evidence_is_given_in_words_and_as_fields(tmp_path):
    found = detect_platform(fastq(tmp_path, "r.fastq", records(100, 8000, Q10)), "reads")
    assert "median read length 8000 bp" in describe(found)
    assert "mean base quality 10.0" in describe(found)
    fields = reads_line_fields(found)
    assert " detected=ont" in fields and " detected_from=reads" in fields
    assert " detected_median_length=8000" in fields and " detected_mean_baseq=10.0" in fields
    assert reads_line_fields(None) == ""


def test_reads_without_qualities_are_said_to_have_none(tmp_path):
    found = detect_platform(fastq(tmp_path, "r.fastq", records(100, 9000, Q0)), "reads")
    assert found.platform == "clr"
    assert "no base qualities" in describe(found)


def test_auto_takes_the_platform_of_each_file(tmp_path):
    short = [
        fastq(tmp_path, "s_R1_001.fastq", records(100, 150, Q40)),
        fastq(tmp_path, "s_R2_001.fastq", records(100, 150, Q40)),
    ]
    long_read = fastq(tmp_path, "lr.fastq", records(100, 9000, Q10))
    samples = pair_read_files(short + [long_read], platform=AUTO)
    assert [(s.name, s.platform, len(s.files)) for s in samples] == [
        ("s", "illumina", 2),
        ("lr", "ont", 1),
    ]
    assert samples[0].detected.platform == "illumina"


def test_auto_does_not_pair_long_reads_named_as_mates(tmp_path):
    """only short reads are paired; two long read files named _R1/_R2 are two samples"""
    files = [
        fastq(tmp_path, "x_R1_001.fastq", records(100, 9000, Q10)),
        fastq(tmp_path, "x_R2_001.fastq", records(100, 9000, Q10)),
    ]
    samples = pair_read_files(files, platform=AUTO)
    assert [(s.name, s.platform) for s in samples] == [("x_R1_001", "ont"), ("x_R2_001", "ont")]


def test_auto_falls_back_to_the_default_platform(tmp_path):
    """a file whose platform cannot be told is typed with the default, not refused"""
    empty = fastq(tmp_path, "e.fastq", "@r\n\n+\n\n")
    samples = pair_read_files([empty], platform=AUTO)
    assert samples[0].platform == DEFAULT_PLATFORM
    assert samples[0].detected.platform is None


def test_a_samplesheet_row_may_say_auto(tmp_path):
    fastq(tmp_path, "one.fastq", records(100, 9000, Q10))
    fastq(tmp_path, "two.fastq", records(100, 150, Q40))
    sheet = tmp_path / "samples.tsv"
    sheet.write_text("sample\treads1\tplatform\na\tone.fastq\tauto\nb\ttwo.fastq\tont\n")
    samples = read_samplesheet(str(sheet))
    assert [(s.name, s.platform) for s in samples] == [("a", "ont"), ("b", "ont")]
    assert samples[1].detected is None, "a row that names its platform is not examined"


def test_a_samplesheet_takes_auto_from_the_command_line(tmp_path):
    fastq(tmp_path, "one.fastq", records(100, 15000, Q40))
    sheet = tmp_path / "samples.tsv"
    sheet.write_text("sample\treads1\na\tone.fastq\n")
    samples = read_samplesheet(str(sheet), platform=AUTO)
    assert [(s.name, s.platform) for s in samples] == [("a", "hifi")]


def test_an_unknown_platform_names_the_choices(tmp_path):
    with pytest.raises(ReadInputError, match="clr"):
        pair_read_files(["a.fastq"], platform="pacbio")
    sheet = tmp_path / "samples.tsv"
    fastq(tmp_path, "one.fastq", records(10, 150, Q40))
    sheet.write_text("sample\treads1\tplatform\na\tone.fastq\tpacbio\n")
    with pytest.raises(ReadInputError, match="auto"):
        read_samplesheet(str(sheet))
