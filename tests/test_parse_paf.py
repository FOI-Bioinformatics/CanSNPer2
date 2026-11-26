"""The minimap2 aligner and its PAF reader. The database stores each SNP as a reference
sequence and a position on it, the coordinates PAF gives for the target, so the reader needs
nothing but the alignment. The cs tag is on the reference forward strand, also for "-"
alignments, and the query position counts down on those."""

import os
import random
import shutil

import pytest

from helpers import make_db
from CanSNPer2.modules.Aligners import Minimap2, get_aligner
from CanSNPer2.modules.CanSNPer2 import CanSNPer2
from CanSNPer2.modules.ParsePAF import ParsePAF


def paf_line(cs, strand="+", tname="chr2", tstart=9, qstart=20, mapq=60, tp="P", qname="ctg1"):
    """one PAF record; the reference span and the query span follow from the cs string"""
    import re

    ref_len = (
        sum(len(m) for m in re.findall(r"=([A-Za-z]+)", cs))
        + cs.count("*")
        + sum(len(m) for m in re.findall(r"-([a-z]+)", cs))
    )
    query_len = (
        sum(len(m) for m in re.findall(r"=([A-Za-z]+)", cs))
        + cs.count("*")
        + sum(len(m) for m in re.findall(r"\+([a-z]+)", cs))
    )
    fields = [
        qname,
        "1000",
        qstart,
        qstart + query_len,
        strand,
        tname,
        "50",
        tstart,
        tstart + ref_len,
        ref_len,
        ref_len,
        mapq,
        "tp:A:" + tp,
        "cs:Z:" + cs,
    ]
    return "\t".join(str(f) for f in fields)


def read(tmp_path, lines, rows, name="a", **kwargs):
    """read PAF lines against the reference REF; rows as for helpers.make_db"""
    db, paf = tmp_path / (name + ".db"), tmp_path / ("REF_" + name + ".paf")
    make_db(db, rows)
    paf.write_text("\n".join(lines) + "\n")
    parser = ParsePAF(**kwargs)
    snps = parser.run(str(paf), reference="REF", database=str(db))
    return snps, {row[0]: row for row in parser.get_snp_info()}


# chr2 from its 10th base: ACGT, then chr2:14 is c in the reference, g in the query
ROWS = [("B.1", "chr2", 14, "C", "G")]


def test_target_coordinates_are_the_database_coordinates(tmp_path):
    snps, info = read(tmp_path, [paf_line("=ACGT*cg=AAAA")], ROWS)
    assert snps == {"B.1": 1}
    assert info["B.1"][:6] == ["B.1", "REF", "chr2", "14", "ctg1", "25"], (
        "qstart 20 + 4 bases before it, 1-based"
    )
    assert info["B.1"][8] == "G"


def test_paf_reads_target_coordinates_without_a_table(tmp_path):
    # no <alignment>.sequences file and no reference file: the target name and the 0-based
    # target start of the record are enough
    snps, info = read(tmp_path, [paf_line("=AC*cg=AAAA")], [("B.1", "chr2", 12, "C", "G")], name="notable")
    assert not (tmp_path / "REF_notable.paf.sequences").exists()
    assert snps == {"B.1": 1}
    assert info["B.1"][2:4] == ["chr2", "12"]


def test_same_position_on_two_sequences_both_called(tmp_path):
    # position 50 on chr1 and position 50 on chr2 are different SNPs
    rows = [("A.1", "chr1", 50, "A", "G"), ("B.2", "chr2", 50, "A", "C")]
    lines = [
        paf_line("=ACGT*ag=AA", tname="chr1", tstart=45, qname="ctg1"),
        paf_line("=ACGT*ac=AA", tname="chr2", tstart=45, qname="ctg2"),
    ]
    snps, info = read(tmp_path, lines, rows, name="two")
    assert snps == {"A.1": 1, "B.2": 1}
    assert info["A.1"][2:6] == ["chr1", "50", "ctg1", "25"]
    assert info["B.2"][2:6] == ["chr2", "50", "ctg2", "25"]
    assert info["A.1"][8] == "G" and info["B.2"][8] == "C"


def test_a_sequence_without_snps_is_skipped(tmp_path):
    snps, info = read(tmp_path, [paf_line("=ACGT*cg=AAAA", tname="plasmid")], ROWS, name="other")
    assert snps == {} and info == {}


def test_query_insertion_and_reverse_strand(tmp_path):
    snps, info = read(tmp_path, [paf_line("=ACGT+tt*cg=AAAA")], ROWS, name="ins")
    assert info["B.1"][5] == "27", "two inserted query bases move the query position, not the reference one"
    # "-": the cs string is still on the reference strand, the query runs down from qend 31
    snps, info = read(tmp_path, [paf_line("=ACGT+tt*cg=AAAA", strand="-")], ROWS, name="rev")
    assert snps == {"B.1": 1}
    assert info["B.1"][5] == "25", "qend 31 - 6 query bases before it"


def test_deletion_at_a_snp_is_no_data(tmp_path):
    snps, info = read(tmp_path, [paf_line("=ACGT-c=AAAA")], ROWS)
    assert snps == {"B.1": 0}
    assert info["B.1"][5] == "N/A" and info["B.1"][8] == "-"


def test_only_primary_alignments_with_enough_mapping_quality_are_read(tmp_path):
    lines = [
        paf_line("=ACGT*cg=AAAA"),
        paf_line("=ACGT*ca=AAAA", tp="S", qname="paralog"),  ## secondary
        paf_line("=ACGT*ct=AAAA", mapq=3, qname="repeat"),
    ]  ## primary, ambiguous placement
    snps, info = read(tmp_path, lines, ROWS)
    assert snps == {"B.1": 1} and info["B.1"][4] == "ctg1"
    snps, info = read(tmp_path, lines, ROWS, name="mapq0", min_mapq=0)
    assert snps == {"B.1": 0} and info["B.1"][8] == "N", (
        "two primary alignments disagree, no evidence either way"
    )


def test_missing_cs_tag_is_an_error(tmp_path):
    with pytest.raises(ValueError, match="cs tag"):
        read(tmp_path, [paf_line("=ACGT*cg=AAAA").rsplit("\t", 1)[0]], ROWS, name="nocs")


def test_minimap2_definition(tmp_path):
    mm2 = get_aligner("minimap2", "/opt/bin/")
    assert isinstance(mm2, Minimap2)
    assert mm2.command("r.fna", "q.fna", "o.paf") == [
        "/opt/bin/minimap2",
        "-c",
        "--cs=long",
        "-x",
        "asm20",
        "-r",
        "1k,10k",
        "-N",
        "50",
        "-p",
        "0.1",
        "-o",
        "o.paf",
        "r.fna",
        "q.fna",
    ]
    assert mm2.alignment_path("/tmp/run", "CO92", "s") == "/tmp/run/CO92_s.paf"
    assert isinstance(mm2.reader(), ParsePAF) and mm2.reader().min_mapq == 20
    fasta = tmp_path / "r.fna"
    fasta.write_text(">c1 x\nACGT\nAC\n>c2\nGG\n")
    mm2.prepare(str(fasta), "q.fna", str(tmp_path / "o.paf"))
    assert not (tmp_path / "o.paf.sequences").exists(), (
        "PAF carries the target coordinates, no sequence table is written"
    )


def revcomp(seq):
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


@pytest.mark.skipif(shutil.which("minimap2") is None, reason="minimap2 is not installed")
def test_minimap2_through_the_whole_alignment_step(tmp_path):
    # a two-sequence reference and a query split into three contigs, one of them reverse
    # complemented and the order shuffled; SNPs on both reference sequences
    rng = random.Random(7)
    chr1 = "".join(rng.choice("ACGT") for _ in range(6000))
    chr2 = "".join(rng.choice("ACGT") for _ in range(4000))
    joined = chr1 + chr2
    ## the query is built from the joined sequence; the SNPs are stored per reference sequence
    snps = {1500: ("B.1", "chr1", 1500), 4500: ("B.2", "chr1", 4500), 6000 + 2000: ("B.3", "chr2", 2000)}
    rows, query = [], list(joined)
    for pos, (snp_id, sequence, seq_pos) in snps.items():
        ref_base = joined[pos - 1]
        alt = "A" if ref_base != "A" else "C"
        derived = snp_id != "B.2"  ## B.2 stays ancestral in the query
        rows.append((snp_id, sequence, seq_pos, ref_base, alt))
        if derived:
            query[pos - 1] = alt
    query = "".join(query)
    contigs = [("q1", query[:3000]), ("q2", revcomp(query[3000:7000])), ("q3", query[7000:])]
    refdir = tmp_path / "refs"
    refdir.mkdir()
    (refdir / "REF.fna").write_text(">chr1\n{0}\n>chr2\n{1}\n".format(chr1, chr2))
    qfile = tmp_path / "sample.fna"
    qfile.write_text(
        "".join(">{0}\n{1}\n".format(name, seq) for name, seq in (contigs[2], contigs[0], contigs[1]))
    )
    db = str(tmp_path / "t.db")
    make_db(db, rows)
    obj = CanSNPer2(
        [str(qfile)],
        workdir=str(tmp_path),
        refdir="refs",
        tmpdir="shared",
        logdir="logs",
        outdir="out",
        database=db,
        aligner="minimap2",
    )
    files = obj.align(str(qfile))
    assert files and all(os.path.isfile(f) for f in files)
    states, info, called = obj.find_snps_multiproc(
        xmfa_obj=obj.aligner.reader(database=db), xmfa_files=files, references=obj.alignment_references
    )
    assert states == {"B.1": 1, "B.2": 2, "B.3": 1}
    assert sorted(called) == ["B.1", "B.3"]
    assert {row[0]: (row[2], row[3]) for row in info} == {
        "B.1": ("chr1", "1500"),
        "B.2": ("chr1", "4500"),
        "B.3": ("chr2", "2000"),
    }
    where = {row[0]: (row[4], int(row[5])) for row in info}
    assert where["B.1"] == ("q1", 1500)
    assert where["B.2"] == ("q2", 7000 - 4500 + 1), "q2 is reverse complemented: query 4500 is its base 2501"
    assert where["B.3"] == ("q3", 8000 - 7000)


def short_form(cs):
    """the cs string in the short form: ":n" replaces each run of identical bases"""
    import re

    return re.sub(r"=([A-Za-z]+)", lambda m: ":" + str(len(m.group(1))), cs)


def with_cs(line, old, new):
    return line.replace("cs:Z:" + old, "cs:Z:" + new)


def test_short_cs_form_gives_the_same_result_as_the_long_form(tmp_path):
    long_cs = "=ACGT*cg=AAAA"
    assert short_form(long_cs) == ":4*cg:4"
    long_snps, long_info = read(tmp_path, [paf_line(long_cs)], ROWS, name="lng")
    short_snps, short_info = read(
        tmp_path, [with_cs(paf_line(long_cs), long_cs, ":4*cg:4")], ROWS, name="shrt"
    )
    assert short_snps == long_snps == {"B.1": 1}
    assert short_info["B.1"][2:] == long_info["B.1"][2:]


@pytest.mark.parametrize("strand", ["+", "-"])
def test_short_cs_form_with_indels_before_the_snp(tmp_path, strand):
    # an insertion and a deletion before the SNPs, a match run after them
    long_cs = "=AC+tt=GT-ga=CC*cg=AA"
    assert short_form(long_cs) == ":2+tt:2-ga:2*cg:2"
    rows = [
        ("B.1", "chr2", 18, "C", "G"),
        ("B.2", "chr2", 12, "G", "A"),
        ("B.3", "chr2", 16, "C", "T"),
        ("B.4", "chr2", 14, "G", "A"),
    ]
    tag = strand.replace("-", "r").replace("+", "f")
    line = paf_line(long_cs, strand=strand)
    long_snps, long_info = read(tmp_path, [line], rows, name="l" + tag)
    short_snps, short_info = read(
        tmp_path, [with_cs(line, long_cs, ":2+tt:2-ga:2*cg:2")], rows, name="s" + tag
    )
    assert set(long_snps) == {"B.1", "B.2", "B.3", "B.4"} and long_snps["B.1"] == 1
    assert short_snps == long_snps
    assert set(short_info) == set(long_info)
    for key in long_info:
        assert short_info[key][2:] == long_info[key][2:]


def test_skipped_records_are_reported_when_none_is_used(tmp_path, caplog):
    import logging

    secondary = paf_line("=ACGT*cg=AAAA", tp="S")
    truncated = "ctg1\t1000\t20"
    with caplog.at_level(logging.WARNING):
        snps, info = read(tmp_path, [secondary, truncated], ROWS, name="skp")
    assert snps == {}
    assert any("2 record" in r.getMessage() for r in caplog.records if r.levelno == logging.WARNING)


def test_no_warning_when_some_records_are_used(tmp_path, caplog):
    import logging

    with caplog.at_level(logging.WARNING):
        read(tmp_path, [paf_line("=ACGT*cg=AAAA"), paf_line("=ACGT*cg=AAAA", tp="S")], ROWS, name="ok")
    assert not any("skipped" in r.getMessage() for r in caplog.records)
