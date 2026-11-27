"""Read typing: the samtools mpileup reader. Reads are mapped to each reference and
samtools mpileup lists the read bases at the canSNP positions; the reader counts them and
decides each canSNP by depth and allele fraction (PileupReader, classify_counts)."""

import pytest

from helpers import make_db
from CanSNPer2.modules.AlignmentReader import NO_DATA, DERIVED, ANCESTRAL, OTHER_BASE
from CanSNPer2.modules.ParsePileup import (
    PLATFORMS,
    PileupReader,
    classify_counts,
    parse_bases,
    read_site_counts,
)


# parse_bases: the base column of mpileup


def test_plain_bases_both_strands():
    assert parse_bases("AaaGg") == {"A": 3, "C": 0, "G": 2, "T": 0}


def test_read_start_with_mapping_quality_and_read_end():
    # ^ is followed by the mapping quality character, which can itself be a base letter or "$"
    assert parse_bases("^AC^$t$A$") == {"A": 1, "C": 1, "G": 0, "T": 1}


def test_insertions_and_deletions_after_a_base_are_skipped():
    assert parse_bases("A+2CCg-12ACGTACGTACGTt") == {"A": 1, "C": 0, "G": 1, "T": 1}


def test_deletion_placeholders_and_n_are_not_counted():
    assert parse_bases("A*#Nn>") == {"A": 1, "C": 0, "G": 0, "T": 0}


def test_dots_and_commas_are_the_reference_base():
    assert parse_bases(".,.G", ref="c") == {"A": 0, "C": 3, "G": 1, "T": 0}
    assert parse_bases(".,", ref="N") == {"A": 0, "C": 0, "G": 0, "T": 0}, (
        "no reference given, nothing to count"
    )


def test_zero_depth():
    assert parse_bases("*") == {"A": 0, "C": 0, "G": 0, "T": 0}
    assert parse_bases("") == {"A": 0, "C": 0, "G": 0, "T": 0}


# classify_counts: the per-site rule


def counts(**kw):
    return {b: kw.get(b, 0) for b in "ACGT"}


@pytest.mark.parametrize(
    "c, label, state",
    [
        (counts(G=10), "derived", DERIVED),
        (counts(C=10), "ancestral", ANCESTRAL),
        (counts(G=9, C=1), "derived", DERIVED),  # 0.9 >= 0.8
        (counts(G=6, C=4), "mixed", NO_DATA),  # neither allele reaches 0.8
        (counts(T=10), "other_base", OTHER_BASE),
        (counts(T=5, G=5), "mixed", NO_DATA),
        (counts(G=4), "low_depth", NO_DATA),  # below min_depth 5
        (counts(), "no_reads", NO_DATA),
    ],
)
def test_classify_counts(c, label, state):
    assert classify_counts(c, ancestral="C", derived="G", min_depth=5, min_allele_fraction=0.8)[:2] == (
        label,
        state,
    )


def test_classify_counts_gives_the_dominant_base_and_fraction():
    label, state, base, depth, fraction = classify_counts(counts(G=9, C=1), "C", "G", 5, 0.8)
    assert (base, depth, fraction) == ("G", 10, 0.9)
    assert classify_counts(counts(G=6, C=4), "C", "G", 5, 0.8)[2] == "N", (
        "an undecided site has no single base"
    )


def test_platform_defaults():
    assert set(PLATFORMS) == {"illumina", "ont", "hifi", "clr"}
    assert PLATFORMS["illumina"]["preset"] == "sr"
    assert PLATFORMS["ont"]["preset"] == "map-ont"
    assert PLATFORMS["hifi"]["preset"] == "map-hifi"
    assert PLATFORMS["clr"]["preset"] == "map-pb"
    for p in PLATFORMS.values():
        assert {"min_mapq", "min_baseq", "min_depth", "min_allele_fraction"} <= set(p)


# PileupReader on mpileup output

ROWS = [
    ("B.1", "chr", 100, "C", "G"),
    ("B.2", "chr", 200, "A", "T"),
    ("B.3", "chr", 300, "G", "A"),
    ("B.4", "chr", 400, "T", "C"),
    ("B.5", "chr", 500, "A", "G"),
]


def pileup(*lines):
    return "\n".join("\t".join(str(f) for f in line) for line in lines) + "\n"


def read(tmp_path, text, rows=ROWS, **kwargs):
    db, path = tmp_path / "a.db", tmp_path / "REF_s1.pileup"
    make_db(db, rows)
    path.write_text(text)
    kwargs.setdefault("min_depth", 5)
    kwargs.setdefault("min_allele_fraction", 0.8)
    reader = PileupReader(**kwargs)
    snps = reader.run(str(path), reference="REF", database=str(db))
    return reader, snps, {row[0]: row for row in reader.get_snp_info()}, path


TEXT = pileup(
    ("chr", 100, "N", 8, "GGGGgggg", "IIIIIIII"),  # derived
    ("chr", 200, "N", 6, "AAAaaa", "IIIIII"),  # ancestral
    ("chr", 300, "N", 10, "GGGGGAAAAA", "I" * 10),  # mixed
    ("chr", 400, "N", 2, "Cc", "II"),  # low depth
    ("chr", 500, "N", 0, "*", "*"),  # no reads
)


def test_states_rows_and_called(tmp_path):
    reader, snps, info, _ = read(tmp_path, TEXT)
    assert snps == {"B.1": DERIVED, "B.2": ANCESTRAL, "B.3": NO_DATA, "B.4": NO_DATA}
    assert reader.get_called_snps() == ["B.1"]
    assert "B.5" not in snps, "a site without reads is no data, as an unaligned position of an assembly"
    assert info["B.1"] == ["B.1", "REF", "chr", "100", "reads", "8", "C", "G", "G", ""]
    assert info["B.3"][8] == "N" and info["B.3"][5] == "10"


def test_site_counts_list_every_site_with_its_label(tmp_path):
    reader, _, _, path = read(tmp_path, TEXT)
    rows = read_site_counts(str(path) + ".counts")
    assert [r["Name"] for r in rows] == ["B.1", "B.2", "B.3", "B.4", "B.5"]
    assert [r["state"] for r in rows] == ["derived", "ancestral", "mixed", "low_depth", "no_reads"]
    assert rows[2]["A"] == "5" and rows[2]["G"] == "5" and rows[2]["depth"] == "10"
    assert rows[0]["fraction_derived"] == "1.00"
    assert reader.mixed_sites == ["B.3"]


def test_a_site_missing_from_the_pileup_is_no_reads(tmp_path):
    reader, snps, _, path = read(tmp_path, pileup(("chr", 100, "N", 8, "GGGGgggg", "IIIIIIII")))
    assert snps == {"B.1": DERIVED}
    assert [r["state"] for r in read_site_counts(str(path) + ".counts")] == ["derived"] + ["no_reads"] * 4


def test_lines_of_other_sequences_and_positions_are_ignored(tmp_path):
    text = pileup(("plasmid", 100, "N", 8, "GGGGGGGG", "I" * 8), ("chr", 101, "N", 8, "GGGGGGGG", "I" * 8))
    _, snps, _, _ = read(tmp_path, text)
    assert snps == {}


def test_thresholds_are_options(tmp_path):
    _, snps, _, _ = read(tmp_path, TEXT, min_depth=2, min_allele_fraction=0.5)
    assert snps["B.4"] == DERIVED, "depth 2 is enough with min_depth 2"
    assert snps["B.3"] == DERIVED, "5 of 10 reaches a fraction of 0.5"


def test_a_malformed_line_is_an_error(tmp_path):
    with pytest.raises(ValueError, match="mpileup"):
        read(tmp_path, "chr\t100\tN\n")
