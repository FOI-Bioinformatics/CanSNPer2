"""The steps of the canSNP file check (SourceCheck.check_snps split into functions, #103),
each tested on its own: the header, one row, the four groups of row checks and the checks that
need the whole file. The messages of a whole file are checked in test_create_input_checks.py."""

import io

from CanSNPer2.modules.DatabaseConnection import CANSNP_SOURCE_COLUMNS, cansnp_row
from CanSNPer2.modules.SourceCheck import (
    MAX_POSITION,
    SnpContext,
    check_snp_bases,
    check_snp_node,
    check_snp_site,
    check_snp_totals,
    check_snp_values,
    snp_header,
    snp_row,
)

ROW = {
    "node": "B.1",
    "rank": "1",
    "reference": "Larkeryd 2014",
    "genome": "SCHU_S4",
    "sequence": "chr1",
    "sequence_position": "100",
    "derived_base": "T",
    "ancestral_base": "C",
    "score": "0.9",
    "typed": "yes",
}


def row(**changes):
    return cansnp_row(dict(ROW, **changes))


def context(**changes):
    fields = dict(
        path="snps.txt",
        tree_nodes=None,
        genomes=None,
        references_path=None,
        tree_name="the tree file",
        existing=None,
    )
    fields.update(changes)
    return SnpContext(**fields)


def handle(*lines):
    return io.StringIO("".join(line + "\n" for line in lines))


def header_line(columns=CANSNP_SOURCE_COLUMNS):
    return "\t".join(columns)


def test_snp_header_counts_the_columns_a_row_must_have():
    problems = []
    header = snp_header("snps.txt", handle(header_line()), problems)
    assert problems == []
    assert header.columns == list(CANSNP_SOURCE_COLUMNS)
    assert header.required == len(CANSNP_SOURCE_COLUMNS)


def test_snp_header_lets_a_row_leave_out_the_trailing_comment_columns():
    problems = []
    header = snp_header(
        "snps.txt", handle(header_line(CANSNP_SOURCE_COLUMNS + ("comment", "notes"))), problems
    )
    assert problems == []
    assert header.required == len(CANSNP_SOURCE_COLUMNS)
    assert len(header.columns) == len(CANSNP_SOURCE_COLUMNS) + 2


def test_snp_header_reports_an_empty_file():
    problems = []
    assert snp_header("snps.txt", handle(""), problems) is None
    assert "the SNP file is empty" in problems[0]


def test_snp_header_reports_a_repeated_column():
    problems = []
    assert snp_header("snps.txt", handle(header_line() + "\tnode"), problems) is None
    assert "more than once" in problems[0]


def test_snp_header_reports_the_joined_position_format():
    problems = []
    columns = tuple("position" if name == "sequence_position" else name for name in CANSNP_SOURCE_COLUMNS)
    assert snp_header("snps.txt", handle(header_line(columns)), problems) is None
    assert "cannot be used" in problems[0]


def test_snp_header_reports_a_missing_column():
    problems = []
    columns = tuple(name for name in CANSNP_SOURCE_COLUMNS if name != "typed")
    assert snp_header("snps.txt", handle(header_line(columns)), problems) is None
    assert "has no column(s) typed" in problems[0]


def test_snp_row_skips_an_empty_line_and_a_comment():
    problems = []
    header = snp_header("snps.txt", handle(header_line()), problems)
    assert snp_row("snps.txt", 2, "\n", header, problems) is None
    assert snp_row("snps.txt", 3, "# a comment\n", header, problems) is None
    assert problems == []


def test_snp_row_reports_a_row_with_too_few_columns():
    problems = []
    header = snp_header("snps.txt", handle(header_line()), problems)
    assert snp_row("snps.txt", 2, "B.1\t1\n", header, problems) is None
    assert problems == ["snps.txt:2: expected 10 columns, found 2"]


def test_snp_row_reads_a_snp_file_row_as_rank_1_typed():
    problems = []
    header = snp_header("snps.txt", handle("snp_id\treference\tgenome\tsequence\tposition"), problems)
    assert header is None  ## position alone is the joined format
    problems = []
    header = snp_header(
        "snps.txt",
        handle("snp_id\treference\tgenome\tsequence\tsequence_position\tderived_base\tancestral_base"),
        problems,
    )
    assert problems == []
    parsed = snp_row("snps.txt", 2, "B.1\tLarkeryd 2014\tSCHU_S4\tchr1\t100\tT\tC\n", header, problems)
    assert problems == []
    assert (parsed["node"], parsed["rank"], parsed["typed"]) == ("B.1", "1", "yes")


def test_check_snp_values_gives_the_rank_and_whether_the_row_is_typed():
    problems = []
    assert check_snp_values(context(), 2, row(), problems) == (1, True)
    assert problems == []
    assert check_snp_values(context(), 2, row(typed="no"), problems) == (1, False)
    assert problems == []


def test_check_snp_values_reports_rank_score_typed_and_reference():
    problems = []
    rank, typed = check_snp_values(
        context(), 2, row(rank="0", score="high", typed="maybe", reference=" "), problems
    )
    assert (rank, typed) == (None, False)
    assert problems == [
        "snps.txt:2: rank '0' is not a positive integer",
        "snps.txt:2: score 'high' is not a number",
        "snps.txt:2: typed 'maybe' is not yes or no",
        "snps.txt:2: reference is empty, it must name where the canSNP is published or NA",
    ]


def test_check_snp_values_allows_an_empty_score():
    problems = []
    assert check_snp_values(context(), 2, row(score=""), problems) == (1, True)
    assert problems == []


def test_check_snp_node_reports_a_node_that_is_not_in_the_tree():
    problems = []
    check_snp_node(context(tree_nodes={"B.2"}), 2, row(), problems, 1, True)
    assert problems == ["snps.txt:2: node 'B.1' does not appear in the tree file"]


def test_check_snp_node_reports_a_rank_and_a_typed_canSNP_the_database_holds():
    problems = []
    existing = ("scheme.db", {("B.1", 1)}, {"B.1"}, set())
    check_snp_node(context(existing=existing), 2, row(), problems, 1, True)
    assert problems == [
        "snps.txt:2: node 'B.1' already has a canSNP of rank 1 in scheme.db",
        "snps.txt:2: node 'B.1' already has a typed canSNP in scheme.db",
    ]


def test_check_snp_node_says_nothing_about_a_rank_that_is_not_a_number():
    problems = []
    existing = ("scheme.db", {("B.1", 1)}, {"B.1"}, set())
    check_snp_node(context(existing=existing), 2, row(), problems, None, False)
    assert problems == []


def test_check_snp_site_gives_the_site_the_database_annotates():
    problems = []
    assert check_snp_site(context(), 2, row(), problems) == ("SCHU_S4", "chr1", 100)
    assert problems == []


def test_check_snp_site_resolves_the_genome_name_to_its_genome():
    problems = []
    site = check_snp_site(
        context(genomes={"GCA_000008985.1": "SCHU_S4"}), 2, row(genome="GCA_000008985.1"), problems
    )
    assert site == ("SCHU_S4", "chr1", 100)
    assert problems == []


def test_check_snp_site_reports_a_genome_that_is_not_listed():
    problems = []
    site = check_snp_site(
        context(genomes={"OSU18": "OSU18"}, references_path="references.txt"), 2, row(), problems
    )
    assert site == ("SCHU_S4", "chr1", 100)
    assert problems == ["snps.txt:2: genome 'SCHU_S4' is not listed in references.txt"]


def test_check_snp_site_reports_an_empty_sequence_and_gives_no_site():
    problems = []
    assert check_snp_site(context(), 2, row(sequence=" "), problems) is None
    assert problems == ["snps.txt:2: sequence is empty, it must name a sequence of the reference genome"]


def test_check_snp_site_reports_a_position_that_is_not_a_positive_integer():
    problems = []
    assert check_snp_site(context(), 2, row(sequence_position="abc"), problems) is None
    assert problems == ["snps.txt:2: sequence_position 'abc' is not a positive integer"]


def test_check_snp_site_reports_a_position_the_database_cannot_hold():
    problems = []
    assert check_snp_site(context(), 2, row(sequence_position=str(MAX_POSITION + 1)), problems) is None
    assert "is larger than" in problems[0]


def test_check_snp_site_reports_a_position_the_database_already_annotates():
    problems = []
    existing = ("scheme.db", set(), set(), {("SCHU_S4", "chr1", 100)})
    assert check_snp_site(context(existing=existing), 2, row(), problems) == ("SCHU_S4", "chr1", 100)
    assert problems == [
        "snps.txt:2: position 100 on sequence 'chr1' of 'SCHU_S4' is already annotated in scheme.db"
    ]


def test_check_snp_bases_reports_a_base_that_is_not_acgt():
    problems = []
    check_snp_bases(context(), 2, row(derived_base="t", ancestral_base="N"), problems)
    assert problems == [
        "snps.txt:2: derived_base 't' is not one of A, C, G, T",
        "snps.txt:2: ancestral_base 'N' is not one of A, C, G, T",
    ]


def test_check_snp_bases_reports_the_same_base_twice():
    problems = []
    check_snp_bases(context(), 2, row(ancestral_base="T"), problems)
    assert problems == ["snps.txt:2: derived and ancestral base are both 'T'"]


def test_check_snp_totals_reports_a_file_without_canSNP_rows():
    problems = []
    check_snp_totals(context(), problems, 0, {}, {}, {})
    assert "has a header but no canSNP rows" in problems[0]


def test_check_snp_totals_reports_ranks_that_do_not_count_from_one():
    problems = []
    check_snp_totals(context(), problems, 2, {"B.1": [1, 3]}, {"B.1": 1}, {})
    assert problems == ["snps.txt: the ranks of node 'B.1' are 1, 3, not 1, 2, 3 ... in order"]


def test_check_snp_totals_reports_a_node_without_and_with_several_typed_canSNPs():
    problems = []
    check_snp_totals(context(), problems, 3, {"B.1": [1], "B.2": [1, 2]}, {"B.2": 2}, {})
    assert problems == [
        "snps.txt: node 'B.1' has no canSNP marked typed yes, CanSNPer2 would not call it",
        "snps.txt: node 'B.2' has 2 canSNPs marked typed yes, CanSNPer2 calls one canSNP per node",
    ]


def test_check_snp_totals_reports_a_position_annotated_twice():
    problems = []
    check_snp_totals(context(), problems, 2, {"B.1": [1]}, {"B.1": 1}, {("SCHU_S4", "chr1", 100): 2})
    assert problems[-1] == ("snps.txt: position 100 on sequence 'chr1' of 'SCHU_S4' is annotated 2 times")
