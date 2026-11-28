"""The order of the result rows and of the derived canSNPs (#126). The rows arrive from one
worker process per alignment, in the order the workers happen to report them, so two runs of the
same input wrote the rows of <sample>_not_called.txt, the derived_cansnps of the sample JSON and
the data of summary.html in different orders when the machine was loaded. They are put in a fixed
order where they are collected."""

import json

from helpers import make_db
from CanSNPer2.modules.CanSNPer2 import CanSNPer2, sort_site_rows
from test_typing_inputs import project, write


def row(name, reference="REF", sequence="chr1", position="100"):
    """a result row in the columns of SNP_FILE_COLUMNS"""
    return [name, reference, sequence, position, "q1", "7", "C", "T", "T", ""]


def test_the_rows_are_ordered_by_genome_sequence_and_position():
    rows = [
        row("B.2", "OSU18", "CP000437.1", "19938"),
        row("A.1", "SCHUS4.2", "AJ749949.2", "83976"),
        row("B.1", "OSU18", "CP000437.1", "1710718"),
        row("B.9", "OSU18", "chrII", "7"),
    ]
    assert [r[0] for r in sort_site_rows(rows)] == ["B.2", "B.1", "B.9", "A.1"]


def test_the_position_is_ordered_as_a_number():
    rows = [row("A.1", position="9"), row("A.2", position="10"), row("A.3", position="100")]
    assert [r[0] for r in sort_site_rows(list(reversed(rows)))] == ["A.1", "A.2", "A.3"]


def test_two_rows_of_one_position_are_ordered_by_node():
    rows = [row("B.2"), row("B.1")]
    assert [r[0] for r in sort_site_rows(rows)] == ["B.1", "B.2"]


def test_a_row_without_a_number_for_a_position_is_kept():
    rows = [row("A.2", position="100"), row("A.1", position="N/A")]
    ordered = sort_site_rows(rows)
    assert len(ordered) == 2 and {r[0] for r in ordered} == {"A.1", "A.2"}


class Reader(object):
    """an AlignmentReader whose workers report their rows in the order the references are given,
    so that the collected order can be scrambled on purpose"""

    def __init__(self, database=None, **kwargs):
        self.database = database
        self.SNPS = {}
        self.SNP_info = []
        self.called_snps = []
        self.reference = "REF"

    def run(self, xmfa_file, *args, **kwargs):
        import os

        self.reference = os.path.basename(xmfa_file).split("_")[0]
        position = {"CREF": "300", "AREF": "100", "BREF": "200"}[self.reference]
        node = {"CREF": "C.1", "AREF": "A.1", "BREF": "B.1"}[self.reference]
        self.SNPS = {node: 1}
        self.SNP_info = [[node, self.reference, "chr", position, "q", "1", "C", "T", "T", ""]]
        self.called_snps = [node]

    def get_snps(self):
        return self.SNPS

    def get_snp_info(self):
        return self.SNP_info

    def get_called_snps(self):
        return self.called_snps


def test_the_collected_rows_do_not_depend_on_which_worker_reports_first(tmp_path):
    """the workers are started in the order of the alignment files and report when they are done;
    the collected rows and derived canSNPs are in a fixed order whichever that is"""
    db = str(tmp_path / "s.db")
    make_db(db, [("A.1", "chr", 100, "C", "T"), ("B.1", "chr", 200, "C", "T"), ("C.1", "chr", 300, "C", "T")])
    alignments = [str(tmp_path / "{0}_q.paf".format(ref)) for ref in ("CREF", "AREF", "BREF")]
    for path in alignments:
        open(path, "w").write("")
    obj = CanSNPer2(["q.fna"], workdir=str(tmp_path), tmpdir="tmp", logdir="logs", outdir="out")
    _snps, snp_info, called = obj.find_snps_multiproc(xmfa_obj=Reader(database=db), xmfa_files=alignments)
    assert [r[0] for r in snp_info] == ["A.1", "B.1", "C.1"], "by reference genome, not by worker"
    assert called == ["A.1", "B.1", "C.1"]


ROWS = [
    ["B.1", "REF", "chr1", "200", "q1", "7", "C", "T", "T", ""],
    ["A.1", "REF", "chr1", "100", "q1", "7", "C", "T", "T", ""],
    ["B.2", "REF", "chr1", "300", "q1", "7", "C", "T", "C", ""],
]


def typed(tmp_path, scrambled):
    """Type one sample with given site states and rows, without an aligner. The stub stands for
    the collection of the rows, so it hands over what find_snps_multiproc hands over: the rows and
    the derived canSNPs in the fixed order."""
    q = write(tmp_path / "x.fasta")
    obj, db = project(tmp_path, [q], min_required_hits=1, summary=True)
    obj.align = lambda query, references=(): ["REF_x.paf"]
    obj.find_snps_multiproc = lambda **kwargs: (
        {"A.1": 1, "B.1": 1},
        sort_site_rows(scrambled),
        sorted(["B.1", "A.1"]),
    )
    obj.run(database=db)
    return obj


def test_the_not_called_file_and_the_sample_json_follow_that_order(tmp_path):
    typed(tmp_path, ROWS)
    out = tmp_path / "out" / "samples"
    rows = [
        line.split("\t")[0]
        for line in (out / "x_not_called.txt").read_text().splitlines()[1:]
        if not line.startswith("#")
    ]
    assert rows == ["B.2", "B.2", "C.1"], (
        "the canSNP with the ancestral base, then the two without alignment data, which the tree "
        "lists after the table rows"
    )
    called = [
        line.split("\t")[0]
        for line in (out / "x_snps.txt").read_text().splitlines()[1:]
        if not line.startswith("#")
    ]
    assert called == ["A.1", "B.1"], "the called rows of a sample without a call follow the rows"
    document = json.loads((out / "x.json").read_text())
    assert document["derived_cansnps"] == ["A.1", "B.1"], "the order of the collection, kept"
