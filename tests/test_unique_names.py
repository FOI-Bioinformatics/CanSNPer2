"""Names that must be unique.

- Samples are named after their file name. Two different files with the same name
  (a/x.fasta and b/x.fasta) were typed as one sample: the second was skipped as already done
  and the summary reported it with the call of the first. A result file from an earlier run
  was also used for any later input with the same file name.
- Sequences are identified by the first word of their header. A reference with two sequences
  of the same name made the minimap2 reader lose SNPs and read others from the wrong sequence;
  a query with duplicate names gave positions in the wrong contig."""

import os
import sqlite3

import pytest

import helpers
from CanSNPer2.DownloadGenomes import DownloadGenomes
from CanSNPer2.modules.Aligners import Minimap2, ProgressiveMauve
from CanSNPer2.modules.CanSNPer2 import CanSNPer2, CanSNPer2Error
from CanSNPer2.modules.ReferenceCheck import check_reference, duplicate_names, sequence_names


class InstalledMauve(ProgressiveMauve):
    """the run now stops up front when the aligner is not installed; these tests are about what
    happens after that check, so the program is reported as found"""

    def available(self):
        return True


class InstalledMinimap2(Minimap2):
    def available(self):
        return True


GENOME = ">chr1 test\n" + "ACGT" * 250 + "\n"


def make_db(path):
    """one SNP on chr1, the sequence of GENOME, on the root of the tree (a node linked to
    itself); typing refuses SNPs on nodes that are not in the tree"""
    helpers.make_db(path, [("S.1", "chr1", 100, "A", "T")])
    conn = sqlite3.connect(path)
    conn.execute("INSERT INTO tree (parent, child, rank_i) VALUES (1, 1, 1)")
    conn.commit()
    conn.close()


def project(tmp_path, queries, **kwargs):
    db = str(tmp_path / "t.db")
    make_db(db)
    refdir = tmp_path / "refs"
    refdir.mkdir(exist_ok=True)
    (refdir / "REF.fna").write_text(GENOME)
    kwargs.setdefault("aligner", InstalledMauve())
    return CanSNPer2(
        queries,
        workdir=str(tmp_path),
        refdir="refs",
        tmpdir="shared",
        logdir="logs",
        outdir="out",
        database=db,
        **kwargs,
    ), db


def write(path, text=GENOME):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return str(path)


def test_different_files_with_the_same_name_stop_before_anything_is_aligned(tmp_path):
    a = write(tmp_path / "a" / "x.fasta")
    b = write(tmp_path / "b" / "x.fasta", GENOME.replace("ACGT", "ACGA", 1))
    obj, db = project(tmp_path, [a, b], keep_going=True)
    with pytest.raises(CanSNPer2Error, match="x = .*a/x.fasta, .*b/x.fasta"):
        obj.run(database=db)
    # a refused run removes its run directory; if it is kept, it holds no alignment
    assert not os.path.isdir(obj.tmpdir) or not [f for f in os.listdir(obj.tmpdir) if f.endswith(".xmfa")], (
        "nothing was aligned"
    )


def test_the_same_file_given_twice_is_typed_once(tmp_path, caplog):
    a = write(tmp_path / "a" / "x.fasta")
    obj, db = project(tmp_path, [])
    assert obj.check_sample_names([a, a, os.path.join(tmp_path, "a", ".", "x.fasta")]) == [a]
    assert "more than once" in caplog.text


def result_file(obj, name, source_line):
    os.makedirs(obj.results_dir, exist_ok=True)
    path = os.path.join(obj.results_dir, name + "_snps.txt")
    with open(path, "w") as out:
        out.write("SNP path: S.1\nFinal SNP: S.1 found/depth: 1/1\n" + source_line)
    return path


def test_an_earlier_result_is_used_only_for_the_same_input(tmp_path):
    q = write(tmp_path / "in" / "x.fasta")
    obj, db = project(tmp_path, [q], keep_going=True)
    sha = obj.file_sha256(q)
    # made from the same input, here at another path: used again
    result_file(obj, "x", "Query file: /elsewhere/x.fasta sha256:{0}\n".format(sha))
    obj.run(database=db)
    assert "x" in obj.reused and "x" not in obj.not_typed
    # made from a different input: not used, the sample fails with the reason
    obj, db = project(tmp_path, [q], keep_going=True)
    result_file(obj, "x", "Query file: /elsewhere/x.fasta sha256:{0}\n".format("0" * 64))
    obj.run(database=db)
    assert "x" not in obj.reused
    assert "exists for a different input (/elsewhere/x.fasta)" in obj.not_typed["x"]
    # from an older version, which did not record its input: used again as before
    obj, db = project(tmp_path, [q], keep_going=True)
    result_file(obj, "x", "")
    obj.run(database=db)
    assert "x" in obj.reused


def test_the_result_records_its_input(tmp_path):
    q = write(tmp_path / "in" / "x.fasta")
    obj, db = project(tmp_path, [q])
    obj.query_sha256 = obj.file_sha256(q)
    line = obj.query_source_line(q)
    assert line == "Query file: {0} sha256:{1}".format(os.path.abspath(q), obj.query_sha256)
    path = result_file(obj, "x", line + "\n")
    assert obj.read_result_source(path) == (os.path.abspath(q), obj.query_sha256)
    assert obj.read_result_file(path) == ("S.1", "S.1")


def test_a_failed_sample_does_not_show_an_old_result_in_the_summary(tmp_path, monkeypatch):
    # the summary read the result file of every skipped sample, also of one skipped because
    # its alignment failed, which then showed the call an earlier run had left under its name
    q = write(tmp_path / "in" / "x.fasta")
    obj, db = project(tmp_path, [q], keep_going=True)
    result_file(obj, "x", "")
    obj.skipped["x"] = "progressiveMauve produced no alignments"
    monkeypatch.setattr(obj, "create_tree", lambda *args, **kwargs: None)
    obj.print_summary()
    assert summary_rows(obj) == ["sample\tfinal_snp\tsupport\tsnp_path", "x\tNA\tNA\tNA"]
    obj.reused.add("x")
    obj.print_summary()
    assert summary_rows(obj) == ["sample\tfinal_snp\tsupport\tsnp_path", "x\tS.1\t1/1\tS.1"]


def summary_rows(obj):
    """the lines of snps_summary.txt after its # lines on the run"""
    lines = open(os.path.join(obj.outdir, "snps_summary.txt")).read().splitlines()
    return [line for line in lines if not line.startswith("#")]


def test_query_sequence_names_must_be_unique(tmp_path):
    obj, db = project(tmp_path, [])
    dup = write(tmp_path / "dup.fasta", ">contig1 a\nACGT\n>contig1 b\nACGT\n")
    with pytest.raises(CanSNPer2Error, match="more than one sequence named contig1"):
        obj.clean_input_fasta(dup, "query")
    empty = write(tmp_path / "empty.fasta", ">\nACGT\n>contig2\nACGT\n")
    with pytest.raises(CanSNPer2Error, match="1 sequence\\(s\\) without a name"):
        obj.clean_input_fasta(empty, "query")
    fine = write(tmp_path / "fine.fasta", ">contig1 a\nACGT\n>contig2 a\nACGT\n")
    assert obj.clean_input_fasta(fine, "query") == fine


def test_reference_sequence_names_must_be_unique(tmp_path):
    q = write(tmp_path / "in" / "x.fasta")
    obj, db = project(tmp_path, [q])
    (tmp_path / "refs" / "REF.fna").write_text(
        ">chr\n" + "ACGT" * 100 + "\n>chr second\n" + "ACGT" * 100 + "\n"
    )
    with pytest.raises(CanSNPer2Error, match="reference REF .* more than one sequence named chr"):
        obj.align(q)
    assert sequence_names(str(tmp_path / "refs" / "REF.fna")) == ["chr", "chr"]
    assert duplicate_names(["a", "b"]) is None
    assert check_reference(str(tmp_path / "refs" / "REF.fna"), "REF", None).startswith("reference REF")


def test_a_downloaded_reference_with_duplicate_names_is_reported(tmp_path):
    import gzip

    obj = DownloadGenomes.__new__(DownloadGenomes)
    obj.source, obj.directory, obj.database = "genbank", str(tmp_path / "refs"), None
    os.makedirs(obj.directory)
    data = gzip.compress(b">chr\nACGT\n>chr\nACGT\n")

    def retrieve(link, target, timeout=None):
        with open(target, "wb") as out:
            out.write(data)

    obj.retrieve = retrieve
    results = []

    class Q(list):
        def put(self, item):
            results.append(item)

    obj.download("GCA_000007845.1", "ASM784v1", "Ames", directory=obj.directory, results=Q())
    refid, ok, reason = results[0]
    assert not ok and "more than one sequence named chr" in reason
    assert not os.path.lexists(os.path.join(obj.directory, "Ames.fna"))


def test_an_earlier_result_is_used_only_with_the_same_database_and_settings(tmp_path):
    q = write(tmp_path / "in" / "x.fasta")
    obj, db = project(tmp_path, [q], keep_going=True)
    source = "Query file: /elsewhere/x.fasta sha256:{0}\n".format(obj.file_sha256(q))
    same = obj.settings_line() + "\n"
    assert same.startswith("Typing settings: database sha256:")
    # same database and settings: used again
    result_file(obj, "x", source + same)
    obj.run(database=db)
    assert "x" in obj.reused
    # another database: the sample is refused with the reason
    obj, db = project(tmp_path, [q], keep_going=True)
    other = same.replace("database sha256:", "database sha256:0")
    result_file(obj, "x", source + other)
    obj.run(database=db)
    assert "x" not in obj.reused
    assert "different database or settings" in obj.not_typed["x"]
    # another threshold: refused as well
    obj, db = project(tmp_path, [q], keep_going=True, min_required_hits=2)
    result_file(obj, "x", source + same.replace("min_required_hits=2", "min_required_hits=3"))
    obj.run(database=db)
    assert "different database or settings" in obj.not_typed["x"]
    # another aligner: refused as well, the alignments differ (minimap2 is the default
    # since 2.6.0, so a result from an earlier progressiveMauve run must not pass as current)
    assert "aligner=progressiveMauve" in same
    obj, db = project(tmp_path, [q], keep_going=True, aligner=InstalledMinimap2())
    result_file(obj, "x", source + same)
    obj.run(database=db)
    assert "x" not in obj.reused
    assert "different database or settings" in obj.not_typed["x"]
    # --rerun replaces it
    obj, db = project(tmp_path, [q], keep_going=True, rerun=True)
    result_file(obj, "x", source + other)
    obj.run(database=db)
    assert "different database or settings" not in obj.not_typed.get("x", "")


def test_a_reference_file_placed_by_hand_is_used_without_a_download(tmp_path):
    """ref get with a plain <genome>.fna already in the folder downloaded the genome into
    source/ anyway, kept the file, reported the genome as not fetched and exited 1"""
    obj = DownloadGenomes.__new__(DownloadGenomes)
    obj.source, obj.directory, obj.database = "genbank", str(tmp_path / "refs"), None
    os.makedirs(obj.directory)
    (tmp_path / "refs" / "Ames.fna").write_text(">chr\n" + "ACGT" * 10 + "\n")

    def retrieve(link, target, timeout=None):
        raise AssertionError("downloaded although the file is there")

    obj.retrieve = retrieve
    results = []

    class Q(list):
        def put(self, item):
            results.append(item)

    snps = [("S.1", "chr", 2, "A", "C")]
    obj.download("GCA_000007845.1", "ASM784v1", "Ames", directory=obj.directory, results=Q(), snps=snps)
    assert results == [("Ames", True, "")]
    assert not os.path.exists(os.path.join(obj.directory, "source"))
    # a placed file that does not fit the database is reported, not replaced
    (tmp_path / "refs" / "Ames.fna").write_text(">other\nACGT\n")
    results.clear()
    obj.download("GCA_000007845.1", "ASM784v1", "Ames", directory=obj.directory, results=Q(), snps=snps)
    ((refid, ok, reason),) = results
    assert not ok and "Ames.fna" in reason and "chr" in reason
    assert (tmp_path / "refs" / "Ames.fna").read_text() == ">other\nACGT\n"
