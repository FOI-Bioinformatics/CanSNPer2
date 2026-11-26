"""Follow-up to the 2.3.0 stress test: database connections of the trees, an overall download
time limit, a FASTA file without a header, the aligner check for a run that only reuses results,
the unused ReadTaxonomy module, and a SNP on an unknown genome in the loader."""

import gc
import os
import time

import pytest

from CanSNPer2 import DownloadGenomes as dg
from CanSNPer2.modules.CanSNPer2 import CanSNPer2Error
from test_aligners import ToyAligner, toy_case


def open_fds():
    return len(os.listdir("/dev/fd"))


def test_the_tree_of_each_sample_closes_its_database_connection(tmp_path):
    obj, db, query = toy_case(tmp_path)
    gc.collect()
    gc.disable()
    try:
        before = open_fds()
        for i in range(40):
            obj.create_tree([], "s{0}".format(i), [], False, min_required_hits=1)
        obj.add_sample_report("s0", "clean", "A.1", "A.1", [])
        obj.write_dashboard()
        after = open_fds()
    finally:
        gc.enable()
    assert after - before < 3, (before, after)


class Dribble(object):
    """a response that never ends and sends one byte every 0.05 s, below any per-read timeout"""

    def __enter__(self):
        return self

    def __exit__(self, *args):
        return False

    def read(self, n=-1):
        time.sleep(0.05)
        return b"A"


def test_a_download_that_keeps_sending_slowly_is_stopped(tmp_path, monkeypatch):
    monkeypatch.setattr(dg, "urlopen", lambda link, timeout=None: Dribble())
    obj = dg.DownloadGenomes.__new__(dg.DownloadGenomes)
    start = time.time()
    with pytest.raises(TimeoutError, match="did not finish within"):
        obj.retrieve("https://example.invalid/x.fna.gz", str(tmp_path / "x"), timeout=5, deadline=0.5)
    assert time.time() - start < 3


def sleeper():
    time.sleep(30)


def test_a_worker_that_runs_too_long_is_stopped():
    p = dg._MP_CTX.Process(target=sleeper)
    p.start()
    start = time.time()
    late = dg.join_workers([p], deadline=0.5)
    assert late == [p] and not p.is_alive()
    assert time.time() - start < 5


def test_a_fasta_file_without_a_header_is_named_as_such(tmp_path):
    obj, db, query = toy_case(tmp_path)
    nohdr = tmp_path / "nohdr.fasta"
    nohdr.write_text("ACGTACGTACGTNNACGT\n" * 500)
    obj.query = [str(nohdr)]
    with pytest.raises(CanSNPer2Error, match="no FASTA header") as info:
        obj.run(database=db)
    assert len(str(info.value)) < 400


def test_a_long_list_of_missing_queries_is_shortened(tmp_path):
    obj, db, query = toy_case(tmp_path)
    obj.query = [str(tmp_path / "missing{0}.fna".format(i)) for i in range(100)]
    with pytest.raises(CanSNPer2Error) as info:
        obj.run(database=db)
    assert len(str(info.value)) < 1000 and "100" in str(info.value)


class Missing(ToyAligner):
    def available(self):
        return False


def test_a_run_that_only_reuses_results_needs_no_aligner(tmp_path):
    obj, db, query = toy_case(tmp_path, min_required_hits=1, export=True)
    obj.run(database=db)
    assert "sample" in obj.typed
    from CanSNPer2.modules.CanSNPer2 import CanSNPer2

    again = CanSNPer2(
        [str(query)],
        workdir=str(tmp_path),
        refdir="refs",
        tmpdir="shared",
        logdir="logs",
        outdir="out",
        database=db,
        aligner=Missing(),
        min_required_hits=1,
        export=True,
    )
    again.run(database=db)
    assert "sample" in again.reused
    rerun = CanSNPer2(
        [str(query)],
        workdir=str(tmp_path),
        refdir="refs",
        tmpdir="shared",
        logdir="logs",
        outdir="out",
        database=db,
        aligner=Missing(),
        min_required_hits=1,
        rerun=True,
        export=True,
    )
    with pytest.raises(CanSNPer2Error, match="was not found"):
        rerun.run(database=db)


def test_the_unused_flextaxd_reader_is_gone():
    import CanSNPer2.modules as modules

    assert not os.path.exists(os.path.join(os.path.dirname(modules.__file__), "ReadTaxonomy.py"))


def test_the_loader_refuses_a_snp_on_an_unknown_genome(tmp_path):
    import helpers
    from CanSNPer2.modules.DatabaseConnection import SourceFileError
    from CanSNPer2.modules.InitiateDatabase import CanSNPDatabase

    paths = helpers.write_sources(tmp_path / "src", ["A.1\t1\tPaper\tREF1\tCP000001.1\t100\tT\tC\t\tyes"])
    snp = open(paths["snp"]).read().replace("\tREF1\t", "\tOTHER\t")
    open(paths["snp"], "w").write(snp)
    database = CanSNPDatabase.__new__(CanSNPDatabase)
    import sqlite3

    db = str(tmp_path / "t.db")
    helpers.make_db(db, [("A.1", "CP000001.1", 100, "C", "T")], genome="REF1")
    database.conn = sqlite3.connect(db)
    database.genomes = {}
    database.get_nodes = lambda: {"A.1": 1}
    with pytest.raises(SourceFileError, match="OTHER"):
        database.load_annotation_file(paths["snp"])
