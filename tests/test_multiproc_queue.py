"""Regression tests for issue #41: SNP calling with a large custom database must
not hang the run.

Two directions of the queue handling are checked:
  * the tens of thousands of items that a >2000 SNP database puts on the output
    queues must all arrive, no matter how many alignments run in parallel;
  * a worker that stops delivering must be stopped and reported, never joined
    forever (that infinite join is what made the whole command freeze)."""

import os
import signal
import time

import pytest

from helpers import make_db
from CanSNPer2.modules.CanSNPer2 import CanSNPer2, CanSNPer2Error
from CanSNPer2.modules.AlignmentReader import SEQUENCES_SUFFIX
from CanSNPer2.modules.ParsePAF import ParsePAF
from CanSNPer2.modules.ParseXMFA import ParseXMFA


def make_big_db(path, nsnp):
    """a database with nsnp custom SNPs registered against one reference (REF), all on its
    sequence chr"""
    positions = list(range(10, 10 + 40 * nsnp, 40))
    make_db(path, [("C.{0}".format(n + 1), "chr", pos, "A", "T") for n, pos in enumerate(positions)])
    return positions


def make_paf(path, positions):
    """one minimap2 alignment covering every registered SNP, the query carries the derived base;
    the reader of a PAF file is the AlignmentReader the queue handling is tested with"""
    last = positions[-1] + 10
    cs, at = [], 1
    for pos in positions:
        if pos > at:
            cs.append("=" + "C" * (pos - at))
        cs.append("*at")
        at = pos + 1
    cs.append("=" + "C" * (last - at + 1))
    fields = [
        "query",
        last,
        0,
        last,
        "+",
        "chr",
        last,
        0,
        last,
        last,
        last,
        60,
        "tp:A:P",
        "cs:Z:" + "".join(cs),
    ]
    with open(path, "w") as handle:
        handle.write("\t".join(str(f) for f in fields) + "\n")


def make_xmfa(path, positions):
    """one progressiveMauve block covering every registered SNP, the query carries the derived
    base; the reference sequence table that ProgressiveMauve.prepare writes is written next
    to it"""
    last = positions[-1] + 10
    ref = ["C"] * last
    target = list(ref)
    for pos in positions:
        ref[pos - 1] = "A"
        target[pos - 1] = "T"
    with open(path, "w") as handle:
        handle.write("#FormatVersion Mauve1\n============\n")
        for seq in (ref, target):
            handle.write(
                "> 1:1-{0} + seq\n{1}\n".format(
                    last, "\n".join("".join(seq[i : i + 60]) for i in range(0, last, 60))
                )
            )
    with open(path + SEQUENCES_SUFFIX, "w") as handle:
        handle.write("chr\t{0}\n".format(last))


## the readers the queue handling is tested with: (reader class, alignment writer, file ending)
READERS = [
    pytest.param(ParsePAF, make_paf, "paf", id="minimap2"),
    pytest.param(ParseXMFA, make_xmfa, "xmfa", id="progressiveMauve"),
]


@pytest.mark.parametrize("reader, write_alignment, suffix", READERS)
def test_large_custom_database_is_collected_completely(tmp_path, reader, write_alignment, suffix):
    """>2000 SNPs per reference over four alignments: every queue item must be taken,
    the queues may not fill up and drop workers (issue #41)"""
    nsnp, nalign = 2500, 4
    db = str(tmp_path / "big.db")
    positions = make_big_db(db, nsnp)
    alignments = []
    for n in range(nalign):
        alignment = str(tmp_path / "REF_big{0}.{1}".format(n, suffix))
        write_alignment(alignment, positions)
        alignments.append(alignment)

    obj = CanSNPer2(["query.fna"], workdir=str(tmp_path), tmpdir="tmp", logdir="logs", outdir="out")
    snps, snp_info, called = obj.find_snps_multiproc(xmfa_obj=reader(database=db), xmfa_files=alignments)
    assert len(snps) == nsnp
    assert len(snp_info) == nsnp * nalign
    assert len(called) == nsnp * nalign


class StallingAlignment(object):
    """one alignment delivers its SNPs, the other one never comes back, which is what
    a worker stuck inside the database or blocked on a full queue looks like"""

    def __init__(self, stop_file):
        self.stop_file = stop_file
        self.SNPS = {}
        self.SNP_info = []
        self.called_snps = []
        self.reference = "REF"

    def run(self, xmfa_file, *args, **kwargs):
        self.reference = os.path.basename(xmfa_file).split("_")[0]
        if "stall" in os.path.basename(xmfa_file):
            deadline = time.time() + 120
            while not os.path.exists(self.stop_file) and time.time() < deadline:
                time.sleep(0.2)
            raise RuntimeError("the stalled worker was released")
        self.SNPS = {"C.1": 1}
        self.SNP_info = [["C.1", self.reference, "10", "A", "T", "T"]]
        self.called_snps = ["C.1"]

    def get_snps(self):
        return self.SNPS

    def get_snp_info(self):
        return self.SNP_info

    def get_called_snps(self):
        return self.called_snps


def test_stalled_worker_is_stopped_and_reported(tmp_path):
    """a worker that stops delivering must not be waited for forever: the run ends as a
    reported failure of that sample instead of freezing the whole command (issue #41)"""
    stop_file = str(tmp_path / "release_stalled_worker")
    alignments = []
    for name in ("REF_ok.xmfa", "REF_stall.xmfa"):
        alignment = str(tmp_path / name)
        open(alignment, "w").write("#FormatVersion Mauve1\n")
        alignments.append(alignment)

    obj = CanSNPer2(
        ["query.fna"],
        workdir=str(tmp_path),
        tmpdir="tmp",
        logdir="logs",
        outdir="out",
        queue_idle_timeout=1,
        worker_join_timeout=1,
    )

    def watchdog(signum, frame):
        raise AssertionError(
            "find_snps_multiproc never returned: the stalled worker was joined forever (issue #41)"
        )

    old_handler = signal.signal(signal.SIGALRM, watchdog)
    signal.setitimer(signal.ITIMER_REAL, 60)
    start = time.time()
    try:
        with pytest.raises(CanSNPer2Error) as error:
            obj.find_snps_multiproc(xmfa_obj=StallingAlignment(stop_file), xmfa_files=alignments)
        elapsed = time.time() - start
        assert elapsed < 30, "the stalled worker was waited for {0:.0f} seconds".format(elapsed)
        assert "Only 1 of 2" in str(error.value)
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, old_handler)
        with open(stop_file, "w") as handle:  ## let a stalled child that is still around leave
            handle.write("")


class HalfBrokenAlignment(object):
    """one alignment parses, the other one stops half way through and leaves its partial
    result behind, which is what a database or file error looks like from the queue side"""

    def __init__(self):
        self.SNPS = {}
        self.SNP_info = []
        self.called_snps = []
        self.reference = "REF"

    def run(self, xmfa_file, *args, **kwargs):
        self.reference = os.path.basename(xmfa_file).split("_")[0]
        self.SNPS = {"C.1": 1}
        self.SNP_info = [["C.1", self.reference, "10", "A", "T", "T"]]
        self.called_snps = ["C.1"]
        if "broken" in os.path.basename(xmfa_file):
            raise ValueError("the alignment could not be read to its end")

    def get_snps(self):
        return self.SNPS

    def get_snp_info(self):
        return self.SNP_info

    def get_called_snps(self):
        return self.called_snps


def test_failed_alignment_is_not_typed_from(tmp_path):
    """a worker whose parse raised used to send ordinary finish signals with its partial
    results, and the sample was typed from a half-read alignment (issue #50)"""
    alignments = []
    for name in ("REF_ok.xmfa", "REF_broken.xmfa"):
        alignment = str(tmp_path / name)
        open(alignment, "w").write("#FormatVersion Mauve1\n")
        alignments.append(alignment)

    obj = CanSNPer2(["query.fna"], workdir=str(tmp_path), tmpdir="tmp", logdir="logs", outdir="out")
    with pytest.raises(CanSNPer2Error) as error:
        obj.find_snps_multiproc(xmfa_obj=HalfBrokenAlignment(), xmfa_files=alignments)
    message = str(error.value)
    assert "failed during SNP calling" in message
    assert "alignment could not be read to its end" in message


@pytest.mark.parametrize("reader, write_alignment, suffix", READERS)
def test_worker_opens_its_own_database_connection(tmp_path, reader, write_alignment, suffix):
    """an AlignmentReader carries the database path and connects inside run(), so the parent
    never holds a cursor that its forked workers would read through (issue #49)"""
    db = str(tmp_path / "one.db")
    positions = make_big_db(db, 3)
    alignment = str(tmp_path / "REF_one.{0}".format(suffix))
    write_alignment(alignment, positions)

    parser = reader(database=db)
    assert parser.database is False, "the parent connected to the database it hands to the workers"
    obj = CanSNPer2(["query.fna"], workdir=str(tmp_path), tmpdir="tmp", logdir="logs", outdir="out")
    snps, snp_info, called = obj.find_snps_multiproc(xmfa_obj=parser, xmfa_files=[alignment])
    assert len(snps) == 3 and len(called) == 3
    assert parser.database is False, "only the worker may connect, the parent object stays unconnected"
