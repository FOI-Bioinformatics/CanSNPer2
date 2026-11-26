"""The alignment-only command (runMauve) must leave alignments that CanSNPer2
--skip_alignment can read: progressiveMauve counts reference positions along the joined
reference sequences, so the reference sequence table has to be written next to every
xmfa and copied to the output directory with it.

progressiveMauve is replaced by a small script that writes an xmfa file, so the test
runs without the program."""

import os
import stat
import sys

from CanSNPer2.modules.AlignmentReader import SEQUENCES_SUFFIX, read_sequence_table
from CanSNPer2.modules.runMauve import RunMauve

STUB = """#!{python}
import sys
args = sys.argv[1:]
output = args[args.index("--output") + 1]
with open(output, "w") as out:
	out.write("#FormatVersion Mauve1\\n")
"""


def stub_mauve(folder):
    """an executable that writes the --output file the way progressiveMauve would"""
    path = os.path.join(str(folder), "progressiveMauve")
    with open(path, "w") as out:
        out.write(STUB.format(python=sys.executable))
    os.chmod(path, os.stat(path).st_mode | stat.S_IXUSR)
    return path


def test_run_mauve_writes_the_sequence_table_next_to_each_alignment(tmp_path):
    refdir = tmp_path / "refs"
    refdir.mkdir()
    (refdir / "REF_A.fna").write_text(
        ">chr1 first\n" + "A" * 60 + "\n" + "C" * 40 + "\n>chr2\n" + "G" * 30 + "\n"
    )
    (refdir / "REF_B.fna").write_text(">contig\n" + "T" * 50 + "\n")
    query = tmp_path / "sample.fna"
    query.write_text(">q\nACGT\n")
    outdir = tmp_path / "align"
    runner = RunMauve(
        [str(query)],
        mauve_path=stub_mauve(tmp_path),
        tmpdir=str(tmp_path / "tmp"),
        refdir=str(refdir),
        logdir=str(tmp_path / "logs"),
        outdir=str(outdir),
    )
    aligned = runner.run()
    assert sorted(os.path.basename(f) for f in aligned) == ["REF_A_sample.xmfa", "REF_B_sample.xmfa"]
    assert read_sequence_table(str(outdir / ("REF_A_sample.xmfa" + SEQUENCES_SUFFIX))) == [
        ("chr1", 100),
        ("chr2", 30),
    ]
    assert read_sequence_table(str(outdir / ("REF_B_sample.xmfa" + SEQUENCES_SUFFIX))) == [("contig", 50)]
