"""The alignment step is defined per aligner (CanSNPer2/modules/Aligners.py) and the
alignment files are read by an AlignmentReader; CanSNPer2 itself only runs the commands
and hands the files to the reader. These tests cover the progressiveMauve definition and
run a second, toy aligner through the whole step to show that another aligner needs no
change in CanSNPer2.py."""
import os
import sqlite3
import sys

import pytest

from CanSNPer2.modules.Aligners import ALIGNERS, DEFAULT_ALIGNER, Aligner, ProgressiveMauve, get_aligner
from CanSNPer2.modules.AlignmentReader import AlignmentReader
from CanSNPer2.modules.CanSNPer2 import AlignmentError, CanSNPer2, MauveError
from CanSNPer2.modules.ParseXMFA import ParseXMFA


def make_db(path, genome, rows):
	'''rows: (snp_id, position, ancestral_base, derived_base) against reference genome'''
	conn = sqlite3.connect(path)
	conn.executescript("""
		CREATE TABLE nodes (id integer PRIMARY KEY, name text NOT NULL);
		CREATE TABLE tree (parent integer, child integer, rank_i integer);
		CREATE TABLE snp_annotation (node_id integer PRIMARY KEY, snp_id varchar(6), position integer,
			ancestral_base varchar(1), derived_base varchar(1), reference varchar(20),
			date datetime, genome_i integer);
		CREATE TABLE snp_references (id integer PRIMARY KEY, genome varchar(6), strain varchar(6));
	""")
	conn.execute("INSERT INTO snp_references (id, genome, strain) VALUES (5, ?, ?)", (genome, genome))
	for i, (snp_id, pos, anc, der) in enumerate(rows):
		conn.execute("INSERT INTO nodes (id, name) VALUES (?, ?)", (i + 1, snp_id))
		conn.execute("INSERT INTO snp_annotation (node_id, snp_id, position, ancestral_base, derived_base, reference, date, genome_i)"
			" VALUES (?, ?, ?, ?, ?, 'test', '06/11/26', 5)", (i + 1, snp_id, pos, anc, der))
	conn.commit()
	conn.close()


def test_progressivemauve_is_the_default_aligner():
	assert DEFAULT_ALIGNER == "progressiveMauve"
	assert isinstance(get_aligner(), ProgressiveMauve)
	assert isinstance(get_aligner("progressivemauve"), ProgressiveMauve), "names are matched without case"
	assert set(ALIGNERS) == {"progressiveMauve", "minimap2", "mummer4"}
	with pytest.raises(ValueError, match="Unknown aligner"):
		get_aligner("nosuchaligner")


def test_aligner_path_is_the_program_or_its_folder(tmp_path):
	assert ProgressiveMauve().executable() == "progressiveMauve"
	assert ProgressiveMauve(str(tmp_path)).executable() == str(tmp_path / "progressiveMauve")
	'''the old mauve_path was a prefix ending in "/"'''
	assert ProgressiveMauve("/opt/mauve/").executable() == "/opt/mauve/progressiveMauve"
	assert ProgressiveMauve("/opt/mauve/progressiveMauve-2.4").executable() == "/opt/mauve/progressiveMauve-2.4"


def test_progressivemauve_definition(tmp_path):
	mauve = ProgressiveMauve()
	assert mauve.command("r.fna", "q.fna", "out.xmfa") == ["progressiveMauve", "--output", "out.xmfa", "r.fna", "q.fna"]
	assert mauve.alignment_path("/tmp/run", "FSC200", "sample") == "/tmp/run/FSC200_sample.xmfa"
	assert mauve.log_path("logs", "FSC200", "sample") == "logs/FSC200_sample.mauve.log"
	assert isinstance(mauve.reader(), ParseXMFA)
	assert mauve.exit_status(0) == ("ok", 0)
	assert mauve.exit_status(11) == ("retry", 11)
	assert mauve.exit_status(-6) == ("retry", 6)
	assert mauve.exit_status(1) == ("error", 1)
	query = tmp_path / "q.fna"
	query.write_text(">q\nAC-GT\n")
	retry = mauve.prepare_retry(str(query), str(tmp_path))
	assert open(retry).read() == ">q\nACNGT\n"


def test_mauve_alignment_files_and_reference_names(tmp_path):
	# the reference of each alignment is kept from the reference file, not parsed back out of
	# the alignment file name, which broke for reference names with "_" in them
	refdir = tmp_path / "refs"
	refdir.mkdir()
	(refdir / "Ames_ancestor.fna").write_text(">chr\nACGT\n")
	query = tmp_path / "sample.fna"
	query.write_text(">q\nACGT\n")
	obj = CanSNPer2([str(query)], workdir=str(tmp_path), refdir="refs", tmpdir="shared", logdir="logs", outdir="out", skip_alignment=True)
	files = obj.align(str(query), references=["Ames_ancestor.fna"])
	assert files == [os.path.join(obj.tmpdir, "Ames_ancestor_sample.xmfa")]
	assert obj.alignment_references == ["Ames_ancestor"]
	assert obj.skip_mauve is True and obj.xmfa_files == files, "old attribute names still work"


'''A toy aligner: the "alignment" is a list of "position<TAB>query base" lines for an
ungapped query of the same length as the reference, written by a python one-liner.'''
TOY_SCRIPT = (
	"import sys\n"
	"ref, query, out, code = sys.argv[1:5]\n"
	"seq = lambda p: ''.join(l.strip() for l in open(p) if not l.startswith('>'))\n"
	"open(out, 'w').write(''.join('{0}\\t{1}\\n'.format(i + 1, b) for i, b in enumerate(seq(query))))\n"
	"sys.exit(int(code))\n")


class ToyReader(AlignmentReader):
	def read_alignment(self, path):
		rows = [line.split() for line in open(path)]
		for position in self.snps_in(1, len(rows)):
			pos, base = rows[position - 1]
			self.observe(int(pos), base, int(pos), len(rows), query_name="q", query_pos=pos)


class ToyAligner(Aligner):
	name = "toy"
	program = sys.executable
	suffix = "tsv"

	def __init__(self, exitcode=0):
		super(ToyAligner, self).__init__()
		self.exitcode = exitcode

	def command(self, reference_file, query_file, output):
		return [self.executable(), "-c", TOY_SCRIPT, reference_file, query_file, output, str(self.exitcode)]

	def reader(self, **kwargs):
		return ToyReader(**kwargs)


def toy_case(tmp_path, exitcode=0, **kwargs):
	db = str(tmp_path / "toy.db")
	make_db(db, "REF_A", [("A.1", 3, "G", "T"), ("A.2", 6, "C", "A")])
	refdir = tmp_path / "refs"
	refdir.mkdir()
	(refdir / "REF_A.fna").write_text(">chr\nACGTACGT\n")
	query = tmp_path / "sample.fna"
	query.write_text(">q\nACTTACGT\n")
	obj = CanSNPer2([str(query)], workdir=str(tmp_path), refdir="refs", tmpdir="shared", logdir="logs", outdir="out",
		database=db, aligner=ToyAligner(exitcode), **kwargs)
	return obj, db, query


def test_another_aligner_runs_through_the_whole_alignment_step(tmp_path):
	obj, db, query = toy_case(tmp_path)
	files = obj.align(str(query))
	assert files == [os.path.join(obj.tmpdir, "REF_A_sample.tsv")]
	assert os.path.isfile(os.path.join(obj.logdir, "REF_A_sample.toy.log"))
	assert not [f for f in os.listdir(obj.tmpdir) if f.startswith(("ref_", "query_"))], "inputs are only staged for aligners that need it"
	snps, snp_info, called = obj.find_snps_multiproc(xmfa_obj=obj.aligner.reader(database=db), xmfa_files=files,
		references=obj.alignment_references)
	assert snps == {"A.1": 1, "A.2": 2}
	assert called == ["A.1"]
	assert [row[:5] + row[7:] for row in snp_info] == [
		["A.1", "REF_A", "3", "q", "3", "T", "edge"],
		["A.2", "REF_A", "6", "q", "6", "C", "edge"]]


def test_a_failing_aligner_stops_the_run_or_skips_the_sample(tmp_path):
	obj, db, query = toy_case(tmp_path, exitcode=3)
	with pytest.raises(AlignmentError, match="toy processes did not run correctly"):
		obj.align(str(query))
	assert MauveError is AlignmentError, "the old exception name is still raised"
	keep = tmp_path / "keep"
	keep.mkdir()
	obj, db, query = toy_case(keep, exitcode=3, keep_going=True)
	assert obj.align(str(query)) == []
