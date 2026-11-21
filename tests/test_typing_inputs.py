"""Checks of the typing command before and around a run: the printed call, the query list,
the aligner and the output folders, simultaneous creation of folders, the SNP path of a
result file and the query fasta."""
import gzip
import os
import sqlite3
import threading

import pytest

import helpers
from CanSNPer2.modules.CanSNPer2 import CanSNPer2, CanSNPer2Error
from CanSNPer2.modules.Aligners import Aligner, Minimap2


class Installed(Minimap2):
	'''minimap2 reported as found, so that runs reach the steps after the aligner check'''
	def available(self):
		return True


GENOME = ">chr1 test\n" + "ACGT" * 250 + "\n"


def make_db(path):
	"""root -> A.1 -> (B.1, B.2) and root -> C.1"""
	conn = sqlite3.connect(path)
	conn.executescript(helpers.SNP_ANNOTATION_TABLE + helpers.REFERENCES_TABLE + """
		CREATE TABLE nodes (id integer PRIMARY KEY, name text NOT NULL);
		CREATE TABLE tree (parent integer, child integer, rank_i integer, unique (parent, child, rank_i));
		INSERT INTO snp_references (id, genome, strain) VALUES (5, 'REF', 'ref strain');
		INSERT INTO nodes (id, name) VALUES (1, 'root'), (2, 'A.1'), (3, 'B.1'), (4, 'B.2'), (5, 'C.1');
		INSERT INTO tree (parent, child, rank_i) VALUES (1, 1, 1), (1, 2, 1), (2, 3, 1), (2, 4, 1), (1, 5, 1);
	""")
	for node_id, snp_id, pos in ((2, "A.1", 10), (3, "B.1", 20), (4, "B.2", 30), (5, "C.1", 40)):
		conn.execute("INSERT INTO snp_annotation (node_id, snp_id, sequence, sequence_position, ancestral_base, derived_base,"
			" reference, date, genome_i) VALUES (?, ?, 'chr1', ?, 'C', 'T', 'REF', '06/11/26', 5)", (node_id, snp_id, pos))
	conn.commit()
	conn.close()


def project(tmp_path, queries, **kwargs):
	db = str(tmp_path / "t.db")
	if not os.path.exists(db):
		make_db(db)
	refdir = tmp_path / "refs"
	refdir.mkdir(exist_ok=True)
	(refdir / "REF.fna").write_text(GENOME)
	kwargs.setdefault("aligner", Installed())
	kwargs.setdefault("export", True)
	return CanSNPer2(queries, workdir=str(tmp_path), refdir="refs", tmpdir="shared", logdir="logs", outdir="out",
		database=db, **kwargs), db


def write(path, text=GENOME, mode="w"):
	path.parent.mkdir(parents=True, exist_ok=True)
	with open(path, mode) as out:
		out.write(text)
	return str(path)


def fake_alignment(obj, called, states, query):
	"""replace the aligning and reading steps: the called SNPs and states are given"""
	obj.align = lambda q, references=[]: ["REF_x.paf"]
	rows = [[name, "REF", "chr1", "1", "chr1", "1", "C", "T", "T", ""] for name in called]
	obj.find_snps_multiproc = lambda **kwargs: (states, rows, list(called))


# 1. the call is printed without result files

def test_the_call_is_printed_with_no_snpfiles(tmp_path, capsys):
	q = write(tmp_path / "x.fasta")
	obj, db = project(tmp_path, [q], export=False)
	fake_alignment(obj, ["A.1", "B.1", "B.2"], {"A.1": 1, "B.1": 1, "B.2": 1}, q)
	obj.min_required_hits = 2
	obj.run(database=db)
	assert "x: " in capsys.readouterr().out
	assert not os.path.exists(os.path.join(obj.outdir, "x_snps.txt"))


# 2. query list detection

def test_a_list_file_is_read_by_content(tmp_path):
	a = write(tmp_path / "a.fasta")
	lst = write(tmp_path / "list.txt", a + "\n")
	obj, db = project(tmp_path, [lst])
	assert obj.resolve_queries([lst]) == [a]


def test_a_fasta_named_txt_is_a_fasta(tmp_path):
	q = write(tmp_path / "x.txt")
	obj, db = project(tmp_path, [q])
	assert obj.resolve_queries([q]) == [q]


def test_a_list_with_other_arguments_is_refused(tmp_path):
	a = write(tmp_path / "a.fasta")
	lst = write(tmp_path / "list.txt", a + "\n")
	obj, db = project(tmp_path, [a, lst])
	with pytest.raises(CanSNPer2Error, match="list.txt.*only query argument"):
		obj.resolve_queries([a, lst])
	with pytest.raises(CanSNPer2Error, match="only query argument"):
		obj.resolve_queries([lst, a])


def test_an_empty_file_is_a_query_not_a_list(tmp_path):
	# it is refused as a FASTA file, with a message that says it is empty (see the query fasta checks)
	empty = write(tmp_path / "list.txt", "\n\n")
	obj, db = project(tmp_path, [empty])
	assert obj.resolve_queries([empty]) == [empty]


# 3. the aligner

class Missing(Aligner):
	name = "missing"
	program = "no_such_aligner_program"


def test_a_missing_aligner_stops_before_the_first_sample(tmp_path):
	q = write(tmp_path / "x.fasta")
	obj, db = project(tmp_path, [q], aligner=Missing())
	with pytest.raises(CanSNPer2Error, match="missing.*not found.*--aligner_path.*CanSNPer2-test"):
		obj.run(database=db)
	assert not obj.not_typed


def test_skip_alignment_needs_no_aligner(tmp_path):
	q = write(tmp_path / "x.fasta")
	obj, db = project(tmp_path, [q], aligner=Missing(), skip_alignment=True, keep_going=True)
	obj.check_aligner()  # does not raise


# 4. output folders

def test_an_outdir_that_is_a_file_is_refused(tmp_path):
	write(tmp_path / "out", "x")
	with pytest.raises(CanSNPer2Error, match="--outdir .*out.* is not a directory"):
		project(tmp_path, [])


@pytest.mark.skipif(os.geteuid() == 0, reason="root can write anywhere")
def test_a_read_only_outdir_is_refused(tmp_path):
	out = tmp_path / "out"
	out.mkdir()
	out.chmod(0o500)
	try:
		with pytest.raises(CanSNPer2Error, match="--outdir .*not writable"):
			project(tmp_path, [])
	finally:
		out.chmod(0o700)


def test_a_logdir_that_is_a_file_is_refused(tmp_path):
	write(tmp_path / "logs", "x")
	with pytest.raises(CanSNPer2Error, match="--logdir .*not a directory"):
		project(tmp_path, [])


def test_a_tmpdir_that_is_a_file_is_refused(tmp_path):
	write(tmp_path / "shared", "x")
	with pytest.raises(CanSNPer2Error, match="--tmpdir .*not a directory"):
		project(tmp_path, [])


# 5. folders created by two runs at once

def test_folders_that_appear_during_creation_are_not_an_error(tmp_path, monkeypatch):
	real = os.path.exists
	made = []
	def late(path):
		# another run creates the folder just after this one saw it missing
		if path.endswith(("out", "logs", "shared")) and path not in made:
			made.append(path)
			os.makedirs(path, exist_ok=True)
			return False
		return real(path)
	monkeypatch.setattr(os.path, "exists", late)
	project(tmp_path, [])
	assert len(made) == 3 or len(made) >= 3
