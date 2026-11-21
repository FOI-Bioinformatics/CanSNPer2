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
from CanSNPer2.modules.Aligners import Aligner

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
	kwargs.setdefault("aligner", "minimap2")
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
