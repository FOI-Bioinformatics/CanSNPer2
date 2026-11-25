"""Regression tests for the database build and the input side of the pipeline.

#34/#32: ReadTaxonomyCanSNPer used to inherit a node named "root" plus a self-link from
ReadTaxonomyBase.__init__, so every created database had an orphan node "root" above its
own tree root and the run warned "#Node has no snp_annotation root" for it. The base node
creation now lives in parse_taxonomy, which the CanSNPer reader overrides.
#31: a tree file with ";" instead of tabs builds the same database, no warning flood.
#22: # comment lines and a trailing free-text comment column in the snp annotation input.
#26: characters that are not ACGTN are replaced with N before anything is aligned.
#12: a progressiveMauve warning names the genome and the reference it failed for.
#33: a reference folder may hold more genomes than the database uses."""
import os
import sqlite3
import sys

import pytest

import helpers
from CanSNPer2.modules.InitiateDatabase import CanSNPDatabase, CreateSNPDatabase
from CanSNPer2.modules.CanSNPer2 import CanSNPer2, CanSNPer2Error
from CanSNPer2.modules.DatabaseConnection import CanSNPdbFunctions


REFERENCES_FILE = ("genome\tstrain\tgenbank_id\trefseq_id\tassembly_name\n"
                   "REF1\tref one strain\tGCA_000000001.1\tGCF_000000001.1\tasm1\n")

def tree_file(separator):
	return "R\nR{0}A.1\nR{0}A.1{0}B.1\n".format(separator)

def annotation_file(extra_comment_column=True):
	head = "\t".join(helpers.CANSNP_HEADER)
	lines = [head + ("\tcomment" if extra_comment_column else "")]
	if extra_comment_column:
		lines.append("A.1\t1\tBirdsell2014\tREF1\tCP000001.1\t100\tT\tC\t\tyes\ta named branch")
	else:
		lines.append("A.1\t1\tBirdsell2014\tREF1\tCP000001.1\t100\tT\tC\t\tyes")
	return "\n".join(lines) + "\n"


def build_db(tmp_path, name, tree, annotation, references=REFERENCES_FILE):
	db = str(tmp_path / name)
	header, _, body = annotation.partition("\n")
	rows = body.split("\n")
	if rows and rows[-1] == "":
		rows.pop()
	paths = helpers.write_sources(tmp_path, rows, tree=tree, references=references,
		name=name.rsplit(".", 1)[0], header=header)
	CanSNPDatabase(database=db, create=True, tree=paths["tree"], annotation=paths["snp"],
	               references=paths["references"], source_type="CanSNPer", verbose=False)
	return db


def tables(db):
	conn = sqlite3.connect(db)
	nodes = conn.execute("SELECT id,name FROM nodes ORDER BY id").fetchall()
	links = conn.execute("SELECT parent,child FROM tree ORDER BY parent,child").fetchall()
	ann = conn.execute("SELECT snp_id,sequence,sequence_position FROM snp_annotation ORDER BY snp_id").fetchall()
	conn.close()
	return nodes, links, ann


def test_created_database_has_one_root_not_an_extra_named_root(tmp_path):
	# #34, #32: the tree root from the file is node 1 with the only self-link; there is no
	# second node named "root" sitting above the tree with no snp_annotation
	db = build_db(tmp_path, "rootfix.db", tree_file("\t"), annotation_file())
	nodes, links, ann = tables(db)
	assert nodes == [(1, "R"), (2, "A.1"), (3, "B.1")], nodes
	assert "root" not in [name for _, name in nodes], "no phantom root node"
	assert (1, 1) in links and [lc for lc in links if lc[0] == lc[1]] == [(1, 1)]
	assert ann == [("A.1", "CP000001.1", 100)]


def test_semicolon_separated_tree_builds_the_same_database(tmp_path, caplog):
	# #31: ";" instead of tabs used to render "#Node has no snp_annotation" warnings for
	# every node; the CanSNPer reader accepts both, and both builds must be identical
	tab_db = build_db(tmp_path, "tab.db", tree_file("\t"), annotation_file())
	semi_db = build_db(tmp_path, "semi.db", tree_file(";"), annotation_file())
	assert tables(tab_db) == tables(semi_db)
	# and the tree reads back without a flood of annotation warnings: only the tree root
	# itself (no snp_annotation by design) may warn
	from CanSNPer2.modules.NewickTree import NewickTree
	nt = NewickTree(semi_db, outdir=str(tmp_path))
	assert nt.nodeDict["root"].name == "-R-"  # the unannotated tree root stays a label


def test_annotation_file_comments_and_trailing_comment_column(tmp_path, caplog):
	# #22: # lines are comments, blank lines are ignored, and a row may leave out the
	# trailing free-text comment column; a row short in a real column is refused
	# (tests/test_create_input_checks.py)
	ann = ("\t".join(helpers.CANSNP_HEADER) + "\tcomment\n"
	       "#A.1 is one of the major branches\n"
	       "\n"
	       "A.1\t1\tBirdsell2014\tREF1\tCP000001.1\t100\tT\tC\t\tyes\ta named branch\n"
	       "B.1\t1\tVogler2009\tREF1\tCP000001.1\t200\tT\tG\t\tyes\n")
	db = build_db(tmp_path, "comment.db", tree_file("\t"), ann)
	nodes, links, rows = tables(db)
	assert [r[0] for r in rows] == ["A.1", "B.1"], rows


def clean_case(obj, tmp_path, content):
	src = tmp_path / "genome.fna"
	src.write_text(content)
	cleaned = obj.clean_input_fasta(str(src), "query")
	return src, cleaned


def test_ambiguous_nucleotides_and_gaps_are_replaced_before_alignment(tmp_path):
	# #26 (Y, M, R ...) and the unseen dashes of #28: everything that is not ACGT becomes
	# N before mauve ever sees the file, position by position
	obj = CanSNPer2([], workdir=str(tmp_path))
	src, cleaned = clean_case(obj, tmp_path,
		">seq one\nACGTacgt-.*RYSMWK\nACGTN\n")
	text = open(cleaned).read()
	assert cleaned != str(src)
	assert text == ">seq one\nACGTacgtNNNNNNNNN\nACGTN\n", repr(text)
	# an already clean file is used as it is, no copy written
	src2, cleaned2 = clean_case(obj, tmp_path, ">ok\nACGTNacgt\n")
	assert cleaned2 == str(src2)


def test_mauve_warnings_name_genome_and_reference(tmp_path, caplog):
	# #12: five unnamed "exitcode: 11" warnings never told which genome failed against
	# which reference; the names come from the staged mauve command arguments, and they
	# must stay attached to the right alignment while the process containers shrink
	import logging
	obj = CanSNPer2([], workdir=str(tmp_path), aligner="progressiveMauve")
	logdir = tmp_path / "mauve_logs"
	logdir.mkdir()
	fast = [sys.executable, "-c", "import sys;sys.exit(11)", "ref_OSU18.fna", "query_genfast.fna"]
	slow = [sys.executable, "-c", "import time;time.sleep(1.5);import sys;sys.exit(11)",
	        "ref_LVS.fna", "query_genslow.fna"]
	ret = obj.run_mauve([fast, slow], [str(logdir / "a.log"), str(logdir / "b.log")])
	assert ret == 11
	warnings = [rec.getMessage() for rec in caplog.records if rec.levelno == logging.WARNING]
	assert any("genfast" in m and "OSU18" in m for m in warnings), warnings
	assert any("genslow" in m and "LVS" in m for m in warnings), warnings


def test_extra_references_in_the_folder_are_tolerated(tmp_path):
	# #33: the reference folder may hold genomes the database does not use, they are
	# filtered out and the run proceeds with the ones the database knows
	db = str(tmp_path / "refs.db")
	CreateSNPDatabase(db, create=True)
	conn = CanSNPdbFunctions(db)
	for node_id, genome in enumerate(("OSU18", "LVS"), start=1):
		genome_i = conn.insert({"genome": genome, "strain": genome}, "snp_references")
		## only reference genomes with SNPs are aligned against
		conn.insert({"node_id": node_id, "snp_id": "S.{0}".format(node_id), "sequence": genome, "sequence_position": 1,
			"ancestral_base": "A", "derived_base": "T", "genome_i": genome_i}, "snp_annotation")
	conn.commit()
	conn.disconnect()
	refdir = tmp_path / "refdir"
	refdir.mkdir()
	for name in ("OSU18", "LVS", "NOTINTHEDATABASE"):
		(refdir / (name + ".fna")).write_text(">{}\nACGT\n".format(name))
	obj = CanSNPer2([], database=db, refdir=str(refdir), workdir=str(tmp_path))
	import logging
	logging.disable(logging.CRITICAL)
	try:
		selected = obj.get_references()
	finally:
		logging.disable(logging.NOTSET)
	assert selected == ["LVS.fna", "OSU18.fna"], selected


# sequence and position on that sequence ------------------------------------------------

def test_get_snps_keys_by_sequence(tmp_path):
	# two SNPs at position 50, one on each chromosome, are two different sites
	from CanSNPer2.modules.DatabaseConnection import XMFAFunctions
	db = helpers.make_db(tmp_path / "two.db", [("A.1", "chrI", 50, "C", "T"), ("A.2", "chrII", 50, "G", "A"),
		("A.3", "chrII", 20, "A", "C")])
	conn = XMFAFunctions(db)
	snps, positions = conn.get_snps("REF")
	conn.disconnect()
	assert snps == {("chrI", 50): ("C", "T", "A.1"), ("chrII", 50): ("G", "A", "A.2"), ("chrII", 20): ("A", "C", "A.3")}
	assert positions == {"chrI": [50], "chrII": [20, 50]}


def test_get_snps_refuses_the_same_site_twice(tmp_path):
	# the read-time guard in XMFAFunctions.get_results is now unreachable through a database the loaders built
	with pytest.raises(sqlite3.IntegrityError):
		helpers.make_db(tmp_path / "dup.db", [("A.1", "chrII", 50, "C", "T"), ("A.2", "chrII", 50, "G", "A")])


def test_joined_database_is_refused(tmp_path):
	from CanSNPer2.modules.DatabaseConnection import check_database, ConnectionError
	db = helpers.make_db(tmp_path / "joined.db", [("A.1", "chrI", 50, "C", "T")], joined=True)
	with pytest.raises(ConnectionError) as error:
		check_database(db)
	assert "db convert" in str(error.value)
	assert check_database(db, current_format=False)
	assert check_database(helpers.make_db(tmp_path / "current.db", [("A.1", "chrI", 50, "C", "T")]))


@pytest.mark.parametrize("module,arguments", [
	("CanSNPer2.CanSNPerTree", ["query.fna", "--logdir", "logs"]),
	("CanSNPer2.DownloadGenomes", ["--logs", "logs"]),
	("CanSNPer2.SNPDatabase", ["--export", "--logs", "logs"]),
])
def test_cli_refuses_joined_database(tmp_path, module, arguments):
	import subprocess
	db = helpers.make_db(tmp_path / "joined.db", [("A.1", "chrI", 50, "C", "T")], joined=True)
	env = dict(os.environ)
	root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
	env["PYTHONPATH"] = root + os.pathsep + env.get("PYTHONPATH", "")
	result = subprocess.run([sys.executable, "-m", module, "--database", db] + arguments,
		capture_output=True, text=True, cwd=str(tmp_path), env=env)
	assert result.returncode != 0, result.stdout + result.stderr
	assert "db convert" in result.stderr, (result.stdout, result.stderr)


def run_database_cli(tmp_path, arguments, code=None):
	'''run CanSNPer2-database in a subprocess from tmp_path; code, when given, is Python run
		before main() in the same process, to make a part of the build fail'''
	import subprocess
	env = dict(os.environ)
	root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
	env["PYTHONPATH"] = root + os.pathsep + env.get("PYTHONPATH", "")
	if code is None:
		command = [sys.executable, "-m", "CanSNPer2.SNPDatabase"]
	else:
		command = [sys.executable, "-c", code + "\nimport sys\nsys.argv[0] = 'CanSNPer2-database'\n"
			"from CanSNPer2.SNPDatabase import main\nmain()\n"]
	return subprocess.run(command + arguments + ["--logs", "logs"],
		capture_output=True, text=True, cwd=str(tmp_path), env=env)


def create_arguments(tmp_path, db):
	paths = helpers.write_sources(tmp_path, ["A.1\t1\tBirdsell2014\tREF1\tCP000001.1\t100\tT\tC\t\tyes",
		"B.1\t1\tBirdsell2014\tREF1\tCP000001.1\t200\tG\tA\t\tyes"])
	return ["--create", "--database", db, "--tree", paths["tree"], "--annotation", paths["snp"],
		"--references", paths["references"]]


def row_counts(db):
	conn = sqlite3.connect(db)
	try:
		return tuple(conn.execute("SELECT count(*) FROM {0}".format(table)).fetchone()[0]
			for table in ("nodes", "tree", "snp_annotation", "snp_references"))
	finally:
		conn.close()


def test_create_refuses_an_existing_joined_database(tmp_path):
	'''--create into a database from before 2.3.0 added the tree and the references a second
		time, failed on every SNP insert with a warning only, and ended with exit status 0'''
	db = helpers.make_db(tmp_path / "joined.db", [("A.1", "chrI", 50, "C", "T")], joined=True)
	before = row_counts(db)
	result = run_database_cli(tmp_path, create_arguments(tmp_path, db))
	assert result.returncode == 1, result.stdout + result.stderr
	assert "db convert" in result.stderr, (result.stdout, result.stderr)
	assert row_counts(db) == before


def test_create_refuses_an_existing_database(tmp_path):
	'''a second --create into a current-format database doubled its nodes, references and SNPs'''
	db = str(tmp_path / "scheme.db")
	arguments = create_arguments(tmp_path, db)
	first = run_database_cli(tmp_path, arguments)
	assert first.returncode == 0, first.stdout + first.stderr
	before = row_counts(db)
	second = run_database_cli(tmp_path, arguments)
	assert second.returncode == 1, second.stdout + second.stderr
	assert "already exists" in second.stderr, (second.stdout, second.stderr)
	assert row_counts(db) == before


def test_create_fails_when_a_snp_cannot_be_written(tmp_path):
	'''a failed SNP insert was only logged as a warning and --create ended with exit status 0;
		the table is made without its sequence column here, so every SNP insert fails'''
	code = ("import CanSNPer2.modules.InitiateDatabase as initiate\n"
		"initiate.snp_annotation_sql = lambda table='snp_annotation': {0!r}\n").format(
		helpers.JOINED_SNP_ANNOTATION_TABLE)
	db = str(tmp_path / "scheme.db")
	result = run_database_cli(tmp_path, create_arguments(tmp_path, db), code=code)
	assert result.returncode == 1, result.stdout + result.stderr
	assert "The canSNP A.1 (rank 1) could not be written" in result.stderr, (result.stdout, result.stderr)
	assert not os.path.exists(db), "an incomplete database must not stay behind"


@pytest.mark.parametrize("prefix", ["\n", "\n\n", "  \n", "\t\n"])
def test_blank_lines_before_the_root_do_not_change_the_database(tmp_path, prefix):
	# a leading blank line gave the root the name "" and a link to a node id 0 that is not
	# in the nodes table, without a message from the source check
	plain = build_db(tmp_path, "plain.db", tree_file("\t"), annotation_file())
	blank = build_db(tmp_path, "blank.db", prefix + tree_file("\t").replace("\n", "\n\n", 1), annotation_file())
	assert tables(blank) == tables(plain)


def modify_arguments(db, mod_file, parent, **extra):
	return ["-db", str(db), "--mod_file", str(mod_file), "--parent", parent]


@pytest.mark.parametrize("rows", [
	"R\tB.1\n",              # B.1 gets a second parent (it is under A.1)
	"B.1\tA.1\n",            # a cycle through A.1 and B.1
	"A.1\tA.1\n",            # a node linked to itself
	"A.1\tR\n",              # a link to the root
])
def test_modification_that_breaks_the_tree_is_refused(tmp_path, rows):
	# these were written without a message and left a database that could not be exported
	db = build_db(tmp_path, "m.db", tree_file("\t"), annotation_file())
	mod = tmp_path / "mod.txt"
	mod.write_text("parent\tchild\n" + rows)
	before = open(db, "rb").read()
	result = run_database_cli(tmp_path, modify_arguments(db, mod, rows.split("\t")[0]))
	assert result.returncode == 1, result.stderr
	assert open(db, "rb").read() == before, "the database is left as it was"
	assert "tree that is not valid" in result.stderr + result.stdout


def test_modification_that_keeps_a_tree_is_accepted(tmp_path):
	db = build_db(tmp_path, "m.db", tree_file("\t"), annotation_file())
	mod = tmp_path / "mod.txt"
	mod.write_text("parent\tchild\nB.1\tC.1\n")
	result = run_database_cli(tmp_path, modify_arguments(db, mod, "B.1"))
	assert result.returncode == 0, result.stderr
