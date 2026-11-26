"""cansnper2, the command structure of 2.6.0: one command, one verb per job.

Four commands with four spellings of the same options were hard to use and easy to
mix up in scripts. cansnper2 offers the same jobs as verbs (type, selftest, db ...,
scheme ..., ref ...) with one spelling for the shared options, calling the same run
functions as the old commands. These tests keep the new structure consistent: every
verb has a --help that stays inside 80 columns without tabs, no flag of a parser is
the abbreviation of another (the audit's prefix trap), and the verbs that replace an
old invocation produce exactly the same output as the old command."""
import os
import subprocess
import sys

import pytest

import helpers
from CanSNPer2 import __version__
from CanSNPer2.modules.InitiateDatabase import CanSNPDatabase

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TREE = "R\nR\tA.1\nR\tA.1\tA.2\nR\tB.1\n"
ROWS = ["R\t1\tP\tREF1\tCP000001.1\t10\tA\tG\t\tyes",
	"A.1\t1\tP\tREF1\tCP000001.1\t20\tT\tC\t\tyes",
	"A.2\t1\tP\tREF1\tCP000001.1\t30\tG\tA\t\tyes",
	"B.1\t1\tP\tREF1\tCP000001.1\t40\tC\tA\t\tyes"]

VERBS = [["type"], ["selftest"],
	["db", "create"], ["db", "import"], ["db", "export"], ["db", "info"],
	["db", "convert"], ["db", "modify"], ["db", "remove"],
	["scheme", "list"], ["scheme", "get"], ["ref", "get"]]


def run_cli(*args, cwd):
	env = dict(os.environ, PYTHONPATH=ROOT)
	return subprocess.run([sys.executable, "-m", "CanSNPer2.cli"] + list(args),
		capture_output=True, text=True, cwd=str(cwd), env=env)


def run_old(module, *args, cwd):
	env = dict(os.environ, PYTHONPATH=ROOT)
	return subprocess.run([sys.executable, "-m", module] + list(args),
		capture_output=True, text=True, cwd=str(cwd), env=env)


@pytest.fixture
def scheme(tmp_path):
	files = helpers.write_sources(tmp_path / "src", ROWS, tree=TREE)
	db = str(tmp_path / "scheme.db")
	CanSNPDatabase(database=db, create=True, tree=files["tree"], annotation=files["snp"],
		references=files["references"], source_type="CanSNPer", verbose=False)
	return db, files


def test_every_verb_has_a_help(tmp_path):
	for verb in VERBS:
		r = run_cli(*(verb + ["--help"]), cwd=tmp_path)
		assert r.returncode == 0, "{verb} --help: {err}".format(verb=" ".join(verb), err=r.stderr)
		assert "usage: cansnper2 " + " ".join(verb) in r.stdout


def test_helps_stay_inside_eighty_columns(tmp_path):
	env = dict(os.environ, PYTHONPATH=ROOT, COLUMNS="80")
	for verb in VERBS + [[]]:
		r = subprocess.run([sys.executable, "-m", "CanSNPer2.cli"] + verb + ["--help"],
			capture_output=True, text=True, cwd=str(tmp_path), env=env)
		long = [line for line in r.stdout.splitlines() if len(line) > 80]
		assert not long, "{verb} --help lines over 80 columns: {0}".format(
			long[:2], verb=" ".join(verb) or "(top level)")


def test_helps_have_no_tabs(tmp_path):
	for verb in VERBS:
		r = run_cli(*(verb + ["--help"]), cwd=tmp_path)
		assert "\t" not in r.stdout, "{verb} --help contains a raw tab".format(
			verb=" ".join(verb))


def test_no_flag_is_the_abbreviation_of_another():
	'''an option that starts with another option of the same parser silently takes
		over its abbreviations: in the old typing parser --annotation was no longer
		reachable as an abbreviation once --annotations existed. An option extended
		with a hyphen (--aligner and --aligner-path) is safe: the shorter one is only
		matched in full, the extension needs its own part'''
	from CanSNPer2.cli import build_parsers
	parser, parsers = build_parsers()
	'''abbreviation happens one parser at a time, so every parser checks its own
		long options'''
	for name, leaf in [("top level", parser)] + sorted(parsers.items()):
		strings = sorted({option for action in leaf._actions
			for option in action.option_strings if option.startswith("--")})
		for position, option in enumerate(strings):
			for other in strings[position + 1:]:
				if other.startswith(option) and not other.startswith(option + "-"):
					raise AssertionError("{0} is an abbreviation of {1} in the {name} parser".format(
						option, other, name=name))


def test_top_help_lists_every_verb(tmp_path):
	r = run_cli("--help", cwd=tmp_path)
	for verb in VERBS:
		assert verb[-1] in r.stdout, "--help does not list {0}".format(" ".join(verb))


def test_cansnper2_version(tmp_path):
	r = run_cli("--version", cwd=tmp_path)
	assert r.returncode == 0
	assert r.stdout.strip() == "CanSNPer2 - version {0}".format(__version__)


def test_without_a_verb_cansnper2_is_a_usage_error(tmp_path):
	'''a bare cansnper2 has to say which verb was missing and leave a status scripts
		can tell apart from a finished run, the old commands print their help instead'''
	r = run_cli(cwd=tmp_path)
	assert r.returncode == 2
	assert "<command>" in (r.stdout + r.stderr)


# 7
def test_db_info_is_the_same_as_the_old_command(tmp_path, scheme):
	db, _ = scheme
	old = run_old("CanSNPer2.SNPDatabase", "--info", "-db", db,
		"--logs", str(tmp_path / "logs_old"), cwd=tmp_path)
	new = run_cli("db", "info", "-d", db, "--log-dir", str(tmp_path / "logs_new"), cwd=tmp_path)
	assert old.returncode == new.returncode == 0
	assert old.stdout == new.stdout


def test_db_export_newick_is_the_same_as_the_old_command(tmp_path, scheme):
	db, _ = scheme
	old = run_old("CanSNPer2.SNPDatabase", "--export", "--export_format", "newick", "-db", db,
		"-o", str(tmp_path / "old_out"), "--logs", str(tmp_path / "logs_old"), cwd=tmp_path)
	new = run_cli("db", "export", "--format", "newick", "-d", db,
		"-o", str(tmp_path / "new_out"), "--log-dir", str(tmp_path / "logs_new"), cwd=tmp_path)
	assert old.returncode == new.returncode == 0
	old_file = (tmp_path / "old_out" / "scheme_tree.nwk").read_text()
	new_file = (tmp_path / "new_out" / "scheme_tree.nwk").read_text()
	assert old_file == new_file


def test_db_create_is_the_same_as_the_old_create(tmp_path, scheme):
	'''the old create flags become the db create verb, and --annotation becomes
		--cansnps because the file it reads lists canSNPs'''
	_, files = scheme
	old = run_old("CanSNPer2.SNPDatabase", "--create", "-db", str(tmp_path / "old.db"),
		"--tree", files["tree"], "--annotation", files["snp"], "--references", files["references"],
		"--scheme_name", "test_scheme", "--logs", str(tmp_path / "logs_old"), cwd=tmp_path)
	new = run_cli("db", "create", "-d", str(tmp_path / "new.db"),
		"--tree", files["tree"], "--cansnps", files["snp"], "--references", files["references"],
		"--scheme-name", "test_scheme", "--log-dir", str(tmp_path / "logs_new"), cwd=tmp_path)
	assert old.returncode == new.returncode == 0
	old_info = run_old("CanSNPer2.SNPDatabase", "--info", "-db", str(tmp_path / "old.db"),
		"--logs", str(tmp_path / "logs_old"), cwd=tmp_path)
	new_info = run_cli("db", "info", "-d", str(tmp_path / "new.db"),
		"--log-dir", str(tmp_path / "logs_new"), cwd=tmp_path)
	assert old_info.stdout == new_info.stdout


def test_selftest_is_the_same_as_the_old_test(tmp_path):
	old = run_old("CanSNPer2.selftest", cwd=tmp_path)
	new = run_cli("selftest", cwd=tmp_path)
	assert old.returncode == new.returncode
	assert old.stdout == new.stdout


def test_import_is_the_same_as_the_old_import(tmp_path, scheme):
	_, files = scheme
	old = run_old("CanSNPer2.SNPDatabase", "--import", "-db", str(tmp_path / "old.db"),
		"--tree", files["tree"], "--annotation", files["snp"], "--references", files["references"],
		"--logs", str(tmp_path / "logs_old"), cwd=tmp_path)
	new = run_cli("db", "import", "-d", str(tmp_path / "new.db"),
		"--tree", files["tree"], "--cansnps", files["snp"], "--references", files["references"],
		"--log-dir", str(tmp_path / "logs_new"), cwd=tmp_path)
	assert old.returncode == new.returncode == 0
	assert "Imported" in new.stdout
	old_info = run_old("CanSNPer2.SNPDatabase", "--info", "-db", str(tmp_path / "old.db"),
		"--logs", str(tmp_path / "logs_old"), cwd=tmp_path)
	new_info = run_cli("db", "info", "-d", str(tmp_path / "new.db"),
		"--log-dir", str(tmp_path / "logs_new"), cwd=tmp_path)
	assert old_info.stdout == new_info.stdout


# read typing and local references through the verbs

import shutil
import sqlite3

READ_OPTIONS = ["--platform", "--samplesheet", "--min-depth", "--min-allele-fraction", "--min-mapq",
	"--min-baseq", "--mixed-sites", "--threads", "--max-reads", "--local-references", "--local-rounds",
	"--screen-reads", "--fetch-local"]
GENOMES = "genome\tnode\trole\nGCA_000000002.1\tA.2\trepresentative\nGCA_000000003.1\tB.1\trepresentative\n"


def test_type_offers_the_read_options(tmp_path):
	r = run_cli("type", "--help", cwd=tmp_path)
	for option in READ_OPTIONS:
		assert option in r.stdout, "cansnper2 type --help does not offer " + option


def test_db_create_and_import_load_a_genomes_file(tmp_path, scheme):
	_, files = scheme
	(tmp_path / "g.txt").write_text(GENOMES)
	for verb in ("create", "import"):
		db = str(tmp_path / (verb + ".db"))
		r = run_cli("db", verb, "-d", db, "--tree", files["tree"], "--cansnps", files["snp"],
			"--references", files["references"], "--genomes", str(tmp_path / "g.txt"),
			"--log-dir", str(tmp_path / "logs"), cwd=tmp_path)
		assert r.returncode == 0, r.stderr
		rows = set(sqlite3.connect(db).execute("SELECT genome, role FROM genomes"))
		assert rows == {("GCA_000000002.1", "representative"), ("GCA_000000003.1", "representative")}


def test_ref_get_and_scheme_get_offer_local_references(tmp_path):
	for verb in (["ref", "get"], ["scheme", "get"]):
		r = run_cli(*(verb + ["--help"]), cwd=tmp_path)
		assert "--local-references" in r.stdout, " ".join(verb)


def test_fetch_local_needs_local_references(tmp_path, scheme):
	db, _ = scheme
	(tmp_path / "r_R1.fq").write_text("@r/1\nACGT\n+\nIIII\n")
	r = run_cli("type", "-d", db, "--fetch-local", str(tmp_path / "r_R1.fq"),
		"--log-dir", str(tmp_path / "logs"), cwd=tmp_path)
	assert r.returncode == 2 and "--local-references" in r.stderr


@pytest.mark.skipif(shutil.which("minimap2") is None or shutil.which("samtools") is None,
	reason="minimap2 and samtools are needed for read typing")
def test_type_reads_is_the_same_as_the_old_command(tmp_path):
	'''reads of a genome typed with cansnper2 type give the call and the site counts of
		the old command'''
	rng = __import__("random").Random(3)
	ref = [rng.choice("ACGT") for _ in range(6000)]
	sites = {"R": 1000, "A.1": 2000, "A.2": 3000, "B.1": 4000}
	sample = list(ref)
	rows = []
	for node, pos in sites.items():
		derived = "ACGT"[("ACGT".index(ref[pos - 1]) + 1) % 4]
		if node != "B.1":
			sample[pos - 1] = derived
		## node rank reference genome sequence sequence_position derived_base ancestral_base score typed
		rows.append("{n}\t1\tP\tREF1\tCP000001.1\t{p}\t{d}\t{a}\t\tyes".format(n=node, p=pos, d=derived, a=ref[pos - 1]))
	files = helpers.write_sources(tmp_path / "src", rows, tree=TREE)
	db = str(tmp_path / "reads.db")
	CanSNPDatabase(database=db, create=True, tree=files["tree"], annotation=files["snp"],
		references=files["references"], source_type="CanSNPer", verbose=False)
	(tmp_path / "refs").mkdir()
	(tmp_path / "refs" / "REF1.fna").write_text(">CP000001.1\n" + "".join(ref) + "\n")
	r1, r2 = helpers.simulate_reads("".join(sample), 30, seed=5)
	(tmp_path / "s_R1.fq").write_text(r1)
	(tmp_path / "s_R2.fq").write_text(r2)
	reads = [str(tmp_path / "s_R1.fq"), str(tmp_path / "s_R2.fq")]
	old = run_old("CanSNPer2.CanSNPerTree", "-db", db, "--refdir", "refs", "-o", "old", "--min_required_hits", "1",
		"--logdir", "logs_old", *reads, cwd=tmp_path)
	new = run_cli("type", "-d", db, "--ref-dir", "refs", "-o", "new", "--min-required-hits", "1",
		"--log-dir", "logs_new", *reads, cwd=tmp_path)
	assert old.returncode == new.returncode == 0, new.stderr
	assert "A.2" in new.stdout and old.stdout.splitlines()[-1] == new.stdout.splitlines()[-1]
	strip = lambda p: [l for l in open(p) if not l.startswith("#")]
	assert strip(tmp_path / "old" / "samples" / "s_site_counts.tsv") == strip(tmp_path / "new" / "samples" / "s_site_counts.tsv")


def test_errors_go_to_standard_error(tmp_path):
	'''a script that reads stdout must not take an error for output: type without a database
		and scheme list of a missing folder printed their errors to stdout'''
	r = run_cli("type", "x.fasta", cwd=tmp_path)
	assert r.returncode == 1 and "A database is required" in r.stderr and r.stdout == ""
	r = run_cli("scheme", "list", "--data-dir", str(tmp_path / "none"), cwd=tmp_path)
	assert r.returncode == 1 and r.stderr.strip() and r.stdout == ""


def test_error_messages_are_not_quoted():
	'''the message of a CanSNPer2 error was printed as its repr, in quotes'''
	from CanSNPer2.modules.CanSNPer2 import CanSNPer2Error
	assert str(CanSNPer2Error("x.fasta is gzip-compressed")) == "x.fasta is gzip-compressed"
