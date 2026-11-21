"""Every bad scheme, database or reference input ends in one message, a non-zero exit and no traceback.

Each case below gives a fixture builder, the command and its arguments, the expected exit
code and a fragment of the expected message. Every case is run as the installed command
would run it, in its own Python process, and checks
- the exit code;
- the message fragment(s) in stdout or stderr;
- that no Python traceback was printed (except where --debug asks for one);
- where the case says so, that a database is byte-identical after the run, or that no
  database or partial (.part) file was left behind.

A few cases replace a function inside the process (inject) to make a write fail, as a full
disk or a locked file would; the replacement code runs before the command starts.

Rows of the table, by group:
1. modification of an existing database (--mod_file, --annotation, --references, --remove)
2. --create and --import when a tree, node or link insert fails; SourceCheck of the
   reference accessions
3. CanSNPer2-download, with --scheme and with --database
4. typing start-up (CanSNPer2)
5. the top-level handler of every command for an unexpected exception"""
import glob
import hashlib
import json
import os
import sqlite3
import subprocess
import sys
from collections import namedtuple

import pytest

import helpers
from CanSNPer2.modules.InitiateDatabase import CanSNPDatabase

TESTS = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(TESTS)

MODULES = {
	"CanSNPer2": "CanSNPer2.CanSNPerTree",
	"CanSNPer2-database": "CanSNPer2.SNPDatabase",
	"CanSNPer2-download": "CanSNPer2.DownloadGenomes",
	"CanSNPer2-test": "CanSNPer2.selftest",
}

## the scheme that the modification and creation cases start from
TREE = "R\nR\tA.1\nR\tA.1\tA.2\nR\tB.1\n"
SNP_HEADER = "\t".join(helpers.SNP_HEADER)
OLD_SNP_HEADER = "snp_id\tstrain\treference\tgenome\tposition\tderived_base\tancestral_base"
SNP_ROWS = ["R\ts1\tPaper2020\tREF1\tCP000001.1\t10\tA\tG",
	"A.1\ts1\tPaper2020\tREF1\tCP000001.1\t100\tT\tC",
	"A.2\ts1\tPaper2021\tREF1\tCP000001.1\t200\tG\tA",
	"B.1\ts2\tPaper2020\tREF2\tCP000002.1\t400\tA\tC"]
REFERENCES_HEADER = "genome\tstrain\tgenbank_id\trefseq_id\tassembly_name\n"
REFERENCES = (REFERENCES_HEADER +
	"REF1\ts1\tGCA_000000001.1\tGCF_000000001.1\tASM1v1\n"
	"REF2\ts2\tGCA_000000002.1\tGCF_000000002.1\tASM2v1\n")
## a SNP on the node C.1 that the modification file adds below R
NEW_SNP = "C.1\ts1\tPaper2022\tREF1\tCP000001.1\t500\tA\tC"
MOD_FILE = "child\tparent\nC.1\tR\n"

## a reference genome of 10,000 bases on one sequence, chr1
GENOME = ">chr1 test\n" + "ACGT" * 2500 + "\n"

MISSING_ACCESSION = "REF: the database has no genbank_id for this genome"

BOOM = "def boom(*args, **kwargs):\n\traise RuntimeError('boom')\n"
TOP_LEVEL = "CanSNPer2 stopped: RuntimeError: boom; run with --debug for details"


def run(command, args, cwd, inject=""):
	'''run command with args in a new Python process, after the code in inject'''
	code = "import sys\nsys.argv = [{command!r}] + sys.argv[1:]\n{inject}\nfrom {module} import main\nmain()\n".format(
		command=command, inject=inject, module=MODULES[command])
	env = dict(os.environ)
	env["PYTHONPATH"] = os.pathsep.join([ROOT, TESTS] + ([env["PYTHONPATH"]] if env.get("PYTHONPATH") else []))
	return subprocess.run([sys.executable, "-c", code] + [str(arg) for arg in args],
		capture_output=True, text=True, cwd=str(cwd), env=env, timeout=300)


def digest(path):
	with open(path, "rb") as handle:
		return hashlib.sha256(handle.read()).hexdigest()


# fixture builders ---------------------------------------------------------------------

def scheme_db(tmp_path, references=REFERENCES):
	'''a database built from the scheme above, and its source files'''
	files = helpers.write_sources(tmp_path / "src", SNP_ROWS, tree=TREE, references=references, header=SNP_HEADER)
	db = str(tmp_path / "scheme.db")
	CanSNPDatabase(database=db, create=True, tree=files["tree"], annotation=files["snp"],
		references=files["references"], source_type="CanSNPer", verbose=False)
	return db, files


def write(path, text):
	path = str(path)
	with open(path, "w") as out:
		out.write(text)
	return path


def modification(snp_rows=(NEW_SNP,), header=SNP_HEADER, references=None, mod=MOD_FILE, parent="R"):
	'''builder: an existing database, a modification file that adds C.1 below R, and an
		annotation file (and a references file) for the new node'''
	def build(tmp_path):
		db, _files = scheme_db(tmp_path)
		args = ["--database", db, "--logs", tmp_path / "logs"]
		if mod is not None:
			args += ["--mod_file", write(tmp_path / "mod.txt", mod)]
		if parent is not None:
			args += ["--parent", parent]
		expect = []
		if snp_rows is not None:
			args += ["--annotation", write(tmp_path / "new_snp.txt", header + "\n" + "".join(row + "\n" for row in snp_rows))]
		if references is not None:
			path = write(tmp_path / "new_references.txt", references)
			args += ["--references", path]
			expect.append(path)
		return {"args": args, "unchanged": [db], "paths": expect}
	return build


def removal(parent):
	def build(tmp_path):
		db, _files = scheme_db(tmp_path)
		return {"args": ["--database", db, "--remove", "--parent", parent, "--logs", tmp_path / "logs"], "unchanged": [db]}
	return build


def creation(option="--create", references=REFERENCES):
	'''builder: the source files of the scheme, and a database path that does not exist yet'''
	def build(tmp_path):
		files = helpers.write_sources(tmp_path / "src", SNP_ROWS, tree=TREE, references=references, header=SNP_HEADER)
		db = tmp_path / "new.db"
		return {"args": [option, "--database", db, "--tree", files["tree"], "--annotation", files["snp"],
			"--references", files["references"], "--logs", tmp_path / "logs"],
			"absent": [str(db)], "paths": [files["references"]]}
	return build


def repo_files(name="test_scheme", status="verified", manifest=True, snp_header=SNP_HEADER, hashed=("snp", "tree", "references")):
	'''the files of a repository laid out as CanSNPer2-data with one scheme'''
	folder = "schemes/{n}/{n}".format(n=name)
	rows = SNP_ROWS if snp_header == SNP_HEADER else [
		"\t".join(row.split("\t")[:4] + row.split("\t")[5:]) for row in SNP_ROWS]
	files = {
		folder + "_snp.txt": snp_header + "\n" + "".join(row + "\n" for row in rows),
		folder + "_tree.txt": TREE,
		folder + "_references.txt": REFERENCES,
		"schemes/{n}/README.md".format(n=name): "# Test\n\n| Field | Value |\n|---|---|\n| Status | {s} |\n".format(s=status),
	}
	if manifest:
		files["database/{n}.manifest.json".format(n=name)] = json.dumps({"sources": {
			folder + "_" + role + ".txt": hashlib.sha256(files[folder + "_" + role + ".txt"].encode()).hexdigest()
			for role in hashed}})
	return files


## replaces the GitHub fetcher with the files of fake_repo.json in the working directory
FAKE_GITHUB = """import json
from test_fetch_scheme import FakeRepo
from CanSNPer2.modules import FetchScheme
FetchScheme.GitHubSchemeSource.__init__.__defaults__ = ("master", FetchScheme.DATA_REPO, FakeRepo(json.load(open("fake_repo.json"))))
"""


def scheme(github=True, extra=(), **repo):
	'''builder: CanSNPer2-download --scheme of a scheme from a fake GitHub repository, or
		from a local directory (github=False)'''
	def build(tmp_path):
		files = repo_files(**repo)
		db = tmp_path / "out" / "test_scheme.db"
		args = ["--scheme", "test_scheme", "--database", db, "--no_references", "--logs", tmp_path / "logs"] + list(extra)
		if github:
			write(tmp_path / "fake_repo.json", json.dumps(files))
		else:
			for path, text in files.items():
				os.makedirs(os.path.dirname(str(tmp_path / "data" / path)), exist_ok=True)
				write(tmp_path / "data" / path, text)
			args += ["--data_dir", tmp_path / "data"]
		return {"args": args, "absent_glob": [str(tmp_path / "out" / "*.part")]}
	return build


def typing_db(path, rows=None, genomes=("REF",)):
	'''a database for typing: one SNP per genome (make_db rows) on chr1, the nodes in a tree
		below the first one, and accessions for every genome'''
	if rows is None:
		rows = [("S.{0}".format(i), "chr1", 501 + 1000 * i, "A", "T") for i in range(len(genomes) + 2)]
	helpers.make_db(path, rows, genome=genomes[0])
	conn = sqlite3.connect(path)
	conn.execute("UPDATE snp_references SET genbank_id='GCA_000000001.1', assembly_name='ASM1v1'")
	for i, genome in enumerate(genomes[1:], start=2):
		conn.execute("INSERT INTO snp_references (id, genome, strain, genbank_id, assembly_name) VALUES (?, ?, ?, ?, ?)",
			(helpers.GENOME_I + i, genome, genome, "GCA_{0:09d}.1".format(i), "ASM{0}v1".format(i)))
		## the last SNPs are placed on the other genomes
		conn.execute("UPDATE snp_annotation SET genome_i=? WHERE node_id=?", (helpers.GENOME_I + i, len(rows) - i + 2))
	conn.execute("INSERT INTO tree (parent, child, rank_i) VALUES (1, 1, 1)")
	for node_id in range(2, len(rows) + 1):
		conn.execute("INSERT INTO tree (parent, child, rank_i) VALUES (1, ?, 1)", (node_id,))
	conn.commit()
	conn.close()
	return path


def typing(sql=(), genomes=("REF", "REF2"), files=("REF", "REF2"), rows=None, extra=()):
	'''builder: CanSNPer2 on one query against a database changed by the sql statements,
		with reference files for the genomes in files'''
	def build(tmp_path):
		db = typing_db(str(tmp_path / "typing.db"), rows=rows, genomes=genomes)
		conn = sqlite3.connect(db)
		for statement in sql:
			conn.execute(statement)
		conn.commit()
		conn.close()
		refdir = tmp_path / "refs"
		refdir.mkdir()
		for name in files:
			write(refdir / (name + ".fna"), GENOME)
		query = write(tmp_path / "sample.fna", GENOME)
		return {"args": ["--database", db, "--refdir", refdir, "--workdir", tmp_path, "--tmpdir", tmp_path / "tmp",
			"--logdir", tmp_path / "logs", "--outdir", tmp_path / "out", query] + list(extra), "unchanged": [db]}
	return build


def download(sql=(), rows=None, extra=()):
	'''builder: CanSNPer2-download --database of a database changed by the sql statements'''
	def build(tmp_path):
		db = typing_db(str(tmp_path / "download.db"), rows=rows, genomes=("REF",))
		conn = sqlite3.connect(db)
		for statement in sql:
			conn.execute(statement)
		conn.commit()
		conn.close()
		return {"args": ["--database", db, "--outdir", tmp_path / "refs", "--logs", tmp_path / "logs"] + list(extra)}
	return build


def exported(extra=()):
	'''builder: CanSNPer2-database --export of a valid database'''
	def build(tmp_path):
		db, _files = scheme_db(tmp_path)
		return {"args": ["--database", db, "--export", "-o", tmp_path / "export", "--logs", tmp_path / "logs"] + list(extra)}
	return build


def nothing(extra=()):
	def build(tmp_path):
		return {"args": list(extra)}
	return build


# injected failures --------------------------------------------------------------------

## a write to snp_annotation fails, as on a full disk
FAIL_SNP_INSERT = """import sqlite3
from CanSNPer2.modules import DatabaseConnection
_insert = DatabaseConnection.TaxonomyDBConnection.insert
def insert(self, data, table):
	if table == "snp_annotation":
		return sqlite3.OperationalError("database or disk is full")
	return _insert(self, data, table)
DatabaseConnection.TaxonomyDBConnection.insert = insert
"""

FAIL_LINK_INSERT = """import sqlite3
from CanSNPer2.modules import DatabaseConnection
DatabaseConnection.DatabaseFunctions.add_link = lambda self, *args, **kwargs: sqlite3.OperationalError("database or disk is full")
"""

FAIL_NODE_INSERT = """import sqlite3
from CanSNPer2.modules import DatabaseConnection
DatabaseConnection.DatabaseFunctions.add_node = lambda self, *args, **kwargs: sqlite3.OperationalError("database or disk is full")
"""

## the finished database cannot be moved into place
FAIL_REPLACE = """import os
_replace = os.replace
def replace(source, target, *args, **kwargs):
	if str(target).endswith(".db"):
		raise PermissionError(13, "Permission denied", str(target))
	return _replace(source, target, *args, **kwargs)
os.replace = replace
"""

FAIL_BUILD_SQLITE = """import sqlite3
from CanSNPer2.modules import InitiateDatabase
def locked(self, *args, **kwargs):
	raise sqlite3.OperationalError("database is locked")
InitiateDatabase.CanSNPDatabase.__init__ = locked
"""


# the table ----------------------------------------------------------------------------

## build(tmp_path) returns {"args": [...], and optionally "unchanged": [files that must stay
## byte-identical], "absent": [files that must not exist afterwards], "absent_glob": [patterns
## that must match nothing], "paths": [values for {0}, {1} ... in the fragments]}
Case = namedtuple("Case", "id command build exit fragments inject forbidden traceback")


def case(id, command, build, exit, fragments, inject="", forbidden=(), traceback=False):
	if isinstance(fragments, str):
		fragments = [fragments]
	return Case(id, command, build, exit, list(fragments), inject, list(forbidden), traceback)


CASES = [
	# 1. modification of an existing database ------------------------------------------
	case("modify: SNP on an unknown node", "CanSNPer2-database",
		modification(snp_rows=[NEW_SNP.replace("C.1", "X.9", 1)]), 1, "snp_id 'X.9' does not appear in the tree"),
	case("modify: SNP on an unknown genome", "CanSNPer2-database",
		modification(snp_rows=[NEW_SNP.replace("REF1", "REF9")]), 1, "genome 'REF9' is not listed"),
	case("modify: short SNP row", "CanSNPer2-database",
		modification(snp_rows=["C.1\ts1\tPaper2022\tREF1\tCP000001.1"]), 1, "expected 8 columns, found 5"),
	case("modify: SNP file with the old header", "CanSNPer2-database",
		modification(header=OLD_SNP_HEADER, snp_rows=["C.1\ts1\tPaper2022\tREF1\t500\tA\tC"]), 1, "--convert"),
	case("modify: SNP on a node that already has one", "CanSNPer2-database",
		modification(snp_rows=[NEW_SNP, "A.1\ts1\tPaper2022\tREF1\tCP000001.1\t900\tA\tC"]), 1,
		"snp_id 'A.1' already has a SNP in"),
	case("modify: SNP insert fails", "CanSNPer2-database",
		modification(), 1, ["The SNP C.1 could not be written", "database or disk is full"], inject=FAIL_SNP_INSERT),
	case("modify: reference genome already in the database", "CanSNPer2-database",
		modification(snp_rows=None, references=REFERENCES_HEADER + "REF2\ts2\tGCA_000000002.1\tGCF_000000002.1\tASM2v1\n"), 1,
		"genome 'REF2' is already in"),
	case("modify: assembly accession of another genome in the database", "CanSNPer2-database",
		modification(snp_rows=None, references=REFERENCES_HEADER + "REF3\ts3\tGCA_000000002.1\t\tASM3v1\n"), 1,
		["{0}:2: genbank_id 'GCA_000000002.1' of genome 'REF3' is already used by genome 'REF2' in", "two genomes cannot share one"]),
	case("modify: short reference row", "CanSNPer2-database",
		modification(snp_rows=None, references=REFERENCES_HEADER + "REF3\ts3\n"), 1, "{0}:2: expected 5 columns, found 2"),
	case("modify: reference row with an empty genbank_id", "CanSNPer2-database",
		modification(snp_rows=None, references=REFERENCES_HEADER + "REF3\ts3\t\tGCF_000000003.1\tASM3v1\n"), 1,
		"{0}:2: genbank_id '' is not an assembly accession"),
	case("modify: unknown --parent", "CanSNPer2-database",
		modification(parent="Z.9"), 1, "The selected parent node (Z.9) could not be found"),
	case("modify: no --parent", "CanSNPer2-database",
		modification(parent=None), 1, "give the branch root with --parent"),
	case("modify: modification file with a malformed row", "CanSNPer2-database",
		modification(mod="child\tparent\nC.1\tR\tfamily\textra\n"), 1, "Malformed row in the modification file"),
	case("modify: modification file without child and parent headers", "CanSNPer2-database",
		modification(mod="a\tb\nC.1\tR\n"), 1, "does not contain proper headers"),
	case("modify: modification file with one column", "CanSNPer2-database",
		modification(mod="child\nC.1\n"), 1, "must contain two or three columns"),
	case("remove: unknown --parent", "CanSNPer2-database",
		removal("Z.9"), 1, "The node Z.9 could not be found"),
	case("remove: the root", "CanSNPer2-database",
		removal("R"), 1, "R is a root of the tree"),

	# 2. --create and --import ----------------------------------------------------------
	case("create: link insert fails", "CanSNPer2-database",
		creation(), 1, ["could not be written", "database or disk is full"], inject=FAIL_LINK_INSERT),
	case("create: node insert fails", "CanSNPer2-database",
		creation(), 1, ["could not be written", "database or disk is full"], inject=FAIL_NODE_INSERT),
	case("import: link insert fails", "CanSNPer2-database",
		creation(option="--import"), 1, ["could not be written", "database or disk is full"], inject=FAIL_LINK_INSERT),
	case("create: reference row with an empty genbank_id", "CanSNPer2-database",
		creation(references=REFERENCES.replace("GCA_000000002.1", "")), 1, "{0}:3: genbank_id '' is not an assembly accession"),
	case("create: reference row with a genbank_id that is not an accession", "CanSNPer2-database",
		creation(references=REFERENCES.replace("GCA_000000002.1", "CP000002.1")), 1,
		"{0}:3: genbank_id 'CP000002.1' is not an assembly accession"),

	# 3. CanSNPer2-download ---------------------------------------------------------------
	case("download --scheme: branch in the old format", "CanSNPer2-download",
		scheme(snp_header=OLD_SNP_HEADER), 1,
		"scheme test_scheme on FOI-Bioinformatics/CanSNPer2-data@master uses positions on the joined reference "
		"(before 2.3.0); use CanSNPer2 2.2 for it, or a branch converted for 2.3.0", inject=FAKE_GITHUB),
	case("download --scheme: local scheme in the old format", "CanSNPer2-download",
		scheme(github=False, snp_header=OLD_SNP_HEADER), 1, "uses positions on the joined reference (before 2.3.0)"),
	case("download --scheme: the database cannot be moved into place", "CanSNPer2-download",
		scheme(github=False), 1, ["could not be written", "Permission denied"], inject=FAIL_REPLACE),
	case("download --scheme: sqlite error while building", "CanSNPer2-download",
		scheme(github=False), 1, ["could not be built", "database is locked"], inject=FAIL_BUILD_SQLITE),
	case("download --scheme: verified scheme without a manifest", "CanSNPer2-download",
		scheme(manifest=False), 1, ["test_scheme.manifest.json was not found", "--allow_unverified"], inject=FAKE_GITHUB),
	case("download --scheme: verified scheme with a file missing from the manifest", "CanSNPer2-download",
		scheme(hashed=("snp", "tree")), 1, ["records no hash for schemes/test_scheme/test_scheme_references.txt",
		"--allow_unverified"], inject=FAKE_GITHUB),
	case("download --scheme: verified scheme without a manifest, --allow_unverified", "CanSNPer2-download",
		scheme(manifest=False, extra=["--allow_unverified"]), 0, "cannot be checked", inject=FAKE_GITHUB),
	case("download --scheme: draft scheme without a manifest", "CanSNPer2-download",
		scheme(manifest=False, status="draft"), 0, "cannot be checked", inject=FAKE_GITHUB),
	case("download --scheme: local verified scheme without a manifest", "CanSNPer2-download",
		scheme(github=False, manifest=False), 0, "cannot be checked"),
	case("download --database: no reference genomes", "CanSNPer2-download",
		download(sql=["DELETE FROM snp_references"]), 1, "nothing to download"),
	case("download --database: NULL accession", "CanSNPer2-download",
		download(sql=["UPDATE snp_references SET genbank_id=NULL"]), 1, MISSING_ACCESSION,
		forbidden=["ended without reporting"]),
	case("download --database: empty accession", "CanSNPer2-download",
		download(sql=["UPDATE snp_references SET genbank_id=''"]), 1, MISSING_ACCESSION,
		forbidden=["ended without reporting"]),
	case("download --database: empty accession and assembly name", "CanSNPer2-download",
		download(sql=["UPDATE snp_references SET genbank_id='', assembly_name=''"]), 1, MISSING_ACCESSION,
		forbidden=["ended without reporting", "CanSNPer2 stopped"]),
	case("download --database: SNPs cannot be read", "CanSNPer2-download",
		download(sql=["UPDATE snp_annotation SET sequence_position='abc' WHERE node_id=2"]), 1,
		"The SNP positions of"),

	# 4. typing start-up ------------------------------------------------------------------
	case("typing: a database genome has no reference file", "CanSNPer2",
		typing(files=("REF",)), 1, ["REF2", "CanSNPer2-download --database"]),
	case("typing: no reference file matches the database", "CanSNPer2",
		typing(files=("OTHER",)), 1, ["REF, REF2", "CanSNPer2-download --database"]),
	case("typing: SNPs cannot be read", "CanSNPer2",
		typing(sql=["UPDATE snp_annotation SET sequence_position='abc' WHERE node_id=2"]), 1,
		"The SNP positions of"),
	case("typing: database without SNPs", "CanSNPer2",
		typing(sql=["DELETE FROM snp_annotation"]), 1, "has no SNPs"),
	case("typing: database without reference genomes", "CanSNPer2",
		typing(sql=["DELETE FROM snp_references"]), 1, "lists no reference genomes"),
	case("typing: SNP on a node missing from nodes", "CanSNPer2",
		typing(sql=["DELETE FROM nodes WHERE id=3"]), 1, ["S.2", "not in the nodes or tree"]),
	case("typing: SNP on a node missing from the tree", "CanSNPer2",
		typing(sql=["DELETE FROM tree WHERE child=3"]), 1, ["S.2", "not in the nodes or tree"]),

	# 5. the top-level handler -----------------------------------------------------------
	case("top level: CanSNPer2-database", "CanSNPer2-database", exported(), 2, TOP_LEVEL,
		inject=BOOM + "import CanSNPer2.SNPDatabase\nCanSNPer2.SNPDatabase.export_source = boom\n"),
	case("top level: CanSNPer2-database --debug", "CanSNPer2-database", exported(extra=["--debug"]), 2,
		["Traceback", "CanSNPer2 stopped: RuntimeError: boom"],
		inject=BOOM + "import CanSNPer2.SNPDatabase\nCanSNPer2.SNPDatabase.export_source = boom\n", traceback=True),
	case("top level: CanSNPer2-database --create removes the started database", "CanSNPer2-database", creation(), 2,
		TOP_LEVEL, inject=BOOM + "from CanSNPer2.modules import InitiateDatabase\n"
		"InitiateDatabase.CanSNPDatabase.load_cansnp_annotation = boom\n"),
	case("top level: CanSNPer2-download", "CanSNPer2-download", download(), 2, TOP_LEVEL,
		inject=BOOM + "import CanSNPer2.DownloadGenomes\nCanSNPer2.DownloadGenomes.DownloadGenomes.run = boom\n"),
	case("top level: CanSNPer2", "CanSNPer2", typing(), 2, TOP_LEVEL,
		inject=BOOM + "import CanSNPer2.modules.CanSNPer2\nCanSNPer2.modules.CanSNPer2.CanSNPer2.run = boom\n"),
	case("top level: CanSNPer2-test", "CanSNPer2-test", nothing(), 2, TOP_LEVEL,
		inject=BOOM + "import shutil\nshutil.which = boom\n"),
]


@pytest.mark.parametrize("case", CASES, ids=[c.id for c in CASES])
def test_bad_input(tmp_path, case):
	built = case.build(tmp_path)
	unchanged = {path: digest(path) for path in built.get("unchanged", [])}
	result = run(case.command, built["args"], tmp_path, inject=case.inject)
	output = result.stdout + result.stderr
	assert result.returncode == case.exit, output
	for fragment in case.fragments:
		fragment = fragment.format(*built.get("paths", []))
		assert fragment in output, "{0!r} not in the output:\n{1}".format(fragment, output)
	for fragment in case.forbidden:
		assert fragment not in output, output
	if not case.traceback:
		assert "Traceback" not in result.stderr, result.stderr
	for path, before in unchanged.items():
		assert digest(path) == before, "{0} was changed by a run that failed".format(path)
	for path in built.get("absent", []):
		assert not os.path.exists(path), "{0} was left behind".format(path)
	for pattern in built.get("absent_glob", []):
		assert glob.glob(pattern) == [], "partial files were left behind"
	if case.exit != 0:
		## no copy of a modified database is left beside it
		assert glob.glob(str(tmp_path / ".*.part")) == [] and glob.glob(str(tmp_path / "*.part")) == []


# the cases that must still work -------------------------------------------------------

def test_a_valid_modification_replaces_the_database(tmp_path):
	built = modification(references=REFERENCES_HEADER + "REF3\ts3\tGCA_000000003.1\tGCF_000000003.1\tASM3v1\n")(tmp_path)
	db = built["unchanged"][0]
	mode = os.stat(db).st_mode
	result = run("CanSNPer2-database", built["args"], tmp_path)
	assert result.returncode == 0, result.stdout + result.stderr
	conn = sqlite3.connect(db)
	assert conn.execute("SELECT sequence, sequence_position FROM snp_annotation WHERE snp_id='C.1'").fetchall() == [("CP000001.1", 500)]
	assert conn.execute("SELECT genome FROM snp_references ORDER BY genome").fetchall() == [("REF1",), ("REF2",), ("REF3",)]
	conn.close()
	assert os.stat(db).st_mode == mode, "the database keeps its mode"
	assert sorted(os.listdir(tmp_path)) == sorted(["src", "scheme.db", "logs", "mod.txt", "new_snp.txt", "new_references.txt"])


def test_a_valid_removal_replaces_the_database(tmp_path):
	built = removal("A.1")(tmp_path)
	db = built["unchanged"][0]
	result = run("CanSNPer2-database", built["args"], tmp_path)
	assert result.returncode == 0, result.stdout + result.stderr
	conn = sqlite3.connect(db)
	assert sorted(n for (n,) in conn.execute("SELECT name FROM nodes")) == ["B.1", "R"]
	conn.close()


def test_nodes_without_snps_are_allowed_at_typing_start_up(tmp_path):
	'''the start-up checks only concern SNPs on missing nodes; a node without a SNP is fine'''
	from CanSNPer2.modules.CanSNPer2 import CanSNPer2
	db = typing_db(str(tmp_path / "typing.db"), genomes=("REF",))
	conn = sqlite3.connect(db)
	conn.execute("INSERT INTO nodes (id, name) VALUES (99, 'EMPTY')")
	conn.execute("INSERT INTO tree (parent, child, rank_i) VALUES (1, 99, 1)")
	conn.commit()
	conn.close()
	(tmp_path / "refs").mkdir()
	write(tmp_path / "refs" / "REF.fna", GENOME)
	obj = CanSNPer2([], workdir=str(tmp_path), refdir="refs", database=db)
	assert obj.check_start() == ["REF.fna"]


def test_a_repeated_tree_path_adds_no_second_node(tmp_path):
	'''a line of the tree file that repeats a known path, or ends at an inner node, used to add
		that node a second time under the same name'''
	files = helpers.write_sources(tmp_path / "src", SNP_ROWS, tree=TREE + "R\tA.1\nR;A.1;A.2\n", references=REFERENCES,
		header=SNP_HEADER)
	db = str(tmp_path / "scheme.db")
	CanSNPDatabase(database=db, create=True, tree=files["tree"], annotation=files["snp"],
		references=files["references"], source_type="CanSNPer", verbose=False)
	conn = sqlite3.connect(db)
	assert sorted(n for (n,) in conn.execute("SELECT name FROM nodes")) == ["A.1", "A.2", "B.1", "R"]
	assert conn.execute("SELECT count(*) FROM tree").fetchone()[0] == 4
	conn.close()


def test_two_genomes_with_the_same_assembly_name_each_get_their_file(tmp_path, monkeypatch):
	'''accessions and assembly names shared one dictionary, so the second genome with an assembly
		name replaced the first: the download ended with exit 0 and only one <genome>.fna'''
	import gzip
	from CanSNPer2.DownloadGenomes import DownloadGenomes
	db = typing_db(str(tmp_path / "download.db"), genomes=("REF", "REF2"))
	conn = sqlite3.connect(db)
	conn.execute("UPDATE snp_references SET assembly_name='ASM1v1'")
	conn.commit()
	conn.close()
	def retrieve(self, link, target, timeout=None):
		with open(target, "wb") as out:
			out.write(gzip.compress(GENOME.encode()))
	monkeypatch.setattr(DownloadGenomes, "retrieve", retrieve)
	refs = tmp_path / "refs"
	assert DownloadGenomes(db, directory=str(refs)).run() == 0
	assert os.path.exists(str(refs / "REF.fna")) and os.path.exists(str(refs / "REF2.fna"))


def test_a_database_given_as_a_symbolic_link_keeps_the_link(tmp_path):
	'''the modified copy replaced the link with a regular file and left the file it pointed to unchanged'''
	built = modification()(tmp_path)
	db = built["unchanged"][0]
	link = tmp_path / "link.db"
	os.symlink(db, str(link))
	args = [str(link) if arg == db else arg for arg in built["args"]]
	result = run("CanSNPer2-database", args, tmp_path)
	assert result.returncode == 0, result.stdout + result.stderr
	assert os.path.islink(str(link)) and os.path.realpath(str(link)) == os.path.realpath(db)
	conn = sqlite3.connect(db)
	assert conn.execute("SELECT count(*) FROM snp_annotation WHERE snp_id='C.1'").fetchone()[0] == 1
	conn.close()


@pytest.mark.skipif(hasattr(os, "geteuid") and os.geteuid() == 0, reason="root can write to any folder")
def test_a_database_in_a_folder_that_cannot_be_written_is_reported(tmp_path):
	'''the copy is made beside the database; a folder without write permission ended with exit 2'''
	built = modification()(tmp_path)
	db = built["unchanged"][0]
	folder = tmp_path / "readonly"
	folder.mkdir()
	moved = str(folder / "scheme.db")
	os.rename(db, moved)
	args = [moved if arg == db else arg for arg in built["args"]]
	before = digest(moved)
	os.chmod(str(folder), 0o555)
	try:
		result = run("CanSNPer2-database", args, tmp_path)
	finally:
		os.chmod(str(folder), 0o755)
	output = result.stdout + result.stderr
	assert result.returncode == 1, output
	assert "Permission denied" in output and "was not changed" in output, output
	assert "Traceback" not in result.stderr
	assert digest(moved) == before and os.listdir(str(folder)) == ["scheme.db"]
