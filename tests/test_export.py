"""CanSNPer2-database --export writes the source files that --create reads.

The export used to write flextaxd's node and link tables: the annotation file held only node
ids and names, without positions or bases, and the tree file had the header "child parent"
over rows that list the parent first. A database could not be rebuilt from it."""
import os
import sqlite3
import subprocess
import sys

import pytest

import helpers
from CanSNPer2.modules.ExportDatabase import ExportError, export_source
from CanSNPer2.modules.InitiateDatabase import CanSNPDatabase
from CanSNPer2.modules.DatabaseConnection import ANNOTATION_SOURCE_COLUMNS
from CanSNPer2.modules.SourceCheck import check_sources

TREE = "R\nR\tA.1\nR\tA.1\tA.2\nR\tA.1\tA.3\nR\tB.1\nR\tB.1\tB.2\n"
SNPS = ("\t".join(helpers.SNP_HEADER) + "\tdate\n"
	"R\tsp\tPaper2020\tREF1\tCP000001.1\t10\tA\tG\t01/02/20\n"
	"A.1\tsp\tPaper2020\tREF1\tCP000001.1\t100\tT\tC\t01/02/20\n"
	"A.2\tsp\tPaper2021\tREF2\tCP000002.1\t200\tG\tA\t03/04/21\n"
	"A.3\tsp\tPaper2021\tREF1\tCP000001.1\t300\tC\tT\t03/04/21\n"
	"B.1\tsp\tPaper2020\tREF2\tCP000002.1\t400\tA\tC\t01/02/20\n")
# REF2 has no assembly name: its row ends with an empty field
REFERENCES = ("genome\tstrain\tgenbank_id\trefseq_id\tassembly_name\n"
	"REF1\tstrain one\tGCA_1.1\tGCF_1.1\tASM1v1\n"
	"REF2\tstrain two\tGCA_2.1\tGCF_2.1\t\n")


def create(db, tree, snp, references):
	CanSNPDatabase(database=db, create=True, tree=tree, annotation=snp, references=references,
		source_type="CanSNPer", verbose=False)


def source_db(tmp_path):
	header, _, body = SNPS.partition("\n")
	files = helpers.write_sources(tmp_path, body.splitlines(), tree=TREE, references=REFERENCES,
		name="src", header=header)
	db = str(tmp_path / "scheme.db")
	create(db, files["tree"], files["snp"], files["references"])
	return db


def content(db):
	conn = sqlite3.connect(db)
	result = {
		"nodes": conn.execute("SELECT id,name FROM nodes ORDER BY id").fetchall(),
		"links": conn.execute("""SELECT c.name, p.name, r.rank FROM tree t JOIN nodes c ON c.id=t.child
			JOIN nodes p ON p.id=t.parent JOIN rank r ON r.rank_i=t.rank_i ORDER BY c.name""").fetchall(),
		"snps": conn.execute("""SELECT a.node_id, a.snp_id, a.sequence, a.sequence_position, a.ancestral_base, a.derived_base,
			a.reference, g.genome FROM snp_annotation a JOIN snp_references g ON g.id=a.genome_i
			ORDER BY a.node_id""").fetchall(),
		"references": conn.execute("""SELECT genome,strain,genbank_id,refseq_id,assembly_name
			FROM snp_references ORDER BY genome""").fetchall(),
	}
	conn.close()
	return result


def test_an_export_rebuilds_the_same_database(tmp_path):
	db = source_db(tmp_path)
	files = export_source(db, str(tmp_path / "export"))
	assert sorted(os.path.basename(f) for f in files.values()) == \
		["scheme_references.txt", "scheme_snp.txt", "scheme_tree.txt"]
	rebuilt = str(tmp_path / "rebuilt.db")
	create(rebuilt, files["tree"], files["snp"], files["references"])
	original = content(db)
	assert len(original["snps"]) == 5 and len(original["references"]) == 2
	assert ("REF2", "strain two", "GCA_2.1", "GCF_2.1", "") in original["references"]
	assert content(rebuilt) == original
	# exporting the rebuilt database gives the same files
	again = export_source(rebuilt, str(tmp_path / "again"), name="scheme")
	for key in files:
		assert open(again[key]).read() == open(files[key]).read()


def test_the_files_have_the_source_format(tmp_path):
	files = export_source(source_db(tmp_path), str(tmp_path / "export"))
	assert open(files["tree"]).read() == TREE
	snp = open(files["snp"]).read().splitlines()
	assert snp[0] == "\t".join(helpers.SNP_HEADER)
	# the strain is that of the SNP's reference genome; the date is not exported
	assert snp[2] == "A.1\tstrain one\tPaper2020\tREF1\tCP000001.1\t100\tT\tC"
	assert snp[3] == "A.2\tstrain two\tPaper2021\tREF2\tCP000002.1\t200\tG\tA"
	assert open(files["references"]).read() == REFERENCES


# the header that CanSNPer2-data/validate.py requires of a scheme's SNP file
CANSNPER2_DATA_SNP_HEADER = ["snp_id", "strain", "reference", "genome", "sequence",
	"sequence_position", "derived_base", "ancestral_base"]


def test_the_snp_header_is_the_one_CanSNPer2_data_requires(tmp_path):
	files = export_source(source_db(tmp_path), str(tmp_path / "export"))
	assert open(files["snp"]).readline().rstrip("\n").split("\t") == CANSNPER2_DATA_SNP_HEADER
	# and it still carries every column the loaders need
	assert set(ANNOTATION_SOURCE_COLUMNS) <= set(CANSNPER2_DATA_SNP_HEADER)


def test_an_export_passes_the_source_checks(tmp_path):
	files = export_source(source_db(tmp_path), str(tmp_path / "export"))
	assert check_sources(tree=files["tree"], annotation=files["snp"], references=files["references"]) == []


def test_an_export_passes_CanSNPer2_data_validate(tmp_path):
	"""run validate.py's SNP, tree and reference checks on an export, when a CanSNPer2-data
		checkout sits next to CanSNPer2 (or CANSNPER2_DATA points to one)"""
	data = os.environ.get("CANSNPER2_DATA", os.path.join(os.path.dirname(__file__), "..", "..", "CanSNPer2-data"))
	script = os.path.join(data, "validate.py")
	if not os.path.exists(script):
		pytest.skip("no CanSNPer2-data checkout")
	import importlib.util
	spec = importlib.util.spec_from_file_location("cansnper2_data_validate", script)
	validate = importlib.util.module_from_spec(spec)
	spec.loader.exec_module(validate)
	if "sequence_position" not in getattr(validate, "SNP_HEADER", ()):
		pytest.skip("the CanSNPer2-data checkout still validates joined positions")
	# a scheme as published there: accessions of the checked form, an assembly name on every
	# reference (an empty last column is trailing whitespace to validate.py) and a SNP on
	# every node (B.2 has none in the shared fixture)
	db = str(tmp_path / "scheme.db")
	files = {}
	references = (REFERENCES.replace("GCA_1.1", "GCA_000000001.1").replace("GCF_1.1", "GCF_000000001.1")
		.replace("GCA_2.1", "GCA_000000002.1").replace("GCF_2.1", "GCF_000000002.1").replace("\t\n", "\tASM2v1\n"))
	for key, text in (("tree", TREE.replace("R\tB.1\tB.2\n", "")), ("snp", SNPS), ("references", references)):
		files[key] = str(tmp_path / ("src_" + key + ".txt"))
		with open(files[key], "w") as out:
			out.write(text)
	create(db, files["tree"], files["snp"], files["references"])
	export_source(db, str(tmp_path / "schemes" / "scheme"))
	scheme = validate.Scheme("scheme")
	scheme.path = lambda relative: str(tmp_path / relative)
	problems, warnings = [], []
	nodes = validate.check_tree(problems, scheme)
	genomes = validate.check_references(problems, scheme)
	validate.check_snp(problems, warnings, scheme, nodes, genomes)
	assert problems == []


def test_a_node_with_two_parents_is_refused(tmp_path):
	db = source_db(tmp_path)
	conn = sqlite3.connect(db)
	a2, b1 = (conn.execute("SELECT id FROM nodes WHERE name=?", (n,)).fetchone()[0] for n in ("A.2", "B.1"))
	conn.execute("INSERT INTO tree (parent, child, rank_i) VALUES (?, ?, 1)", (b1, a2))
	conn.commit()
	conn.close()
	with pytest.raises(ExportError, match="A.2 has more than one parent"):
		export_source(db, str(tmp_path / "export"))
	assert not os.path.exists(tmp_path / "export"), "nothing is written for a tree that cannot be exported"


def test_a_snp_without_its_reference_genome_is_refused(tmp_path):
	db = source_db(tmp_path)
	conn = sqlite3.connect(db)
	conn.execute("UPDATE snp_annotation SET genome_i=99 WHERE snp_id='A.3'")
	conn.commit()
	conn.close()
	with pytest.raises(ExportError, match="1 SNP\\(s\\) point to a reference genome .*A.3"):
		export_source(db, str(tmp_path / "export"))


@pytest.mark.parametrize("export_format", ["source", "tab"])
def test_the_command_line_export(tmp_path, export_format):
	db = source_db(tmp_path)
	out = tmp_path / ("cli_" + export_format)
	result = subprocess.run([sys.executable, "-m", "CanSNPer2.SNPDatabase", "--database", db, "--export",
		"--export_format", export_format, "-o", str(out), "--logs", str(tmp_path / "logs")],
		capture_output=True, text=True, cwd=str(tmp_path))
	assert result.returncode == 0, result.stderr
	assert sorted(os.listdir(out)) == ["scheme_references.txt", "scheme_snp.txt", "scheme_tree.txt"]
	assert "Import back with: CanSNPer2-database --import" in result.stdout
