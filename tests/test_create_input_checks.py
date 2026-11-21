"""CanSNPer2-database --create refuses source files that do not form a valid scheme.

The loaders used to accept almost anything: a position "abc" and bases "a"/"a" were
stored, a SNP on an unknown genome or a duplicated snp_id was dropped with a warning,
a node given two parents became two nodes, and a tree in the wrong layout built a
different tree than intended - all with exit status 0. The rules below follow
CanSNPer2-data/validate.py, so a scheme that passes there also builds here.

Every refusal must name the problem and leave no database file behind; all problems
of a file set are reported together, so they can be fixed in one pass."""
import os
import subprocess
import sys

import sqlite3

import pytest

import helpers
from CanSNPer2.modules.DatabaseConnection import SourceFileError
from CanSNPer2.modules.ImportDatabase import DatabaseImportError, import_from_sources
from CanSNPer2.modules.InitiateDatabase import CanSNPDatabase
from CanSNPer2.modules.SourceCheck import check_sources

TREE = "R\nR\tA.1\nR\tA.1\tA.2\nR\tB.1\n"
SNP_HEADER = "\t".join(helpers.SNP_HEADER) + "\n"
SNP_ROWS = ("R\ts1\tPaper2020\tREF1\tCP000001.1\t10\tA\tG\n"
	"A.1\ts1\tPaper2020\tREF1\tCP000001.1\t100\tT\tC\n"
	"A.2\ts1\tPaper2021\tREF1\tCP000001.1\t200\tG\tA\n"
	"B.1\ts2\tPaper2020\tREF2\tCP000002.1\t400\tA\tC\n")
REFERENCES = ("genome\tstrain\tgenbank_id\trefseq_id\tassembly_name\n"
	"REF1\ts1\tGCA_000000001.1\tGCF_000000001.1\tASM1v1\n"
	"REF2\ts2\tGCA_000000002.1\tGCF_000000002.1\tASM2v1\n")


def write_sources(tmp_path, tree=TREE, snp=SNP_HEADER + SNP_ROWS, references=REFERENCES):
	'''the whole SNP file text is given here, header included, so a test can change any line'''
	header, _, body = snp.partition("\n")
	rows = body.split("\n")
	if rows and rows[-1] == "":
		rows.pop()
	return helpers.write_sources(tmp_path, rows, tree=tree, references=references, header=header)


def create(tmp_path, files, name="scheme.db"):
	db = str(tmp_path / name)
	CanSNPDatabase(database=db, create=True, tree=files["tree"], annotation=files["snp"],
		references=files["references"], source_type="CanSNPer", verbose=False)
	return db


def refused(tmp_path, files):
	'''build the files, expect a refusal, and return its message'''
	db = str(tmp_path / "scheme.db")
	with pytest.raises(SourceFileError) as error:
		create(tmp_path, files)
	assert not os.path.exists(db), "a refused build must not leave a database behind"
	return str(error.value)


def with_snp_row(row):
	return SNP_HEADER + SNP_ROWS + row + "\n"


def test_a_valid_scheme_still_builds(tmp_path):
	db = create(tmp_path, write_sources(tmp_path))
	assert os.path.exists(db)


# SNP file ---------------------------------------------------------------------------

@pytest.mark.parametrize("position", ["abc", "0", "-5", "1.5", "12a", ""])
def test_a_position_must_be_a_positive_integer(tmp_path, position):
	snp = SNP_HEADER + SNP_ROWS.replace("\t100\t", "\t{0}\t".format(position))
	message = refused(tmp_path, write_sources(tmp_path, snp=snp))
	assert "position '{0}' is not a positive integer".format(position) in message


@pytest.mark.parametrize("column,base", [
	("derived_base", "a"), ("derived_base", "AT"), ("derived_base", "N"), ("derived_base", "-"),
	("ancestral_base", "c"), ("ancestral_base", "N"), ("ancestral_base", ""),
])
def test_a_base_must_be_one_of_ACGT(tmp_path, column, base):
	if column == "derived_base":
		row = "A.1\ts1\tPaper2020\tREF1\tCP000001.1\t100\t{0}\tC".format(base)
	else:
		row = "A.1\ts1\tPaper2020\tREF1\tCP000001.1\t100\tT\t{0}".format(base)
	snp = SNP_HEADER + SNP_ROWS.replace("A.1\ts1\tPaper2020\tREF1\tCP000001.1\t100\tT\tC", row)
	message = refused(tmp_path, write_sources(tmp_path, snp=snp))
	assert "{0} '{1}' is not one of A, C, G, T".format(column, base) in message


def test_derived_and_ancestral_base_must_differ(tmp_path):
	snp = SNP_HEADER + SNP_ROWS.replace("\t100\tT\tC", "\t100\tG\tG")
	message = refused(tmp_path, write_sources(tmp_path, snp=snp))
	assert "derived and ancestral base are both 'G'" in message


def test_a_snp_on_a_genome_missing_from_the_references_is_refused(tmp_path):
	snp = SNP_HEADER + SNP_ROWS.replace("B.1\ts2\tPaper2020\tREF2", "B.1\ts2\tPaper2020\tREF9")
	message = refused(tmp_path, write_sources(tmp_path, snp=snp))
	assert "genome 'REF9' is not listed in" in message


def test_a_snp_id_listed_twice_is_refused(tmp_path):
	message = refused(tmp_path, write_sources(tmp_path,
		snp=with_snp_row("A.1\ts1\tPaper2020\tREF1\tCP000001.1\t150\tT\tC")))
	assert "snp_id 'A.1' is listed 2 times" in message


def test_a_position_annotated_twice_on_one_genome_is_refused(tmp_path):
	tree = TREE + "R\tB.1\tB.2\n"
	message = refused(tmp_path, write_sources(tmp_path, tree=tree,
		snp=with_snp_row("B.2\ts1\tPaper2020\tREF1\tCP000001.1\t100\tA\tC")))
	assert "position 100 on sequence 'CP000001.1' of 'REF1' is annotated 2 times" in message


def test_the_same_position_on_two_genomes_is_allowed(tmp_path):
	'''positions are per reference genome, 100 on REF1 and 100 on REF2 are different sites'''
	tree = TREE + "R\tB.1\tB.2\n"
	create(tmp_path, write_sources(tmp_path, tree=tree,
		snp=with_snp_row("B.2\ts2\tPaper2020\tREF2\tCP000002.1\t100\tA\tC")))


# sequence and position on that sequence ------------------------------------------------

def checked(files):
	return check_sources(tree=files["tree"], annotation=files["snp"], references=files["references"])


def test_create_stores_sequence_and_sequence_position(tmp_path):
	rows = [("R", "s1", "Paper2020", "REF1", "chrI", 10, "A", "G"),
		("A.1", "s1", "Paper2020", "REF1", "chrII", 10, "T", "C"),
		("B.1", "s1", "Paper2020", "REF1", "chrII", 250, "G", "A")]
	files = helpers.write_sources(tmp_path, rows)
	db = create(tmp_path, files)
	conn = sqlite3.connect(db)
	stored = conn.execute("SELECT snp_id, sequence, sequence_position FROM snp_annotation ORDER BY snp_id").fetchall()
	columns = [row[1] for row in conn.execute("PRAGMA table_info(snp_annotation)")]
	conn.close()
	assert stored == [("A.1", "chrII", 10), ("B.1", "chrII", 250), ("R", "chrI", 10)]
	assert "sequence" in columns and "sequence_position" in columns
	assert "position" not in columns


def test_old_header_is_refused_with_the_conversion_hint(tmp_path):
	header = ("snp_id", "strain", "reference", "genome", "position", "derived_base", "ancestral_base")
	files = helpers.write_sources(tmp_path, [("R", "s1", "Paper2020", "REF1", 10, "A", "G")], header=header)
	problems = checked(files)
	assert len(problems) == 1, problems
	assert "--convert" in problems[0]


def test_same_position_on_two_sequences_is_allowed(tmp_path):
	rows = [("R", "s1", "Paper2020", "REF1", "chrI", 50, "A", "G"),
		("A.1", "s1", "Paper2020", "REF1", "chrII", 50, "T", "C")]
	assert checked(helpers.write_sources(tmp_path, rows)) == []


def test_same_sequence_and_position_twice_is_refused(tmp_path):
	rows = [("R", "s1", "Paper2020", "REF1", "chrII", 50, "A", "G"),
		("A.1", "s1", "Paper2020", "REF1", "chrII", 50, "T", "C")]
	problems = checked(helpers.write_sources(tmp_path, rows))
	assert any("position 50 on sequence 'chrII' of 'REF1' is annotated 2 times" in p for p in problems), problems


def test_empty_sequence_is_refused(tmp_path):
	rows = [("R", "s1", "Paper2020", "REF1", "chrI", 10, "A", "G"),
		("A.1", "s1", "Paper2020", "REF1", "", 20, "T", "C")]
	files = helpers.write_sources(tmp_path, rows)
	problems = checked(files)
	assert problems == ["{0}:3: sequence is empty, it must name a sequence of the reference genome".format(files["snp"])], problems


def test_a_snp_whose_node_is_not_in_the_tree_is_refused(tmp_path):
	message = refused(tmp_path, write_sources(tmp_path,
		snp=with_snp_row("X.9\ts1\tPaper2020\tREF1\tCP000001.1\t900\tA\tC")))
	assert "snp_id 'X.9' does not appear in" in message


def test_a_short_snp_row_is_refused(tmp_path):
	message = refused(tmp_path, write_sources(tmp_path, snp=with_snp_row("B.1\ts2\tPaper2020")))
	assert "expected 8 columns, found 3" in message


def test_comments_blank_lines_and_a_comment_column_are_still_accepted(tmp_path):
	'''the relaxations of issue #22 stay: # lines, blank lines, and a trailing comment
		column that rows may leave out'''
	snp = (SNP_HEADER.rstrip("\n") + "\tcomment\n"
		+ "# a comment line\n\n"
		+ SNP_ROWS.replace("\t100\tT\tC\n", "\t100\tT\tC\ta named branch\n"))
	create(tmp_path, write_sources(tmp_path, snp=snp))


def test_a_node_without_a_snp_is_allowed(tmp_path):
	'''the tree root and grouping nodes may carry no SNP; that is reported by
		CanSNPer2-data/validate.py for published schemes, not refused by the loader'''
	create(tmp_path, write_sources(tmp_path, tree=TREE + "R\tB.1\tB.2\n"))


# tree file --------------------------------------------------------------------------

def test_a_node_with_two_parents_is_refused(tmp_path):
	tree = TREE + "R\tB.1\tA.2\n"
	message = refused(tmp_path, write_sources(tmp_path, tree=tree))
	assert "node 'A.2' is given more than one parent: A.1, B.1" in message


def test_every_line_must_start_at_the_same_root(tmp_path):
	tree = TREE + "Q\tA.1\n"
	message = refused(tmp_path, write_sources(tmp_path, tree=tree))
	assert "every line must start at the root 'R'" in message and "'Q'" in message


def test_a_child_parent_table_given_as_tree_is_refused(tmp_path):
	'''the --mod_file layout (child, parent per line, with a header) is not a tree file'''
	tree = "child\tparent\nA.1\tR\nA.2\tA.1\nB.1\tR\n"
	message = refused(tmp_path, write_sources(tmp_path, tree=tree))
	assert "every line must start at the root" in message


def test_a_space_separated_tree_is_refused(tmp_path):
	'''spaces are allowed inside node names (Francisellaceae uses "Francisella
		tularensis"), so a space-separated tree reads as one long node name; the SNPs then
		name nodes that are not in the tree'''
	tree = "R\nR A.1\nR A.1 A.2\nR B.1\n"
	message = refused(tmp_path, write_sources(tmp_path, tree=tree))
	assert "does not appear in" in message


def test_an_empty_column_between_nodes_is_refused(tmp_path):
	tree = TREE + "R\t\tB.2\n"
	message = refused(tmp_path, write_sources(tmp_path, tree=tree))
	assert "empty column 2 between populated columns" in message


def test_semicolon_trees_and_repeated_paths_are_still_accepted(tmp_path):
	'''CanSNPer1 trees use ";", and repeating a path that is already known is harmless'''
	tree = TREE.replace("\t", ";") + "R;A.1\n"
	create(tmp_path, write_sources(tmp_path, tree=tree))


# references file --------------------------------------------------------------------

def test_a_genome_listed_twice_in_the_references_is_refused(tmp_path):
	references = REFERENCES + "REF1\ts3\tGCA_000000003.1\tGCF_000000003.1\tASM3v1\n"
	message = refused(tmp_path, write_sources(tmp_path, references=references))
	assert "genome 'REF1' is listed 2 times" in message


def test_a_genbank_id_used_by_two_genomes_is_refused(tmp_path):
	references = REFERENCES + "REF3\ts3\tGCA_000000001.1\tGCF_000000003.1\tASM3v1\n"
	message = refused(tmp_path, write_sources(tmp_path, references=references))
	assert "references.txt:4: genbank_id 'GCA_000000001.1' of genome 'REF3' is already used by genome 'REF1' on line 2" in message


def test_the_old_snp_header_message_is_the_shared_conversion_hint(tmp_path):
	from CanSNPer2.modules.DatabaseConnection import JOINED_FORMAT_HINT
	from CanSNPer2.modules.SourceCheck import check_snps
	old = tmp_path / "old_snp.txt"
	old.write_text("snp_id\tstrain\treference\tgenome\tposition\tderived_base\tancestral_base\n")
	problems = []
	check_snps(str(old), problems)
	assert len(problems) == 1 and JOINED_FORMAT_HINT in problems[0]


def test_a_short_reference_row_is_refused(tmp_path):
	references = REFERENCES + "REF3\ts3\n"
	message = refused(tmp_path, write_sources(tmp_path, references=references))
	assert "expected 5 columns, found 2" in message


def test_an_empty_assembly_name_is_allowed(tmp_path):
	'''exported databases may carry an empty last column (issue of #63)'''
	references = REFERENCES.replace("\tASM2v1\n", "\t\n")
	create(tmp_path, write_sources(tmp_path, references=references))


def test_snps_without_a_references_file_are_refused(tmp_path):
	files = write_sources(tmp_path)
	db = str(tmp_path / "scheme.db")
	with pytest.raises(SourceFileError, match="needs the references file"):
		CanSNPDatabase(database=db, create=True, tree=files["tree"], annotation=files["snp"],
			references=False, source_type="CanSNPer", verbose=False)
	assert not os.path.exists(db)


# reporting --------------------------------------------------------------------------

def test_all_problems_are_reported_together_with_file_and_line(tmp_path):
	snp = SNP_HEADER + SNP_ROWS.replace("\t100\t", "\tabc\t").replace("\t400\tA\tC", "\t400\tA\tA")
	files = write_sources(tmp_path, snp=snp, tree=TREE + "R\tB.1\tA.2\n")
	message = refused(tmp_path, files)
	assert "{0}:3: sequence_position 'abc'".format(files["snp"]) in message
	assert "{0}:5: derived and ancestral base are both 'A'".format(files["snp"]) in message
	assert "node 'A.2' is given more than one parent" in message


def test_import_refuses_the_same_files(tmp_path):
	snp = SNP_HEADER + SNP_ROWS.replace("\t100\t", "\tabc\t")
	files = write_sources(tmp_path, snp=snp)
	target = str(tmp_path / "imported.db")
	with pytest.raises(DatabaseImportError, match="position 'abc'"):
		import_from_sources(files["tree"], files["snp"], files["references"], target)
	assert not os.path.exists(target)


def test_the_command_line_create_exits_1_and_names_the_problem(tmp_path):
	snp = SNP_HEADER + SNP_ROWS.replace("\t100\tT\tC", "\tabc\ta\ta")
	files = write_sources(tmp_path, snp=snp)
	db = tmp_path / "cli.db"
	result = subprocess.run([sys.executable, "-m", "CanSNPer2.SNPDatabase", "--create",
		"--database", str(db), "--tree", files["tree"], "--annotation", files["snp"],
		"--references", files["references"], "--logs", str(tmp_path / "logs")],
		capture_output=True, text=True, cwd=str(tmp_path))
	assert result.returncode == 1, result.stdout + result.stderr
	assert "position 'abc' is not a positive integer" in result.stderr
	assert "derived_base 'a' is not one of A, C, G, T" in result.stderr
	assert not db.exists()


# published schemes ------------------------------------------------------------------

DATA = os.environ.get("CANSNPER2_DATA",
	os.path.join(os.path.dirname(__file__), "..", "..", "CanSNPer2-data"))


def published_schemes():
	folder = os.path.join(DATA, "schemes")
	if not os.path.isdir(folder):
		return []
	names = []
	for name in sorted(os.listdir(folder)):
		snp_file = os.path.join(folder, name, name + "_snp.txt")
		if not os.path.exists(snp_file):
			continue
		with open(snp_file) as snps:
			header = snps.readline().rstrip("\r\n").lstrip("#").split("\t")
		if "sequence_position" not in header:
			'''a checkout with joined positions (the format before 2.3.0) is refused by the input
				checks; it builds once the scheme is converted with CanSNPer2-database --convert'''
			names.append(pytest.param(name, marks=pytest.mark.skip(
				reason="published with joined positions, not yet converted")))
		else:
			names.append(name)
	return names


@pytest.mark.skipif(not published_schemes(), reason="no CanSNPer2-data checkout next to CanSNPer2 (set CANSNPER2_DATA)")
@pytest.mark.parametrize("name", published_schemes())
def test_the_schemes_of_CanSNPer2_data_build(tmp_path, name):
	folder = os.path.join(DATA, "schemes", name)
	files = {key: os.path.join(folder, "{0}_{1}.txt".format(name, key)) for key in ("tree", "snp", "references")}
	create(tmp_path, files)


def test_the_command_line_mod_file_refuses_an_old_header_with_the_conversion_hint(tmp_path):
	'''--mod_file with --annotation reads the SNP file through check_source_headers, not
		check_sources; a file with joined positions must still name the conversion command'''
	db = create(tmp_path, write_sources(tmp_path))
	mod = tmp_path / "mod.txt"
	mod.write_text("child\tparent\nC.1\tR\n")
	old = tmp_path / "old_snp.txt"
	old.write_text("snp_id\tstrain\treference\tgenome\tposition\tderived_base\tancestral_base\n"
		"C.1\ts1\tPaper2020\tREF1\t500\tA\tC\n")
	result = subprocess.run([sys.executable, "-m", "CanSNPer2.SNPDatabase", "--database", db,
		"--mod_file", str(mod), "--parent", "R", "--annotation", str(old), "--logs", str(tmp_path / "logs")],
		capture_output=True, text=True, cwd=str(tmp_path))
	assert result.returncode == 1, result.stdout + result.stderr
	assert "--convert" in result.stderr, result.stderr


def test_check_source_headers_refuses_joined_positions_with_the_conversion_hint():
	from CanSNPer2.modules.DatabaseConnection import ANNOTATION_SOURCE_COLUMNS, JOINED_FORMAT_HINT, check_source_headers
	headers = ["snp_id", "strain", "reference", "genome", "position", "derived_base", "ancestral_base"]
	with pytest.raises(SourceFileError) as error:
		check_source_headers(headers, ANNOTATION_SOURCE_COLUMNS, "old_snp.txt")
	assert JOINED_FORMAT_HINT in str(error.value) and "old_snp.txt" in str(error.value)


## a genome name that is also the accession of another genome, and whitespace in the source files

def test_genome_name_equal_to_another_genomes_accession_is_refused(tmp_path):
	# the loaders kept genome, genbank_id and refseq_id of every genome in one dict, so the
	# name REF2 below resolved to the genome that has it as an accession and the SNPs of
	# genome REF2 were attached to the wrong genome
	references = ("genome\tstrain\tgenbank_id\trefseq_id\tassembly_name\n"
		"REF1\ts1\tGCA_000000001.1\tREF2\tASM1v1\n"
		"REF2\ts2\tGCA_000000002.1\tGCF_000000002.1\tASM2v1\n")
	files = write_sources(tmp_path, references=references)
	message = refused(tmp_path, files)
	assert "REF2" in message and "REF1" in message
	problems = check_sources(tree=files["tree"], annotation=files["snp"], references=files["references"])
	assert any("'REF2'" in p and "'REF1'" in p for p in problems), problems


def test_genome_name_equal_to_its_own_accession_is_accepted(tmp_path):
	references = ("genome\tstrain\tgenbank_id\trefseq_id\tassembly_name\n"
		"REF1\ts1\tGCA_000000001.1\tREF1\tASM1v1\n"
		"REF2\ts2\tGCA_000000002.1\tGCF_000000002.1\tASM2v1\n")
	files = write_sources(tmp_path, references=references)
	assert check_sources(tree=files["tree"], annotation=files["snp"], references=files["references"]) == []


def test_loader_refuses_a_name_that_maps_to_two_genomes(tmp_path):
	# the loader itself, without the source check in front of it
	from CanSNPer2.modules.InitiateDatabase import CreateSNPDatabase
	references = tmp_path / "references.txt"
	references.write_text("genome\tstrain\tgenbank_id\trefseq_id\tassembly_name\n"
		"REF1\ts1\tGCA_000000001.1\tREF2\tASM1v1\n"
		"REF2\ts2\tGCA_000000002.1\tGCF_000000002.1\tASM2v1\n")
	db = str(tmp_path / "x.db")
	CreateSNPDatabase(db, create=True)
	loader = CanSNPDatabase(database=db, create=False, tree=None, annotation=None, references=None,
		source_type="CanSNPer", verbose=False)
	with pytest.raises(SourceFileError) as error:
		loader.load_genome_annotation_file(str(references))
	assert "REF1" in str(error.value) and "REF2" in str(error.value)


def test_modification_refuses_a_name_that_maps_to_two_genomes(tmp_path):
	from CanSNPer2.modules.SourceCheck import check_modification
	files = write_sources(tmp_path)
	db = create(tmp_path, files)
	references = tmp_path / "new_references.txt"
	references.write_text("genome\tstrain\tgenbank_id\trefseq_id\tassembly_name\n"
		"REF3\ts3\tGCA_000000003.1\tREF1\tASM3v1\n")
	problems = check_modification(db, references=str(references))
	assert any("'REF1'" in p and "'REF3'" in p for p in problems), problems
