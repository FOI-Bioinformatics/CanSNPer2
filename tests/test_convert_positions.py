'''CanSNPer2-database --convert: SNP positions along the joined reference sequences (the
format before 2.3.0) are converted to a sequence and a position on that sequence, and each
base is checked against the reference on the way.

The reference used here has chr1 (600 bases) and chr2 (400 bases), so joined position 650
is position 50 of chr2.'''
import gzip
import os
import sqlite3
import subprocess
import sys

import pytest

import helpers
from CanSNPer2.DownloadGenomes import DownloadGenomes
from CanSNPer2.modules.ConvertPositions import (ConversionError, convert_database, convert_table,
	genome_files)
from CanSNPer2.modules.DatabaseConnection import XMFAFunctions, check_database

CHR1 = ("ACGT" * 150)                 ## 600 bases, position p has "ACGT"[(p-1) % 4]
CHR2 = ("TTGCA" * 80)                 ## 400 bases, position p has "TTGCA"[(p-1) % 5]
FASTA = ">chr1 first\n" + CHR1 + "\n>chr2 second\n" + CHR2 + "\n"
REVERSED = ">chr2 second\n" + CHR2 + "\n>chr1 first\n" + CHR1 + "\n"


def base(position):
	'''the base at a joined position of FASTA'''
	return (CHR1 + CHR2)[position - 1]


def write_reference(refdir, genome="REF", text=FASTA):
	os.makedirs(str(refdir), exist_ok=True)
	path = os.path.join(str(refdir), genome + ".fna")
	with open(path, "w") as out:
		out.write(text)
	return path


def write(path, text):
	with open(str(path), "w") as out:
		out.write(text)
	return str(path)


OLD_SNP = ("# a comment line\n"
	"snp_id\tstrain\treference\tgenome\tposition\tderived_base\tancestral_base\n"
	"A.1\tREF\ttest\tREF\t10\tT\t{b10}\n"
	"\n"
	"A.2\tREF\ttest\tREF\t650\t{b650}\tG\t# trailing comment\n").format(b10=base(10), b650=base(650))


def test_converts_snp_file_on_two_sequences(tmp_path):
	refs = {"REF": write_reference(tmp_path / "refs")}
	source = write(tmp_path / "old_snp.txt", OLD_SNP)
	target = str(tmp_path / "out" / "old_snp.txt")
	sequences = set()
	assert convert_table(source, target, refs, sequences=sequences) == 2
	assert sequences == {"chr1", "chr2"}
	with open(target) as fin:
		lines = fin.read().split("\n")
	assert lines[0] == "# a comment line"
	assert lines[1] == "snp_id\tstrain\treference\tgenome\tsequence\tsequence_position\tderived_base\tancestral_base"
	assert lines[2] == "A.1\tREF\ttest\tREF\tchr1\t10\tT\t" + base(10)
	assert lines[3] == ""
	assert lines[4] == "A.2\tREF\ttest\tREF\tchr2\t50\t" + base(650) + "\tG\t# trailing comment"


def test_converts_canSNPs_file_keeping_rank_and_score(tmp_path):
	refs = {"REF": write_reference(tmp_path / "refs")}
	source = write(tmp_path / "old_canSNPs.txt",
		"node\trank\tgenome\tposition\tderived_base\tancestral_base\tscore\n"
		"A.1\t1\tREF\t601\tG\t{0}\t400\n"
		"A.1\t2\tREF\t3\t{1}\tT\t350\n".format(base(601), base(3)))
	target = str(tmp_path / "new_canSNPs.txt")
	assert convert_table(source, target, refs) == 2
	with open(target) as fin:
		assert fin.read() == ("node\trank\tgenome\tsequence\tsequence_position\tderived_base\tancestral_base\tscore\n"
			"A.1\t1\tREF\tchr2\t1\tG\t{0}\t400\n"
			"A.1\t2\tREF\tchr1\t3\t{1}\tT\t350\n".format(base(601), base(3)))


def test_refuses_when_reference_is_in_another_order(tmp_path):
	'''with chr2 first, joined position 650 is chr1:250 and 10 is chr2:10, whose bases differ'''
	refs = {"REF": write_reference(tmp_path / "refs", text=REVERSED)}
	source = write(tmp_path / "old_snp.txt", OLD_SNP)
	target = str(tmp_path / "out" / "old_snp.txt")
	with pytest.raises(ConversionError) as error:
		convert_table(source, target, refs)
	message = str(error.value)
	assert "line 3" in message and "A.1" in message and "chr2:10" in message
	assert "line 5" in message and "A.2" in message and "650" in message and "chr1:250" in message
	assert not os.path.exists(target)


def test_refuses_a_position_beyond_the_end(tmp_path):
	refs = {"REF": write_reference(tmp_path / "refs")}
	source = write(tmp_path / "old_snp.txt",
		"snp_id\tstrain\treference\tgenome\tposition\tderived_base\tancestral_base\nA.1\tREF\ttest\tREF\t1001\tC\tA\n")
	target = str(tmp_path / "new.txt")
	with pytest.raises(ConversionError, match="beyond the end"):
		convert_table(source, target, refs)
	assert not os.path.exists(target)


def old_database(path):
	db = helpers.make_db(path, [("A.1", None, 10, base(10), "T"), ("A.2", None, 650, "G", base(650)),
		("A.3", None, 600, base(600), "A")], joined=True)
	conn = sqlite3.connect(db)
	conn.execute("UPDATE snp_references SET genbank_id = 'GCA_000000001.1', assembly_name = 'ASM1v1'")
	conn.commit()
	conn.close()
	return db


def test_converts_database_and_types_the_same(tmp_path):
	source = old_database(tmp_path / "old.db")
	refs = {"REF": write_reference(tmp_path / "refs")}
	target = str(tmp_path / "out" / "old.db")
	assert convert_database(source, target, refs) == 3
	assert check_database(target)
	conn = XMFAFunctions(target)
	snps, positions = conn.get_snps("REF")
	conn.disconnect()
	assert snps == {("chr1", 10): (base(10), "T", "A.1"), ("chr2", 50): ("G", base(650), "A.2"),
		("chr1", 600): (base(600), "A", "A.3")}
	assert positions == {"chr1": [10, 600], "chr2": [50]}
	with pytest.raises(Exception):
		check_database(source)  ## the old database is left as it was


def test_refuses_a_database_in_another_order(tmp_path):
	source = old_database(tmp_path / "old.db")
	refs = {"REF": write_reference(tmp_path / "refs", text=REVERSED)}
	target = str(tmp_path / "out" / "old.db")
	with pytest.raises(ConversionError, match="A.2"):
		convert_database(source, target, refs)
	assert not os.path.exists(target)


def test_refuses_to_overwrite(tmp_path):
	refs = {"REF": write_reference(tmp_path / "refs")}
	source = write(tmp_path / "old_snp.txt", OLD_SNP)
	with pytest.raises(ConversionError, match="is the input"):
		convert_table(source, source, refs)
	existing = write(tmp_path / "exists.txt", "keep\n")
	with pytest.raises(ConversionError, match="already exists"):
		convert_table(source, existing, refs)
	with open(existing) as fin:
		assert fin.read() == "keep\n"
	db = old_database(tmp_path / "old.db")
	with pytest.raises(ConversionError, match="is the input"):
		convert_database(db, db, refs)
	result = run_cli(tmp_path, "--convert", "--database", db, "--refdir", str(tmp_path / "refs"),
		"-o", str(tmp_path))
	assert result.returncode == 1
	assert "is the input" in result.stderr + result.stdout


def run_cli(tmp_path, *args):
	env = dict(os.environ, PYTHONPATH=os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
	return subprocess.run([sys.executable, "-m", "CanSNPer2.SNPDatabase"] + list(args) + ["--logs", str(tmp_path / "logs")],
		capture_output=True, text=True, cwd=str(tmp_path), env=env)


def test_cli_converts_a_snp_file_and_a_database(tmp_path):
	write_reference(tmp_path / "refs")
	source = write(tmp_path / "old_snp.txt", OLD_SNP)
	result = run_cli(tmp_path, "--convert", "--annotation", source, "--refdir", str(tmp_path / "refs"),
		"-o", str(tmp_path / "out"))
	assert result.returncode == 0, result.stderr
	assert "converted 2 SNPs of old_snp.txt to positions on 2 reference sequences" in result.stdout
	assert os.path.exists(str(tmp_path / "out" / "old_snp.txt"))
	db = old_database(tmp_path / "old.db")
	result = run_cli(tmp_path, "--convert", "--database", db, "--refdir", str(tmp_path / "refs"),
		"-o", str(tmp_path / "out"))
	assert result.returncode == 0, result.stderr
	assert check_database(str(tmp_path / "out" / "old.db"))
	result = run_cli(tmp_path, "--convert", "--database", db, "--refdir", str(tmp_path / "refs"))
	assert result.returncode == 1
	assert "--outdir" in result.stderr + result.stdout


def test_cli_names_one_reference_sequence_in_the_singular(tmp_path):
	write_reference(tmp_path / "refs")
	source = write(tmp_path / "old_snp.txt", "snp_id\tstrain\treference\tgenome\tposition\tderived_base\tancestral_base\n"
		"A.1\tREF\ttest\tREF\t10\tT\t{0}\n".format(base(10)))
	result = run_cli(tmp_path, "--convert", "--annotation", source, "--refdir", str(tmp_path / "refs"),
		"-o", str(tmp_path / "out"))
	assert result.returncode == 0, result.stderr
	assert "converted 1 SNP of old_snp.txt to positions on 1 reference sequence, written to" in result.stdout


def test_downloads_missing_reference_once(tmp_path, monkeypatch):
	fetched = []
	def retrieve(self, link, target, timeout=None):
		fetched.append(link)
		with open(target, "wb") as out:
			out.write(gzip.compress(FASTA.encode()))
	monkeypatch.setattr(DownloadGenomes, "retrieve", retrieve)
	refdir = str(tmp_path / "refs")
	references = {"REF": ("GCA_000000001.1", "ASM1v1")}
	files = genome_files(refdir, ["REF", "REF"], references)
	assert fetched == ["https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/000/001/GCA_000000001.1_ASM1v1/GCA_000000001.1_ASM1v1_genomic.fna.gz"]
	assert files == {"REF": os.path.join(refdir, "REF.fna")}
	assert genome_files(refdir, ["REF"], references) == files
	assert len(fetched) == 1
	with pytest.raises(ConversionError, match="OTHER"):
		genome_files(refdir, ["OTHER"], references)


def test_reads_accessions_from_the_references_file(tmp_path):
	'''the accessions of a genome missing from --refdir come from the --references file'''
	from CanSNPer2.modules.ConvertPositions import references_from_file
	path = write(tmp_path / "old_references.txt",
		"genome\tstrain\tgenbank_id\trefseq_id\tassembly_name\nREF\tREF strain\tGCA_000000001.1\tGCF_000000001.1\tASM1v1\n")
	assert references_from_file(path) == {"REF": ("GCA_000000001.1", "ASM1v1")}


# failures partway through leave no output behind ---------------------------------------

def test_a_failed_table_write_leaves_no_output(tmp_path, monkeypatch):
	import CanSNPer2.modules.ConvertPositions as ConvertPositions
	refs = {"REF": write_reference(tmp_path / "refs")}
	source = write(tmp_path / "old_snp.txt", OLD_SNP)
	def fail(src, dst):
		raise OSError("disk full")
	monkeypatch.setattr(ConvertPositions.os, "replace", fail)
	with pytest.raises(OSError, match="disk full"):
		convert_table(source, str(tmp_path / "out" / "old_snp.txt"), refs)
	assert os.listdir(str(tmp_path / "out")) == []


def test_a_failed_database_rebuild_leaves_no_output(tmp_path, monkeypatch):
	import CanSNPer2.modules.ConvertPositions as ConvertPositions
	source = old_database(tmp_path / "old.db")
	refs = {"REF": write_reference(tmp_path / "refs")}
	monkeypatch.setattr(ConvertPositions, "snp_annotation_sql", lambda table: "CREATE TABLE broken (")
	with pytest.raises(ConversionError, match="nothing was written"):
		convert_database(source, str(tmp_path / "out" / "old.db"), refs)
	assert os.listdir(str(tmp_path / "out")) == []


def test_an_existing_part_file_is_not_overwritten(tmp_path):
	refs = {"REF": write_reference(tmp_path / "refs")}
	source = write(tmp_path / "old_snp.txt", OLD_SNP)
	os.makedirs(str(tmp_path / "out"))
	part = write(tmp_path / "out" / "old_snp.txt.part", "keep\n")
	db_part = write(tmp_path / "out" / "old.db.part", "keep\n")
	assert convert_table(source, str(tmp_path / "out" / "old_snp.txt"), refs) == 2
	assert convert_database(old_database(tmp_path / "old.db"), str(tmp_path / "out" / "old.db"), refs) == 3
	for path in (part, db_part):
		with open(path) as fin:
			assert fin.read() == "keep\n"
	assert sorted(os.listdir(str(tmp_path / "out"))) == ["old.db", "old.db.part", "old_snp.txt", "old_snp.txt.part"]


def test_cli_reports_a_rebuilt_copy_that_fails_the_check(tmp_path, monkeypatch, caplog):
	'''the check of the rebuilt copy and sqlite errors end in one message and exit 1'''
	import CanSNPer2.modules.ConvertPositions as ConvertPositions
	from CanSNPer2 import SNPDatabase
	from CanSNPer2.modules.DatabaseConnection import ConnectionError as DBConnectionError
	db = old_database(tmp_path / "old.db")
	write_reference(tmp_path / "refs")
	real = ConvertPositions.check_database
	def check(path, current_format=True):
		if current_format:
			raise DBConnectionError("{0} is not a CanSNP database".format(path))
		return real(path, current_format=current_format)
	monkeypatch.setattr(ConvertPositions, "check_database", check)
	monkeypatch.setattr(sys, "argv", ["CanSNPer2-database", "--convert", "--database", db, "--refdir",
		str(tmp_path / "refs"), "-o", str(tmp_path / "out"), "--logs", str(tmp_path / "logs")])
	with pytest.raises(SystemExit) as stop:
		SNPDatabase.main()
	assert stop.value.code == 1
	errors = [record.getMessage() for record in caplog.records if record.levelname == "ERROR"]
	assert len(errors) == 1 and "failed the database check" in errors[0]
	assert os.listdir(str(tmp_path / "out")) == []


def test_cli_reports_a_file_that_is_not_a_database(tmp_path):
	write_reference(tmp_path / "refs")
	db = write(tmp_path / "old.db", "this is not an sqlite database\n" * 100)
	result = run_cli(tmp_path, "--convert", "--database", db, "--refdir", str(tmp_path / "refs"),
		"-o", str(tmp_path / "out"))
	assert result.returncode == 1
	assert "Traceback" not in result.stderr
	assert "not a database" in result.stderr
	assert not os.path.exists(str(tmp_path / "out" / "old.db"))


@pytest.mark.parametrize("umask", [0o022, 0o002, 0o077])
def test_outputs_get_the_mode_of_open(tmp_path, umask):
	'''the outputs are written through mkstemp, which gives mode 0600; a converted file is meant
		to be shared, so it gets the mode open() gives: 0666 masked by the umask'''
	refs = {"REF": write_reference(tmp_path / "refs")}
	source_table = write(tmp_path / "old_snp.txt", OLD_SNP)
	source_db = old_database(tmp_path / "old.db")
	previous = os.umask(umask)
	try:
		convert_table(source_table, str(tmp_path / "out" / "old_snp.txt"), refs)
		convert_database(source_db, str(tmp_path / "out" / "old.db"), refs)
	finally:
		os.umask(previous)
	for name in ("old_snp.txt", "old.db"):
		mode = os.stat(str(tmp_path / "out" / name)).st_mode & 0o777
		assert mode == 0o666 & ~umask, (name, oct(mode))
