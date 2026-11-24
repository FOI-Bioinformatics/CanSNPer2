"""Read typing end to end: reads simulated from toy genomes are mapped, piled up and typed
through the whole run, with the result files of assembly typing plus the site counts."""
import gzip
import logging
import random
import shutil
import sqlite3

import pytest

import helpers
from CanSNPer2.modules.CanSNPer2 import CanSNPer2, CanSNPer2Error
from CanSNPer2.modules.ParsePileup import read_site_counts

pytestmark = pytest.mark.skipif(shutil.which("minimap2") is None or shutil.which("samtools") is None,
	reason="minimap2 and samtools are needed for read typing")

# root -> A.1 -> A.2 -> A.3 and root -> B.1 -> B.2 -> B.3, one canSNP each
SITES = {"A.1": 2000, "A.2": 5000, "A.3": 8000, "B.1": 11000, "B.2": 14000, "B.3": 17000}


def genomes():
	rng = random.Random(3)
	ref = "".join(rng.choice("ACGT") for _ in range(20000))
	derived = {name: ("A" if ref[pos - 1] != "A" else "C") for name, pos in SITES.items()}
	def with_snps(prefix):
		seq = list(ref)
		for name, pos in SITES.items():
			if name.startswith(prefix):
				seq[pos - 1] = derived[name]
		return "".join(seq)
	return ref, with_snps("A"), with_snps("B"), derived


def make_scheme(tmp_path):
	ref, a, b, derived = genomes()
	db = str(tmp_path / "toy.db")
	conn = sqlite3.connect(db)
	conn.executescript(helpers.SNP_ANNOTATION_TABLE + helpers.REFERENCES_TABLE + """
		CREATE TABLE nodes (id integer PRIMARY KEY, name text NOT NULL);
		CREATE TABLE tree (parent integer, child integer, rank_i integer, unique (parent, child, rank_i));
		INSERT INTO snp_references (id, genome, strain) VALUES (5, 'REF', 'ref strain');
		INSERT INTO nodes (id, name) VALUES (1, 'root'), (2, 'A.1'), (3, 'A.2'), (4, 'A.3'), (5, 'B.1'), (6, 'B.2'), (7, 'B.3');
		INSERT INTO tree (parent, child, rank_i) VALUES (1, 1, 1), (1, 2, 1), (2, 3, 1), (3, 4, 1), (1, 5, 1), (5, 6, 1), (6, 7, 1);""")
	for node, name in enumerate(["A.1", "A.2", "A.3", "B.1", "B.2", "B.3"], start=2):
		pos = SITES[name]
		conn.execute("INSERT INTO snp_annotation (node_id, snp_id, sequence, sequence_position, ancestral_base, derived_base,"
			" reference, date, genome_i) VALUES (?, ?, 'chr', ?, ?, ?, 'test', '24/11/26', 5)", (node, name, pos, ref[pos - 1], derived[name]))
	conn.commit()
	conn.close()
	refs = tmp_path / "refs"
	refs.mkdir(exist_ok=True)
	(refs / "REF.fna").write_text(">chr\n" + ref + "\n")
	return db, a, b


def write_pair(tmp_path, name, parts, gz=True):
	'''parts: (genome, depth, seed); the reads of all parts in one pair of files'''
	r1 = "".join(helpers.simulate_reads(g, d, seed=s)[0] for g, d, s in parts)
	r2 = "".join(helpers.simulate_reads(g, d, seed=s)[1] for g, d, s in parts)
	paths = []
	for mate, text in ((1, r1), (2, r2)):
		path = tmp_path / "{n}_R{m}_001.fastq{gz}".format(n=name, m=mate, gz=".gz" if gz else "")
		if gz:
			with gzip.open(path, "wt") as out:
				out.write(text)
		else:
			path.write_text(text)
		paths.append(str(path))
	return paths


def run(tmp_path, queries, **kwargs):
	db = kwargs.pop("db")
	obj = CanSNPer2(queries, workdir=str(tmp_path), refdir="refs", tmpdir="tmp", logdir="logs", outdir="out",
		database=db, export=True, **kwargs)
	obj.run(database=db)
	return obj


def comments(path):
	return [line[2:].rstrip("\n") for line in open(path) if line.startswith("# ")]


def test_paired_reads_are_typed_like_the_genome(tmp_path, capsys):
	db, a, _ = make_scheme(tmp_path)
	files = write_pair(tmp_path, "sa", [(a, 30, 1)])
	obj = run(tmp_path, files, db=db)
	assert obj.typed == {"sa": "A.3"}
	assert "sa: A.3" in capsys.readouterr().out
	out = tmp_path / "out" / "samples"
	header = comments(out / "sa_snps.txt")
	assert any(line.startswith("Reads: ") and "platform=illumina" in line and "mixed_sites=0" in line for line in header)
	reads_line = [line for line in header if line.startswith("Reads: ")][0]
	assert "sites_with_reads=6/6" in reads_line
	assert any("platform=illumina min_depth=5" in line for line in header if line.startswith("Typing settings:"))
	rows = read_site_counts(str(out / "sa_site_counts.tsv"))
	assert {r["Name"]: r["state"] for r in rows} == {"A.1": "derived", "A.2": "derived", "A.3": "derived",
		"B.1": "ancestral", "B.2": "ancestral", "B.3": "ancestral"}
	assert all(int(r["depth"]) >= 5 for r in rows)


def test_the_median_depth_counts_only_sites_with_reads(tmp_path):
	# reads of the first half of the genome only: the B sites (11-17 kb) get none
	db, a, _ = make_scheme(tmp_path)
	obj = run(tmp_path, write_pair(tmp_path, "half", [(a[:10000], 30, 15)]), db=db)
	header = comments(tmp_path / "out" / "samples" / "half_snps.txt")
	reads_line = [line for line in header if line.startswith("Reads: ")][0]
	assert "sites_with_reads=3/6" in reads_line
	depth = int(reads_line.split("median_depth=")[1].split()[0])
	assert depth >= 15, "the median of the three sites with reads, not of all six"


def test_too_few_reads_give_no_call_not_a_wrong_one(tmp_path):
	db, a, _ = make_scheme(tmp_path)
	obj = run(tmp_path, write_pair(tmp_path, "thin", [(a, 1, 2)]), db=db)
	assert obj.typed == {} and "thin" in obj.not_typed


def test_an_even_mixture_is_flagged(tmp_path):
	db, a, b = make_scheme(tmp_path)
	obj = run(tmp_path, write_pair(tmp_path, "mix", [(a, 20, 3), (b, 20, 4)]), db=db)
	report = obj.sample_reports["mix"]
	assert any("possibly mixed" in reason and "duplicated region" in reason for reason in report["reasons"])
	assert report["severity"] != "clean"
	header = comments(tmp_path / "out" / "samples" / "mix_snps.txt")
	assert any("mixed_sites=6" in line for line in header)


def test_a_minor_strain_below_the_fraction_is_not_flagged(tmp_path):
	db, a, b = make_scheme(tmp_path)
	obj = run(tmp_path, write_pair(tmp_path, "most", [(a, 40, 5), (b, 2, 6)]), db=db)
	assert obj.typed == {"most": "A.3"}
	assert not any("possibly mixed" in reason for reason in obj.sample_reports["most"]["reasons"])


def test_long_reads_from_a_samplesheet(tmp_path):
	db, _, b = make_scheme(tmp_path)
	reads, _ = helpers.simulate_reads(b, 30, seed=7, length=2000, paired=False)
	(tmp_path / "long.fq").write_text(reads)
	sheet = tmp_path / "samples.tsv"
	sheet.write_text("sample\treads1\tplatform\nlr\tlong.fq\tont\n")
	obj = run(tmp_path, [], db=db, samplesheet=str(sheet))
	assert obj.typed == {"lr": "B.3"}
	settings = [l for l in comments(tmp_path / "out" / "samples" / "lr_snps.txt") if l.startswith("Typing settings:")]
	assert "platform=ont min_depth=5 min_allele_fraction=0.7" in settings[0]


def test_alignments_are_kept_and_read_again(tmp_path):
	db, a, _ = make_scheme(tmp_path)
	files = write_pair(tmp_path, "keep", [(a, 30, 8)])
	first = run(tmp_path, files, db=db, keep_temp=True)
	again = CanSNPer2(files, workdir=str(tmp_path), refdir="refs", tmpdir=first.tmpdir, logdir="logs", outdir="out2",
		database=db, export=True, skip_alignment=True)
	again.run(database=db)
	assert again.typed == first.typed == {"keep": "A.3"}


def test_reads_and_assemblies_are_not_mixed_in_one_run(tmp_path):
	db, a, _ = make_scheme(tmp_path)
	files = write_pair(tmp_path, "r", [(a, 5, 9)])
	(tmp_path / "g.fna").write_text(">c\n" + a + "\n")
	with pytest.raises(CanSNPer2Error, match="reads and assemblies"):
		run(tmp_path, files + [str(tmp_path / "g.fna")], db=db)


def test_an_uncompressed_fastq_is_not_taken_for_a_query_list(tmp_path):
	db, a, _ = make_scheme(tmp_path)
	reads, _ = helpers.simulate_reads(a, 30, seed=10, paired=False, length=150)
	(tmp_path / "single.fastq").write_text(reads)
	obj = run(tmp_path, [str(tmp_path / "single.fastq")], db=db)
	assert obj.typed == {"single": "A.3"}


def test_mask_does_not_apply_to_reads(tmp_path, caplog):
	db, a, _ = make_scheme(tmp_path)
	with caplog.at_level(logging.WARNING):
		obj = run(tmp_path, write_pair(tmp_path, "m", [(a, 30, 11)]), db=db, mask=10)
	assert obj.typed == {"m": "A.3"}
	assert "--mask" in caplog.text


def cli(tmp_path, *args):
	import subprocess, sys
	return subprocess.run([sys.executable, "-m", "CanSNPer2.CanSNPerTree", *args, "--workdir", str(tmp_path),
		"--refdir", "refs", "--tmpdir", "tmp", "-o", "out"], capture_output=True, text=True, cwd=str(tmp_path))


def test_command_line_reads(tmp_path):
	db, a, _ = make_scheme(tmp_path)
	files = write_pair(tmp_path, "c1", [(a, 30, 12)])
	run = cli(tmp_path, "--database", db, *files, "--min_depth", "4", "--threads", "2")
	assert run.returncode == 0, run.stderr
	assert "c1: A.3" in run.stdout
	settings = [l for l in comments(tmp_path / "out" / "samples" / "c1_snps.txt") if l.startswith("Typing settings:")][0]
	assert "min_depth=4" in settings


def test_command_line_samplesheet_without_query_arguments(tmp_path):
	db, _, b = make_scheme(tmp_path)
	reads, _ = helpers.simulate_reads(b, 30, seed=13, length=2000, paired=False)
	(tmp_path / "hifi.fq").write_text(reads)
	(tmp_path / "s.tsv").write_text("sample\treads1\nh1\thifi.fq\n")
	run = cli(tmp_path, "--database", db, "--samplesheet", str(tmp_path / "s.tsv"), "--platform", "hifi")
	assert run.returncode == 0, run.stderr
	assert "h1: B.3" in run.stdout


def test_read_input_flag_is_accepted(tmp_path):
	db, a, _ = make_scheme(tmp_path)
	files = write_pair(tmp_path, "c2", [(a, 30, 14)])
	run = cli(tmp_path, "--database", db, "--read_input", *files)
	assert run.returncode == 0, run.stderr
	assert "c2: A.3" in run.stdout


def test_unknown_platform_is_an_argument_error(tmp_path):
	db, _, _ = make_scheme(tmp_path)
	run = cli(tmp_path, "--database", db, "--platform", "sanger", "x.fq")
	assert run.returncode == 2 and "platform" in run.stderr
