"""Read typing with a local reference, end to end. The scheme reference is far from clade X in
part of the genome: reads of an X genome map to it only in the conserved half, so the first
round misses the canSNPs in the divergent half; the representative of X, lifted and mapped to
in a second round, gives them."""
import random
import shutil
import sqlite3

import pytest

import helpers
from CanSNPer2.modules.CanSNPer2 import CanSNPer2
from CanSNPer2.modules.ParsePileup import read_site_counts
from CanSNPer2.modules.SchemeGenomes import load_genomes

pytestmark = pytest.mark.skipif(shutil.which("minimap2") is None or shutil.which("samtools") is None,
	reason="minimap2 and samtools are needed for read typing")

BLOCK = 1000
# canSNPs: X.1, X.2, X.3 in conserved blocks, X.4, X.5 in divergent blocks; Y.1 on the other branch
SITES = {"X.1": 1500, "X.2": 3500, "X.3": 5500, "X.4": 8500, "X.5": 10500, "Y.1": 13500}


def mutate(seq, rate, rng, keep=()):
	out = list(seq)
	for i in range(len(out)):
		if i + 1 not in keep and rng.random() < rate:
			out[i] = rng.choice("ACGT".replace(out[i], ""))
	return out


def genomes(divergent):
	rng = random.Random(5)
	ref = [rng.choice("ACGT") for _ in range(16000)]
	keep = set(SITES.values())
	x = list(ref)
	for start in range(0, 16000, BLOCK):
		# odd blocks from 8 kb on are divergent, the rest 1 %
		rate = divergent if start >= 8000 and (start // BLOCK) % 2 == 0 else 0.01
		x[start:start + BLOCK] = mutate(ref[start:start + BLOCK], rate, rng, keep={p - start for p in keep if start < p <= start + BLOCK})
	y = mutate(ref, 0.01, rng, keep)
	derived = {}
	for name, pos in SITES.items():
		derived[name] = "A" if ref[pos - 1] != "A" else "C"
		if name.startswith("X"):
			x[pos - 1] = derived[name]
		else:
			y[pos - 1] = derived[name]
	sample = mutate(x, 0.002, rng, keep)  # a genome of clade X, close to its representative
	return "".join(ref), "".join(x), "".join(y), "".join(sample), derived


def make_scheme(tmp_path, divergent=0.15):
	ref, x, y, sample, derived = genomes(divergent)
	db = str(tmp_path / "toy.db")
	conn = sqlite3.connect(db)
	conn.executescript(helpers.SNP_ANNOTATION_TABLE + helpers.REFERENCES_TABLE + """
		CREATE TABLE nodes (id integer PRIMARY KEY, name text NOT NULL);
		CREATE TABLE tree (parent integer, child integer, rank_i integer, unique (parent, child, rank_i));
		INSERT INTO snp_references (id, genome, strain) VALUES (5, 'REF', 'ref');
		INSERT INTO nodes (id, name) VALUES (1, 'root'), (2, 'X.1'), (3, 'X.2'), (4, 'X.3'), (5, 'X.4'), (6, 'X.5'), (7, 'Y.1');
		INSERT INTO tree (parent, child, rank_i) VALUES (1, 1, 1), (1, 2, 1), (2, 3, 1), (3, 4, 1), (4, 5, 1), (5, 6, 1), (1, 7, 1);""")
	for node, name in enumerate(["X.1", "X.2", "X.3", "X.4", "X.5", "Y.1"], start=2):
		pos = SITES[name]
		conn.execute("INSERT INTO snp_annotation (node_id, snp_id, sequence, sequence_position, ancestral_base, derived_base,"
			" reference, date, genome_i) VALUES (?, ?, 'chr', ?, ?, ?, 't', '24/11/26', 5)", (node, name, pos, ref[pos - 1], derived[name]))
	conn.commit()
	(tmp_path / "genomes.txt").write_text("genome\tnode\trole\nLX\tX.5\trepresentative\nLY\tY.1\trepresentative\n")
	load_genomes(conn, str(tmp_path / "genomes.txt"))
	conn.close()
	local = tmp_path / "refs" / "local"
	local.mkdir(parents=True)
	(tmp_path / "refs" / "REF.fna").write_text(">chr\n" + ref + "\n")
	(local / "LX.fna").write_text(">x1\n" + x[:9000] + "\n>x2\n" + x[9000:] + "\n")
	(local / "LY.fna").write_text(">y1\n" + y + "\n")
	r1, r2 = helpers.simulate_reads(sample, 40, seed=21)
	(tmp_path / "s_R1.fq").write_text(r1)
	(tmp_path / "s_R2.fq").write_text(r2)
	return db, [str(tmp_path / "s_R1.fq"), str(tmp_path / "s_R2.fq")]


def run(tmp_path, files, db, out, **kwargs):
	obj = CanSNPer2(files, workdir=str(tmp_path), refdir="refs", tmpdir="tmp", logdir="logs", outdir=out,
		database=db, export=True, **kwargs)
	obj.run(database=db)
	return obj


def comments(path):
	return [line[2:].rstrip("\n") for line in open(path) if line.startswith("# ")]


def test_the_first_round_alone_misses_the_divergent_canSNPs(tmp_path):
	db, files = make_scheme(tmp_path)
	obj = run(tmp_path, files, db, "plain")
	assert obj.typed == {"s": "X.3"}
	rows = {r["Name"]: r for r in read_site_counts(str(tmp_path / "plain" / "samples" / "s_site_counts.tsv"))}
	assert rows["X.4"]["state"] in ("no_reads", "low_depth") and rows["X.5"]["state"] in ("no_reads", "low_depth")


def test_a_local_reference_reaches_the_leaf(tmp_path):
	db, files = make_scheme(tmp_path)
	obj = run(tmp_path, files, db, "local", local_references=True)
	assert obj.typed == {"s": "X.5"}
	out = tmp_path / "local" / "samples"
	rows = {r["Name"]: r for r in read_site_counts(str(out / "s_site_counts.tsv"))}
	assert rows["X.5"]["state"] == "derived" and rows["X.5"]["source"] == "LX"
	assert rows["X.5"]["lifted_position"].startswith("x2:")
	reads = [l for l in comments(out / "s_snps.txt") if l.startswith("Reads:")][0]
	assert "local_reference=LX local_node=X.5" in reads and "conflicts=0" in reads
	settings = [l for l in comments(out / "s_snps.txt") if l.startswith("Typing settings:")][0]
	assert "local_references=on" in settings
	assert any(p.name.endswith(".lift.tsv") for p in (tmp_path / "refs" / "local" / "lift").iterdir()), "the lift is cached"


def test_without_local_genome_files_the_first_round_stands(tmp_path, caplog):
	db, files = make_scheme(tmp_path)
	for f in (tmp_path / "refs" / "local").iterdir():
		f.unlink()
	obj = run(tmp_path, files, db, "nolocal", local_references=True)
	assert obj.typed == {"s": "X.3"}
	assert "CanSNPer2-download" in caplog.text and "--local_references" in caplog.text
	reads = [l for l in comments(tmp_path / "nolocal" / "samples" / "s_snps.txt") if l.startswith("Reads:")][0]
	assert "local_reference=none" in reads


def test_a_sample_close_to_the_scheme_reference_keeps_its_call(tmp_path):
	db, _ = make_scheme(tmp_path)
	_ref, _x, y, _s, _d = genomes(0.15)
	r1, r2 = helpers.simulate_reads(y, 40, seed=22)
	(tmp_path / "y_R1.fq").write_text(r1)
	(tmp_path / "y_R2.fq").write_text(r2)
	files = [str(tmp_path / "y_R1.fq"), str(tmp_path / "y_R2.fq")]
	plain = run(tmp_path, files, db, "yplain", min_required_hits=1)
	local = run(tmp_path, files, db, "ylocal", min_required_hits=1, local_references=True)
	assert plain.typed == local.typed == {"y": "Y.1"}


def test_a_local_genome_that_does_not_fit_its_node_is_not_used(tmp_path, caplog):
	db, files = make_scheme(tmp_path)
	# the file named LX holds the Y genome: its canSNP bases type it to Y.1, not to X.5
	(tmp_path / "refs" / "local" / "LX.fna").write_text((tmp_path / "refs" / "local" / "LY.fna").read_text())
	obj = run(tmp_path, files, db, "wrong", local_references=True)
	assert obj.typed == {"s": "X.3"}, "the first round stands"
	assert "do not type it to its own node" in caplog.text
