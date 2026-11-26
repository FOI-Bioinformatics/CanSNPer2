"""Lifting canSNP positions from the scheme reference to a local reference genome: the local
genome is aligned to the reference with minimap2 asm20, and each canSNP gets the contig,
position and strand of the local genome and the local base (LiftOver)."""

import random
import shutil

import pytest

from CanSNPer2.modules.LiftOver import Lifted, lift, lift_snps

pytestmark = pytest.mark.skipif(shutil.which("minimap2") is None, reason="minimap2 is not installed")

COMP = str.maketrans("ACGT", "TGCA")


def rc(s):
    return s.translate(COMP)[::-1]


def other(base, *avoid):
    return next(b for b in "ACGT" if b != base and b not in avoid)


def genomes():
    """a 20 kb reference and a local genome 2 % different, in three contigs: c1 = local 1-7000,
    c2 = reverse complement of local 7001-14000, c3 = local 14501-20000 (14001-14500 deleted)"""
    rng = random.Random(11)
    ref = [rng.choice("ACGT") for _ in range(20000)]
    local = list(ref)
    sites = {1000, 5000, 9000, 12000, 14200, 17000}
    for i in range(20000):
        if i + 1 not in sites and rng.random() < 0.02:
            local[i] = other(local[i])
    snps = []
    for snp_id, pos, state in (
        ("S1", 1000, "derived"),
        ("S2", 5000, "ancestral"),
        ("S3", 9000, "derived"),
        ("S4", 12000, "ancestral"),
        ("S5", 14200, "derived"),
        ("S6", 17000, "third"),
    ):
        anc = ref[pos - 1]
        der = other(anc)
        if state == "derived":
            local[pos - 1] = der
        elif state == "third":
            local[pos - 1] = other(anc, der)
        snps.append((snp_id, "chr", pos, anc, der))
    ref, local = "".join(ref), "".join(local)
    contigs = {"c1": local[:7000], "c2": rc(local[7000:14000]), "c3": local[14500:]}
    return ref, contigs, snps


def write(tmp_path):
    ref, contigs, snps = genomes()
    (tmp_path / "ref.fna").write_text(">chr\n" + ref + "\n")
    (tmp_path / "local.fna").write_text("".join(">{0}\n{1}\n".format(k, v) for k, v in contigs.items()))
    return str(tmp_path / "ref.fna"), str(tmp_path / "local.fna"), snps


def test_lifted_positions_strands_and_bases(tmp_path):
    ref, local, snps = write(tmp_path)
    lifted = lift_snps(ref, local, snps, str(tmp_path))
    anc = {s[0]: s[3] for s in snps}
    der = {s[0]: s[4] for s in snps}
    assert lifted["S1"] == Lifted("c1", 1000, "+", der["S1"])
    assert lifted["S2"] == Lifted("c1", 5000, "+", anc["S2"])
    assert lifted["S3"] == Lifted("c2", 14000 - 9000 + 1, "-", der["S3"]), "reverse complemented contig"
    assert lifted["S4"] == Lifted("c2", 14000 - 12000 + 1, "-", anc["S4"])
    assert "S5" not in lifted, "the local genome lacks the region"
    assert "S6" not in lifted, "a local base that is neither allele is not lifted"


def test_the_local_base_is_the_base_of_the_local_contig(tmp_path):
    ref, local, snps = write(tmp_path)
    contigs = {}
    name = None
    for line in open(local):
        if line.startswith(">"):
            name = line[1:].strip()
            contigs[name] = ""
        else:
            contigs[name] += line.strip()
    for snp_id, lf in lift_snps(ref, local, snps, str(tmp_path)).items():
        base = contigs[lf.contig][lf.position - 1]
        assert (base if lf.strand == "+" else base.translate(COMP)) == lf.base, snp_id


def test_the_lift_is_cached_and_made_again_when_a_genome_changes(tmp_path):
    ref, local, snps = write(tmp_path)
    cache = tmp_path / "cache"
    first = lift(ref, local, snps, str(cache), str(tmp_path))
    files = list(cache.iterdir())
    assert len(files) == 1 and files[0].name.endswith(".lift.tsv")
    stamp = files[0].stat().st_mtime_ns
    assert lift(ref, local, snps, str(cache), str(tmp_path)) == first
    assert files[0].stat().st_mtime_ns == stamp, "read from the cache, not made again"
    with open(local, "a") as out:
        out.write(">c4\n" + "ACGT" * 50 + "\n")
    assert lift(ref, local, snps, str(cache), str(tmp_path)) == first
    assert files[0].stat().st_mtime_ns != stamp, "the local genome changed, the lift was made again"


def test_a_snp_added_to_the_scheme_makes_the_lift_again(tmp_path):
    ref, local, snps = write(tmp_path)
    cache = tmp_path / "cache"
    lift(ref, local, snps[:3], str(cache), str(tmp_path))
    assert "S4" in lift(ref, local, snps, str(cache), str(tmp_path))
