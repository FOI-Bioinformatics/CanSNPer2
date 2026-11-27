"""db export --format vcf (#98): the canSNPs of a database as sites-only VCF files, one per
reference genome, with the REF base read from the reference file."""

import glob
import os
import shutil
import sqlite3
import subprocess

import pytest

import helpers
from CanSNPer2.modules.ExportVCF import encode
from test_modify_remove import cansnper2

SNP_HEADER = "\t".join(helpers.CANSNP_HEADER)
TREE = "R\nR\tA.1\nR\tA.1\tA.2\nR\tB.1\n"
## node, rank, genome, sequence, position, derived, ancestral, typed; the reference files below
## carry the ancestral base, except at A.1, where REF1 carries the derived base
CANSNPS = [
    ("R", 1, "REF1", "chrA", 10, "T", "C", "yes"),
    ("A.1", 1, "REF1", "chrA", 20, "G", "A", "yes"),
    ("A.1", 2, "REF1", "chrA", 30, "A", "G", "no"),
    ("A.2", 1, "REF1", "chrB", 5, "T", "C", "yes"),
    ("B.1", 1, "REF2", "contig_1", 7, "A", "G", "yes"),
]
REFERENCES = (
    "genome\tstrain\tgenbank_id\trefseq_id\tassembly_name\n"
    "REF1\ts1\tGCA_000000001.1\tGCF_000000001.1\tASM1v1\n"
    "REF2\ts2\tGCA_000000002.1\t\tASM2v1\n"
)
LENGTHS = {"REF1": {"chrA": 60, "chrB": 30}, "REF2": {"contig_1": 40}}


@pytest.fixture(autouse=True)
def in_tmp_path(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)


def reference_files(folder, override=None):
    """REF1.fna and REF2.fna with the ancestral base at every canSNP (the derived one at A.1)"""
    os.makedirs(folder, exist_ok=True)
    seqs = {g: {s: ["T"] * n for s, n in contigs.items()} for g, contigs in LENGTHS.items()}
    for node, rank, genome, seq, pos, derived, ancestral, _typed in CANSNPS:
        seqs[genome][seq][pos - 1] = derived if (node, rank) == ("A.1", 1) else ancestral
    for (genome, seq, pos), base in (override or {}).items():
        seqs[genome][seq][pos - 1] = base
    for genome, contigs in seqs.items():
        with open(os.path.join(folder, genome + ".fna"), "w") as out:
            for name, bases in contigs.items():
                out.write(">{name} test sequence\n{seq}\n".format(name=name, seq="".join(bases)))


@pytest.fixture
def scheme(tmp_path):
    rows = ["{0}\t{1}\tPaper2020\t{2}\t{3}\t{4}\t{5}\t{6}\t\t{7}".format(*row) for row in CANSNPS]
    files = helpers.write_sources(tmp_path / "src", rows, tree=TREE, references=REFERENCES, header=SNP_HEADER)
    db = tmp_path / "s.db"
    status = cansnper2(
        "db", "create", "-d", db, "--tree", files["tree"], "--cansnps", files["snp"],
        "--references", files["references"],
    )  # fmt: skip
    assert status == 0
    reference_files(str(tmp_path / "refs"))
    return db


def records(path):
    lines = open(path).read().splitlines()
    header = [line for line in lines if line.startswith("##")]
    columns = [line for line in lines if line.startswith("#CHROM")]
    body = [line.split("\t") for line in lines if not line.startswith("#")]
    return header, columns, body


def test_one_vcf_per_reference_genome_with_ref_and_alt(scheme, tmp_path):
    assert cansnper2("db", "export", "-d", scheme, "--format", "vcf", "--ref-dir", "refs", "-o", "vcf") == 0
    assert sorted(os.listdir(tmp_path / "vcf")) == ["s_REF1.vcf", "s_REF2.vcf"]
    header, columns, body = records(tmp_path / "vcf" / "s_REF1.vcf")
    assert header[0] == "##fileformat=VCFv4.3"
    assert "##contig=<ID=chrA,length=60>" in header and "##contig=<ID=chrB,length=30>" in header
    assert header.index("##contig=<ID=chrA,length=60>") < header.index("##contig=<ID=chrB,length=30>")
    assert "##cansnper2_reference_genome=REF1 (GCA_000000001.1)" in header
    assert any(line.startswith("##cansnper2_database=s.db sha256:") for line in header)
    assert columns == ["#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"]
    assert body == [
        ["chrA", "10", ".", "C", "T", ".", ".", "AA=C;DA=T;NODE=R;RANK=1;TYPED"],
        ["chrA", "20", ".", "G", "A", ".", ".", "AA=A;DA=G;NODE=A.1;RANK=1;TYPED"],
        ["chrA", "30", ".", "G", "A", ".", ".", "AA=G;DA=A;NODE=A.1;RANK=2"],
        ["chrB", "5", ".", "C", "T", ".", ".", "AA=C;DA=T;NODE=A.2;RANK=1;TYPED"],
    ]
    _h, _c, body = records(tmp_path / "vcf" / "s_REF2.vcf")
    assert body == [["contig_1", "7", ".", "G", "A", ".", ".", "AA=G;DA=A;NODE=B.1;RANK=1;TYPED"]]


def test_every_cansnp_of_the_database_is_in_the_vcf_files(scheme, tmp_path):
    assert cansnper2("db", "export", "-d", scheme, "--format", "vcf", "--ref-dir", "refs", "-o", "vcf") == 0
    exported = set()
    for path in glob.glob(str(tmp_path / "vcf" / "*.vcf")):
        for row in records(path)[2]:
            info = dict(field.split("=") for field in row[7].split(";") if "=" in field)
            exported.add((info["NODE"], int(info["RANK"]), row[0], int(row[1])))
    conn = sqlite3.connect(str(scheme))
    stored = set(conn.execute("SELECT snp_id, rank, sequence, sequence_position FROM snp_annotation"))
    conn.close()
    assert exported == stored


def test_a_missing_reference_file_exports_nothing(scheme, tmp_path):
    os.remove(tmp_path / "refs" / "REF2.fna")
    assert cansnper2("db", "export", "-d", scheme, "--format", "vcf", "--ref-dir", "refs", "-o", "vcf") == 1
    assert not os.path.exists(tmp_path / "vcf") or os.listdir(tmp_path / "vcf") == []


def test_a_missing_reference_file_names_ref_get(scheme, tmp_path, caplog):
    os.remove(tmp_path / "refs" / "REF2.fna")
    cansnper2("db", "export", "-d", scheme, "--format", "vcf", "--ref-dir", "refs", "-o", "vcf")
    assert "REF2" in caplog.text and "cansnper2 ref get" in caplog.text


def test_a_reference_base_that_is_neither_allele_is_refused(scheme, tmp_path, caplog):
    reference_files(str(tmp_path / "refs"), override={("REF1", "chrB", 5): "A"})
    assert cansnper2("db", "export", "-d", scheme, "--format", "vcf", "--ref-dir", "refs", "-o", "vcf") == 1
    assert "A.2" in caplog.text and "neither the ancestral C nor the derived T" in caplog.text
    assert not os.path.exists(tmp_path / "vcf") or os.listdir(tmp_path / "vcf") == []


def test_ref_dir_is_refused_with_another_format(scheme, tmp_path, caplog):
    assert cansnper2("db", "export", "-d", scheme, "--ref-dir", "refs", "-o", "src_out") == 1
    assert "--ref-dir is read by --format vcf only" in caplog.text


def test_existing_vcf_files_are_not_overwritten(scheme, tmp_path):
    os.makedirs(tmp_path / "vcf")
    (tmp_path / "vcf" / "s_REF2.vcf").write_text("mine\n")
    assert cansnper2("db", "export", "-d", scheme, "--format", "vcf", "--ref-dir", "refs", "-o", "vcf") == 1
    assert (tmp_path / "vcf" / "s_REF2.vcf").read_text() == "mine\n"
    assert not (tmp_path / "vcf" / "s_REF1.vcf").exists()


@pytest.mark.parametrize(
    "value,expected", [("T/N.1", "T/N.1"), ("A;B", "A%3BB"), ("a=b,c", "a%3Db%2Cc"), ("x y", "x%20y")]
)
def test_info_values_are_percent_encoded(value, expected):
    assert encode(value) == expected


def bcftools():
    found = shutil.which("bcftools") or os.environ.get("BCFTOOLS")
    return found if found and os.path.exists(found) else None


@pytest.mark.skipif(bcftools() is None, reason="bcftools is not installed")
def test_bcftools_reads_the_vcf(scheme, tmp_path):
    assert cansnper2("db", "export", "-d", scheme, "--format", "vcf", "--ref-dir", "refs", "-o", "vcf") == 0
    out = subprocess.run(
        [
            bcftools(),
            "query",
            "-f",
            "%CHROM\t%POS\t%REF\t%ALT\t%INFO/NODE\t%INFO/TYPED\n",
            str(tmp_path / "vcf" / "s_REF1.vcf"),
        ],
        capture_output=True,
        text=True,
    )
    assert out.returncode == 0, out.stderr
    assert out.stdout.splitlines() == [
        "chrA\t10\tC\tT\tR\t1",
        "chrA\t20\tG\tA\tA.1\t1",
        "chrA\t30\tG\tA\tA.1\t.",  # bcftools prints an absent flag as .
        "chrB\t5\tC\tT\tA.2\t1",
    ]
    check = subprocess.run(
        [
            bcftools(),
            "norm",
            "--check-ref",
            "e",
            "-f",
            str(tmp_path / "refs" / "REF1.fna"),
            str(tmp_path / "vcf" / "s_REF1.vcf"),
        ],
        capture_output=True,
        text=True,
    )
    assert check.returncode == 0, check.stderr
