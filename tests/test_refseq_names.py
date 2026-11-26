"""RefSeq downloads are renamed to the GenBank names the database uses.

A canSNP is stored on a sequence named by its GenBank accession. The RefSeq assembly
(GCF) names its sequences NC_ or NZ_, so the sequences are renamed with the NCBI
assembly report before the reference check."""

import gzip
import os

import pytest

from test_reference_sequences import CHR1, CHR2, PLASMID, Results, fasta, snp_rows
from CanSNPer2.DownloadGenomes import DownloadGenomes, assembly_report_url, genbank_names, rename_sequences

REPORT = (
    "# Assembly name:  ASM712v1\n"
    "# Sequence-Name\tSequence-Role\tAssigned-Molecule\tAssigned-Molecule-Location/Type\tGenBank-Accn\tRelationship\tRefSeq-Accn\tAssembly-Unit\tSequence-Length\tUCSC-style-name\n"
    "chrI\tassembled-molecule\tI\tChromosome\tchrI\t=\tNZ_CP007762.1\tPrimary Assembly\t600\tna\n"
    "chrII\tassembled-molecule\tII\tChromosome\tchrII\t=\tNZ_CP007763.1\tPrimary Assembly\t400\tna\n"
    "pX\tassembled-molecule\tpX\tPlasmid\tna\t<>\tNZ_CP007764.1\tPrimary Assembly\t200\tna\n"
)
GENOME = fasta(("NZ_CP007762.1", CHR1), ("NZ_CP007763.1", CHR2), ("NZ_CP007764.1", PLASMID))


def make(tmp_path, genome, report=REPORT, fail_report=False):
    obj = DownloadGenomes.__new__(DownloadGenomes)
    obj.source, obj.directory, obj.database = "refseq", str(tmp_path / "refs"), None
    fetched = []

    def retrieve(link, target, timeout=None):
        fetched.append(link)
        if link.endswith("_assembly_report.txt"):
            if fail_report:
                raise OSError("report not found")
            with open(target, "w") as out:
                out.write(report)
        else:
            with open(target, "wb") as out:
                out.write(gzip.compress(genome.encode()))

    obj.retrieve = retrieve
    return obj, fetched


def test_assembly_report_url():
    url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/007/125/GCF_000007125.1_ASM712v1/GCF_000007125.1_ASM712v1_genomic.fna.gz"
    assert assembly_report_url(url) == url[: -len("_genomic.fna.gz")] + "_assembly_report.txt"


def test_genbank_names_skips_na():
    names = genbank_names(REPORT)
    assert names == {"NZ_CP007762.1": "chrI", "NZ_CP007763.1": "chrII"}


def test_refseq_download_is_renamed_and_passes_check(tmp_path):
    obj, fetched = make(tmp_path, GENOME)
    results = Results()
    obj.download(
        "GCF_000007125.1",
        "ASM712v1",
        "REF",
        source="refseq",
        directory=obj.directory,
        results=results,
        snps=snp_rows(),
    )
    assert results[0][1] is True, results
    text = open(os.path.join(obj.directory, "REF.fna")).read()
    assert text.startswith(">chrI test\n") and ">chrII test\n" in text
    assert ">pX" not in text and ">NZ_CP007764.1 test\n" in text
    assert any(link.endswith("GCF_000007125.1_ASM712v1_assembly_report.txt") for link in fetched)
    assert not [
        f
        for f in os.listdir(os.path.join(obj.directory, "source"))
        if f.endswith(".tmp") or f.endswith(".part")
    ]


def test_report_row_with_na_keeps_the_refseq_name_and_the_check_names_it(tmp_path):
    """a sequence without a GenBank accession keeps its RefSeq name, so a scheme that
    needs it is refused by the check with that name"""
    obj, _ = make(tmp_path, GENOME)
    results = Results()
    snps = snp_rows() + [("S.x", "pX", 10, "A", "C")]
    obj.download(
        "GCF_000007125.1",
        "ASM712v1",
        "REF",
        source="refseq",
        directory=obj.directory,
        results=results,
        snps=snps,
    )
    assert results[0][1] is False and "pX" in results[0][2]
    genome = os.path.join(obj.directory, "source", "GCF_000007125.1_ASM712v1_genomic.fna")
    assert "NZ_CP007764.1" in open(genome).read()


def test_failed_report_reports_the_genome_as_failed(tmp_path):
    obj, _ = make(tmp_path, GENOME, fail_report=True)
    results = Results()
    obj.download(
        "GCF_000007125.1",
        "ASM712v1",
        "REF",
        source="refseq",
        directory=obj.directory,
        results=results,
        snps=snp_rows(),
    )
    assert results[0][1] is False and "report not found" in results[0][2]
    assert not os.path.exists(os.path.join(obj.directory, "REF.fna"))
    assert not os.path.exists(os.path.join(obj.directory, "source", "GCF_000007125.1_ASM712v1_genomic.fna"))


def test_genbank_download_does_not_fetch_a_report(tmp_path):
    genbank = fasta(("chrI", CHR1), ("chrII", CHR2))
    obj, fetched = make(tmp_path, genbank)
    results = Results()
    obj.download(
        "GCA_000007125.1",
        "ASM712v1",
        "REF",
        source="genbank",
        directory=obj.directory,
        results=results,
        snps=snp_rows(),
    )
    assert results[0][1] is True
    assert len(fetched) == 1 and not any("assembly_report" in link for link in fetched)


def test_rename_sequences_leaves_no_temporary_file_on_error(tmp_path, monkeypatch):
    path = tmp_path / "g.fna"
    path.write_text(GENOME)

    def boom(src, dst):
        raise OSError("no")

    monkeypatch.setattr(os, "replace", boom)
    try:
        rename_sequences(str(path), {"NZ_CP007762.1": "chrI"})
    except OSError:
        pass
    assert os.listdir(str(tmp_path)) == ["g.fna"] and path.read_text() == GENOME


@pytest.mark.parametrize("umask", [0o022, 0o002, 0o077])
def test_renamed_genome_gets_the_mode_of_open(tmp_path, umask):
    """the renamed genome is written through mkstemp, which gives mode 0600; a reference folder
    is often shared, so it gets the mode open() gives: 0666 masked by the umask"""
    path = tmp_path / "g.fna"
    path.write_text(GENOME)
    previous = os.umask(umask)
    try:
        rename_sequences(str(path), {"NZ_CP007762.1": "chrI"})
    finally:
        os.umask(previous)
    assert os.stat(str(path)).st_mode & 0o777 == 0o666 & ~umask


DIFFERENT_REPORT = (
    "# Sequence-Name\tSequence-Role\tAssigned-Molecule\tAssigned-Molecule-Location/Type\tGenBank-Accn\tRelationship\tRefSeq-Accn\tAssembly-Unit\tSequence-Length\tUCSC-style-name\n"
    "chrI\tassembled-molecule\tI\tChromosome\tchrI\t=\tNZ_CP007762.1\tPrimary Assembly\t600\tna\n"
    "chrII\tassembled-molecule\tII\tChromosome\tchrII\t<>\tNZ_CP007763.1\tPrimary Assembly\t400\tna\n"
    "pX\tassembled-molecule\tpX\tPlasmid\tpX\t<>\tNZ_CP007764.1\tPrimary Assembly\t200\tna\n"
)


def test_genbank_names_skips_rows_whose_sequences_differ(caplog):
    """a RefSeq sequence that differs from its GenBank sequence (Relationship <>) keeps its
    name, as positions on one need not hold on the other; one warning names them"""
    import logging

    with caplog.at_level(logging.WARNING):
        names = genbank_names(DIFFERENT_REPORT)
    assert names == {"NZ_CP007762.1": "chrI"}
    warnings = [r.getMessage() for r in caplog.records if r.levelno == logging.WARNING]
    assert len(warnings) == 1, warnings
    assert "NZ_CP007763.1" in warnings[0] and "NZ_CP007764.1" in warnings[0], warnings


def test_genbank_names_without_a_relationship_column():
    report = "# Sequence-Name\tGenBank-Accn\tRefSeq-Accn\nchrI\tchrI\tNZ_CP007762.1\n"
    assert genbank_names(report) == {"NZ_CP007762.1": "chrI"}
