"""db export --format vcf (#98): the canSNPs of a database as sites-only VCF files, one per
reference genome, so that other variant callers can genotype the same sites.

Each file is VCFv4.3 with the sequences of the reference genome as contigs. A record is one
canSNP: CHROM and POS are its reference sequence and position, REF the base of the reference
genome there and ALT the other allele of the canSNP (the derived base when the reference
carries the ancestral one, and the other way round). INFO gives the ancestral (AA) and derived
(DA) base, the node (NODE), the rank of the canSNP in its node (RANK) and, for the canSNP the
scheme types, TYPED. The REF base is read from the reference file, so the files of the reference
genomes are needed (cansnper2 ref get); a reference base that is neither allele of its canSNP
is refused, as the typing does (ReferenceCheck)."""

import datetime
import hashlib
import os
import sqlite3

from CanSNPer2 import __version__
from CanSNPer2.modules.ExportDatabase import ExportError, write_outputs
from CanSNPer2.modules.ReferenceCheck import read_sequences

INFO_HEADER = (
    '##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral allele of the canSNP">',
    '##INFO=<ID=DA,Number=1,Type=String,Description="Derived allele of the canSNP">',
    '##INFO=<ID=NODE,Number=1,Type=String,Description="Node of the scheme tree that the canSNP defines">',
    '##INFO=<ID=RANK,Number=1,Type=Integer,Description="Rank of the canSNP within its node, 1 is the first">',
    '##INFO=<ID=TYPED,Number=0,Type=Flag,Description="The canSNP that CanSNPer2 types for its node">',
)

## characters that VCF 4.3 asks to be percent-encoded in INFO values
_ENCODE = {"%": "%25", ":": "%3A", ";": "%3B", "=": "%3D", ",": "%2C", " ": "%20", "\t": "%09"}


def encode(value):
    """an INFO value, percent-encoded as VCF 4.3 requires"""
    return "".join(_ENCODE.get(ch, ch) for ch in str(value))


def file_sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_cansnps(database):
    """{reference genome: [canSNP rows]} and the accession of each genome; a row is
    (node, rank, typed, sequence, position, ancestral, derived)"""
    conn = sqlite3.connect(database)
    try:
        genomes = conn.execute(
            "SELECT id, genome, genbank_id, refseq_id FROM snp_references ORDER BY id"
        ).fetchall()
        rows = conn.execute("""SELECT g.genome, a.snp_id, a.rank, a.typed, a.sequence, a.sequence_position,
				a.ancestral_base, a.derived_base
			FROM snp_annotation a LEFT JOIN snp_references g ON g.id = a.genome_i""").fetchall()
    finally:
        conn.close()
    without_genome = [row[1] for row in rows if row[0] is None]
    if without_genome:
        raise ExportError(
            "{n} canSNP(s) point to a reference genome that is not in the snp_references table: "
            "{snps}".format(n=len(without_genome), snps=", ".join(map(str, without_genome[:10])))
        )
    accessions = {genome: genbank or refseq for _id, genome, genbank, refseq in genomes}
    snps = {}
    for genome, node, rank, typed, sequence, position, ancestral, derived in rows:
        snps.setdefault(genome, []).append((node, rank, typed, sequence, int(position), ancestral, derived))
    return snps, accessions


def vcf_text(genome, accession, records, database, database_sha256, date):
    """the VCF of one reference genome. records: [(sequence name, length), ...] of the
    reference file and the canSNP rows with their REF and ALT, sorted"""
    contigs, rows = records
    header = [
        "##fileformat=VCFv4.3",
        "##fileDate={date}".format(date=date),
        "##source=CanSNPer2 {version}".format(version=__version__),
        "##reference={genome}.fna".format(genome=genome),
        "##cansnper2_reference_genome={genome}{acc}".format(
            genome=genome, acc=" ({0})".format(accession) if accession else ""
        ),
        "##cansnper2_database={db} sha256:{sha}".format(db=os.path.basename(database), sha=database_sha256),
    ]
    header += [
        "##contig=<ID={name},length={length}>".format(name=name, length=length) for name, length in contigs
    ]
    header += list(INFO_HEADER)
    header.append("\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]))
    lines = []
    for sequence, position, ref, alt, node, rank, typed, ancestral, derived in rows:
        info = "AA={aa};DA={da};NODE={node};RANK={rank}".format(
            aa=ancestral, da=derived, node=encode(node), rank=rank
        )
        if typed not in (0, "0"):  ## as the canSNP export reads the typed column (typed_text)
            info += ";TYPED"
        lines.append("\t".join([sequence, str(position), ".", ref, alt, ".", ".", info]))
    return "\n".join(header + lines) + "\n"


def export_vcf(database, outdir, refdir, name=None):
    """Write <name>_<genome>.vcf for each reference genome of database that has canSNPs into
    outdir and return {genome: path}. refdir holds the reference files <genome>.fna. Nothing is
    written when a file is missing or a canSNP does not match its reference (ExportError)."""
    if name is None:
        name = os.path.basename(database).rsplit(".", 1)[0]
    snps, accessions = read_cansnps(database)
    sha = file_sha256(database)
    date = datetime.date.today().strftime("%Y%m%d")
    contents = {}
    paths = {}
    problems = []
    for genome in sorted(snps):
        reference = os.path.join(refdir, "{genome}.fna".format(genome=genome))
        if not os.path.isfile(reference):
            problems.append(
                "the reference genome {genome} has no file {path}; download the reference genomes with "
                "cansnper2 ref get -d {db} -o {refdir}".format(
                    genome=genome, path=reference, db=database, refdir=refdir
                )
            )
            continue
        sequences = read_sequences(reference)
        lookup = dict(sequences)
        order = {seq_name: i for i, (seq_name, _seq) in enumerate(sequences)}
        rows = []
        for node, rank, typed, sequence, position, ancestral, derived in snps[genome]:
            seq = lookup.get(sequence)
            if seq is None or not 1 <= position <= len(seq):
                problems.append(
                    "{node} (rank {rank}) is at {seq}:{pos}, which {genome}.fna does not have".format(
                        node=node, rank=rank, seq=sequence, pos=position, genome=genome
                    )
                )
                continue
            ref = seq[position - 1]
            if ref == ancestral:
                alt = derived
            elif ref == derived:
                alt = ancestral
            else:
                problems.append(
                    "{node} (rank {rank}) at {seq}:{pos}: {genome}.fna has {ref}, neither the ancestral "
                    "{aa} nor the derived {da} base".format(
                        node=node,
                        rank=rank,
                        seq=sequence,
                        pos=position,
                        genome=genome,
                        ref=ref,
                        aa=ancestral,
                        da=derived,
                    )
                )
                continue
            rows.append((sequence, position, ref, alt, node, rank, typed, ancestral, derived))
        rows.sort(key=lambda row: (order[row[0]], row[1], row[4], row[5]))
        contigs = [(seq_name, len(seq)) for seq_name, seq in sequences]
        path = os.path.join(outdir, "{name}_{genome}.vcf".format(name=name, genome=genome))
        contents[path] = vcf_text(genome, accessions.get(genome), (contigs, rows), database, sha, date)
        paths[genome] = path
    if problems:
        raise ExportError(
            "{n} problem(s), nothing was exported:\n  {list}".format(
                n=len(problems), list="\n  ".join(problems[:20])
            )
        )
    if not contents:
        raise ExportError("{db} has no canSNPs to export".format(db=database))
    write_outputs(contents)
    return paths
