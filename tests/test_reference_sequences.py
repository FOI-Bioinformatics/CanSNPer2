"""References of several sequences (chromosomes, plasmids).

The database stores each SNP as a reference sequence (the first word of its fasta header)
and a position on that sequence. The reference check therefore looks every SNP up on its
own sequence:

- every sequence the database names must be in the file; the order of the sequences and
  sequences the database does not name (plasmids) do not matter;
- each sequence must be at least as long as its largest SNP position;
- the reference base at each SNP position must be the ancestral or the derived base of
  that SNP. A reference that fails this at many SNPs is refused, a few mismatches are
  reported as a warning.

A reference whose sequence names match none of the database (RefSeq names in a file for a
scheme built on GenBank accessions) used to give NA calls without a message, also when the
alignments were read again with --skip_alignment."""
import gzip
import logging
import random

import pytest

from helpers import make_db
from CanSNPer2.DownloadGenomes import DownloadGenomes
from CanSNPer2.modules.CanSNPer2 import CanSNPer2, CanSNPer2Error
from CanSNPer2.modules.DatabaseConnection import CanSNPdbFunctions
from CanSNPer2.modules.ParsePAF import ParsePAF
from CanSNPer2.modules.ReferenceCheck import check_reference, read_sequences, reference_snps

rng = random.Random(7)
CHR1 = "".join(rng.choice("ACGT") for _ in range(600))
CHR2 = "".join(rng.choice("ACGT") for _ in range(400))
PLASMID = "".join(rng.choice("ACGT") for _ in range(200))
SEQUENCES = {"chrI": CHR1, "chrII": CHR2}
POSITIONS = [("chrI", 50), ("chrI", 120), ("chrI", 333), ("chrI", 480), ("chrI", 599),
	("chrII", 1), ("chrII", 100), ("chrII", 250), ("chrII", 399), ("chrII", 400)]
OTHER = {"A": "G", "C": "T", "G": "A", "T": "C"}


def fasta(*records):
	return "".join(">{0} test\n{1}\n".format(name, seq) for name, seq in records)


def snp_rows():
	'''(snp_id, sequence, position, ancestral, derived): the reference carries the ancestral
		base, the derived base is another one'''
	rows = []
	for i, (sequence, pos) in enumerate(POSITIONS):
		base = SEQUENCES[sequence][pos - 1]
		rows.append(("S.{0}".format(i), sequence, pos, base, OTHER[base]))
	return rows


def write(path, text):
	path.write_text(text)
	return str(path)


def neither(sequence, positions):
	'''sequence with a base that is neither the ancestral nor the derived one at positions'''
	seq = list(sequence)
	for pos in positions:
		base = seq[pos - 1]
		seq[pos - 1] = next(b for b in "ACGT" if b not in (base, OTHER[base]))
	return "".join(seq)


def paf_line(sequence, position, query_base, qname="ctg1"):
	'''one primary minimap2 record of 11 reference bases around position, with query_base at it'''
	ref = SEQUENCES.get(sequence, CHR2)
	start = position - 6        ## 0-based start, 5 bases before the SNP
	left, base, right = ref[start:position - 1], ref[position - 1], ref[position:position + 5]
	if query_base == base:
		cs = "=" + left + base + right
	else:
		cs = "={0}*{1}{2}={3}".format(left, base.lower(), query_base.lower(), right)
	fields = [qname, "1000", 100, 111, "+", sequence, "1000", start, start + 11, 11, 11, 60, "tp:A:P", "cs:Z:" + cs]
	return "\t".join(str(f) for f in fields)


def test_reference_snps_are_read_per_sequence(tmp_path):
	db = make_db(tmp_path / "t.db", snp_rows())
	snps = reference_snps(CanSNPdbFunctions(db))
	assert sorted(snps) == ["REF"]
	assert sorted(snps["REF"]) == sorted(snp_rows())


def test_read_sequences_gives_names_and_upper_case_sequences(tmp_path):
	path = write(tmp_path / "r.fna", ">chrI first\nacgt\nAC\n>chrII\nGG\n")
	assert read_sequences(path) == [("chrI", "ACGTAC"), ("chrII", "GG")]


def test_the_right_reference_passes(tmp_path):
	right = write(tmp_path / "right.fna", fasta(("chrI", CHR1), ("chrII", CHR2)))
	assert check_reference(right, "REF", snp_rows()) is None


def type_with_skipped_alignment(tmp_path, reference_text, name):
	'''type a PAF alignment again (--skip_alignment) against a reference file with the given
		text; B.1 (chrII:100) is derived in the query, S.0 (chrI:50) ancestral'''
	work = tmp_path / name
	(work / "refs").mkdir(parents=True)
	(work / "shared").mkdir()
	(work / "refs" / "REF.fna").write_text(reference_text)
	query = work / "sample.fna"
	query.write_text(fasta(("ctg1", CHR1)))
	db = make_db(work / "t.db", snp_rows())
	anc = CHR2[99]
	(work / "shared" / "REF_sample.paf").write_text(paf_line("chrII", 100, OTHER[anc]) + "\n" + paf_line("chrI", 50, CHR1[49], qname="ctg2") + "\n")
	obj = CanSNPer2([str(query)], workdir=str(work), refdir="refs", tmpdir="shared", logdir="logs", outdir="out",
		database=db, aligner="minimap2", skip_alignment=True)
	files = obj.align(str(query))
	states, _info, called = obj.find_snps_multiproc(xmfa_obj=obj.aligner.reader(database=db), xmfa_files=files,
		references=obj.alignment_references)
	return states, called


def test_reference_in_another_order_passes_and_types_the_same(tmp_path):
	right = fasta(("chrI", CHR1), ("chrII", CHR2))
	swapped = fasta(("chrII", CHR2), ("chrI", CHR1))
	assert check_reference(write(tmp_path / "swapped.fna", swapped), "REF", snp_rows()) is None
	assert type_with_skipped_alignment(tmp_path, right, "right") == type_with_skipped_alignment(tmp_path, swapped, "swapped")
	states, called = type_with_skipped_alignment(tmp_path, swapped, "again")
	assert states == {"S.6": 1, "S.0": 2} and called == ["S.6"]


def test_reference_with_an_extra_plasmid_passes(tmp_path):
	path = write(tmp_path / "plasmid.fna", fasta(("pXO1", PLASMID), ("chrI", CHR1), ("chrII", CHR2)))
	assert check_reference(path, "REF", snp_rows()) is None


def test_header_with_description_matches_by_first_word(tmp_path):
	path = write(tmp_path / "desc.fna", ">chrI Brucella chromosome I, complete sequence\n{0}\n>chrII Brucella chromosome II\n{1}\n".format(CHR1, CHR2))
	assert check_reference(path, "REF", snp_rows()) is None


def test_missing_sequence_is_refused_naming_it(tmp_path):
	path = write(tmp_path / "missing.fna", fasta(("chrI", CHR1), ("NZ_CP007763.1", CHR2)))
	problem = check_reference(path, "REF", snp_rows())
	assert problem is not None
	assert "chrII" in problem and "NZ_CP007763.1" in problem
	assert "GenBank" in problem and "cansnper2 ref get" in problem


def test_short_sequence_is_refused(tmp_path):
	path = write(tmp_path / "short.fna", fasta(("chrI", CHR1), ("chrII", CHR2[:300])))
	problem = check_reference(path, "REF", snp_rows())
	assert problem is not None
	assert "chrII" in problem and "300 bases" in problem and "400" in problem


def test_wrong_bases_on_one_sequence_refused_at_fraction(tmp_path):
	'''5 of 10 SNPs, all on chrII, have neither base: more than MISMATCH_FRACTION of the
		SNPs of the reference'''
	wrong = neither(CHR2, [pos for seq, pos in POSITIONS if seq == "chrII"])
	path = write(tmp_path / "wrong.fna", fasta(("chrI", CHR1), ("chrII", wrong)))
	problem = check_reference(path, "REF", snp_rows())
	assert problem is not None
	assert "5 of 10 SNPs" in problem and "reference base" in problem and "chrII" in problem


def test_few_wrong_bases_warn(tmp_path, caplog):
	'''one wrong base in a scheme is a scheme error to report, not a reason to stop typing'''
	path = write(tmp_path / "few.fna", fasta(("chrI", CHR1), ("chrII", neither(CHR2, [250]))))
	with caplog.at_level(logging.WARNING):
		assert check_reference(path, "REF", snp_rows()) is None
	assert "S.7" in caplog.text and "1 of 10 SNPs" in caplog.text


class Results(list):
	def put(self, item):
		self.append(item)


def test_download_refuses_a_genome_without_a_named_sequence(tmp_path):
	obj = DownloadGenomes.__new__(DownloadGenomes)
	obj.source, obj.directory, obj.database = "genbank", str(tmp_path / "refs"), None
	data = gzip.compress(fasta(("NZ_CP007762.1", CHR1), ("NZ_CP007763.1", CHR2)).encode())
	def retrieve(link, target, timeout=None):
		with open(target, "wb") as out:
			out.write(data)
	obj.retrieve = retrieve
	results = Results()
	obj.download("GCA_000007125.1", "ASM712v1", "REF", directory=obj.directory, results=results, snps=snp_rows())
	assert results[0][1] is False and "chrI" in results[0][2] and "chrII" in results[0][2]


def test_skipped_alignment_checks_the_reference_when_it_is_there(tmp_path):
	'''--skip_alignment used to skip the reference check, so a renamed reference gave NA
		calls without a message'''
	renamed = fasta(("NZ_CP007762.1", CHR1), ("NZ_CP007763.1", CHR2))
	with pytest.raises(CanSNPer2Error, match="chrII"):
		type_with_skipped_alignment(tmp_path, renamed, "renamed")


def test_an_alignment_on_unknown_sequences_is_reported(tmp_path, caplog):
	'''no SNP observed, and the alignment names reference sequences the database does not use'''
	db = make_db(tmp_path / "t.db", snp_rows())
	paf = tmp_path / "REF_sample.paf"
	paf.write_text(paf_line("NZ_CP007763.1", 100, "A") + "\n")
	with caplog.at_level(logging.WARNING):
		snps = ParsePAF().run(str(paf), reference="REF", database=db)
	assert snps == {}
	assert "NZ_CP007763.1" in caplog.text and "chrI" in caplog.text and "chrII" in caplog.text
