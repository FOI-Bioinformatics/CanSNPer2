"""Read typing: telling read files from assemblies, grouping paired-end files into samples
and the samplesheet (ReadInput)."""
import gzip

import pytest

from CanSNPer2.modules.ReadInput import (ReadInputError, ReadSample, input_kind, pair_read_files,
	read_samplesheet, read_sample_name)

FASTQ = "@r1\nACGT\n+\nIIII\n"


def write(path, text, gz=False):
	if gz:
		with gzip.open(path, "wt") as out:
			out.write(text)
	else:
		path.write_text(text)
	return str(path)


def test_input_kind_from_content(tmp_path):
	assert input_kind(write(tmp_path / "a.fastq", FASTQ)) == "reads"
	assert input_kind(write(tmp_path / "a.fq.gz", FASTQ, gz=True)) == "reads"
	assert input_kind(write(tmp_path / "odd_name.txt", FASTQ)) == "reads", "content, not the file name"
	assert input_kind(write(tmp_path / "g.fasta", ">c\nACGT\n")) == "fasta"
	assert input_kind(write(tmp_path / "g.fna.gz", ">c\nACGT\n", gz=True)) == "fasta"
	assert input_kind(write(tmp_path / "list.txt", "a.fasta\nb.fasta\n")) is None
	assert input_kind(write(tmp_path / "empty.fq", "")) is None
	assert input_kind(str(tmp_path / "missing.fq")) is None


@pytest.mark.parametrize("name, sample", [
	("B_iv_01_M_S1_R1_001.fastq.gz", "B_iv_01_M_S1"),
	("B_iv_01_M_S1_R2_001.fastq.gz", "B_iv_01_M_S1"),
	("x_R1.fq.gz", "x"),
	("SRR123_1.fastq", "SRR123"),
	("SRR123_2.fastq.gz", "SRR123"),
	("barcode05.fastq.gz", "barcode05"),
	("sample.fq", "sample"),
])
def test_read_sample_name(name, sample):
	assert read_sample_name(name) == sample


def test_pairs_by_mate_names(tmp_path):
	files = [write(tmp_path / n, FASTQ) for n in ("a_R2_001.fastq", "a_R1_001.fastq", "b_1.fq", "b_2.fq", "c.fq")]
	samples = pair_read_files(files)
	assert samples == [
		ReadSample("a", (files[1], files[0]), "illumina"),
		ReadSample("b", (files[2], files[3]), "illumina"),
		ReadSample("c", (files[4],), "illumina"),
	]


def test_a_mate_without_its_partner_is_refused(tmp_path):
	files = [write(tmp_path / "a_R2_001.fastq", FASTQ)]
	with pytest.raises(ReadInputError, match="a_R2_001.fastq.*mate"):
		pair_read_files(files)


def test_long_reads_are_not_paired(tmp_path):
	files = [write(tmp_path / n, FASTQ) for n in ("run_1.fastq", "run_2.fastq")]
	samples = pair_read_files(files, platform="ont")
	assert [s.files for s in samples] == [(files[0],), (files[1],)]
	assert [s.name for s in samples] == ["run_1", "run_2"], "long read files keep their full names"
	assert {s.platform for s in samples} == {"ont"}


def test_two_files_of_one_sample_name_are_refused(tmp_path):
	files = [write(tmp_path / n, FASTQ) for n in ("a.fq", "a.fastq")]
	with pytest.raises(ReadInputError, match="sample name a"):
		pair_read_files(files)


def test_samplesheet(tmp_path):
	r1 = write(tmp_path / "x1.fq", FASTQ)
	write(tmp_path / "x2.fq", FASTQ)
	write(tmp_path / "long.fq", FASTQ)
	sheet = tmp_path / "samples.tsv"
	sheet.write_text("sample\treads1\treads2\tplatform\nS1\tx1.fq\tx2.fq\t\nS2\tlong.fq\t\tont\n")
	samples = read_samplesheet(str(sheet))
	assert samples[0] == ReadSample("S1", (r1, str(tmp_path / "x2.fq")), "illumina"), "relative paths are read from the folder of the sheet"
	assert samples[1] == ReadSample("S2", (str(tmp_path / "long.fq"),), "ont")


def test_samplesheet_without_platform_column_uses_the_default(tmp_path):
	write(tmp_path / "l.fq", FASTQ)
	sheet = tmp_path / "s.tsv"
	sheet.write_text("sample\treads1\nL\tl.fq\n")
	assert read_samplesheet(str(sheet), platform="hifi")[0].platform == "hifi"


@pytest.mark.parametrize("body, message", [
	("S1\tnone.fq\n", "none.fq.*not found"),
	("S1\tx1.fq\nS1\tx1.fq\n", "S1.*more than once"),
	("S1\t\n", "no read file"),
	("S1\tx1.fq\t\tpacbio\n", "platform pacbio"),
	("", "no samples"),
	("lab/x7\tx1.fq\n", "sample name 'lab/x7' is not allowed"),
	("..\tx1.fq\n", "sample name '..' is not allowed"),
	(".hidden\tx1.fq\n", "sample name '.hidden' is not allowed"),
	("a b\tx1.fq\n", "sample name 'a b' is not allowed"),
])
def test_samplesheet_errors(tmp_path, body, message):
	write(tmp_path / "x1.fq", FASTQ)
	sheet = tmp_path / "s.tsv"
	sheet.write_text("sample\treads1\treads2\tplatform\n" + body)
	with pytest.raises(ReadInputError, match=message):
		read_samplesheet(str(sheet))


def test_samplesheet_needs_its_header(tmp_path):
	sheet = tmp_path / "s.tsv"
	sheet.write_text("S1\tx.fq\n")
	with pytest.raises(ReadInputError, match="header"):
		read_samplesheet(str(sheet))
