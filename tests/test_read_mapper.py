"""Read typing: the read mapper, minimap2 with a read preset, samtools sort and samtools
mpileup at the canSNP positions, run as one command per reference (Aligners.ReadMapper)."""
import pytest

from CanSNPer2.modules.Aligners import ALIGNERS, ReadMapper, write_sites_bed
from CanSNPer2.modules.ParsePileup import PileupReader


def script(cmd):
	assert cmd[:3] == ["bash", "-o", "pipefail"] and cmd[3] == "-c"
	return cmd[4]


def test_the_read_mapper_is_not_an_aligner_choice():
	assert "minimap2-reads" not in ALIGNERS, "it is chosen by the input, not by --aligner"


@pytest.mark.parametrize("platform, preset", [("illumina", "sr"), ("ont", "map-ont"), ("hifi", "map-hifi")])
def test_preset_per_platform(platform, preset):
	s = script(ReadMapper(platform=platform).command("ref.fna", ("r.fq",), "out/REF_s.pileup"))
	assert "minimap2 -a -x {p} ".format(p=preset) in s


def test_paired_reads_and_the_pipeline():
	s = script(ReadMapper(threads=3).command("ref dir/ref.fna", ("r 1.fq", "r2.fq"), "out/REF_s.pileup"))
	assert "minimap2 -a -x sr -t 3 'ref dir/ref.fna' 'r 1.fq' r2.fq" in s, "paths are quoted"
	assert "| samtools sort -@ 3 -o out/REF_s.pileup.bam -" in s
	assert "samtools index out/REF_s.pileup.bam" in s
	assert "samtools mpileup -l out/REF_s.pileup.sites.bed -a -B -A -d 0 -q 20 -Q 20 -o out/REF_s.pileup out/REF_s.pileup.bam" in s


def test_quality_thresholds_follow_the_platform_unless_given():
	assert "-q 20 -Q 10 " in script(ReadMapper(platform="ont").command("r.fna", ("x.fq",), "o"))
	assert "-q 30 -Q 7 " in script(ReadMapper(platform="ont", min_mapq=30, min_baseq=7).command("r.fna", ("x.fq",), "o"))


def test_reader_thresholds_follow_the_platform_unless_given():
	r = ReadMapper(platform="ont").reader()
	assert isinstance(r, PileupReader) and (r.min_depth, r.min_allele_fraction) == (10, 0.7)
	r = ReadMapper(platform="ont", min_depth=3, min_allele_fraction=0.9).reader()
	assert (r.min_depth, r.min_allele_fraction) == (3, 0.9)


def test_settings_text_names_everything_that_changes_a_call():
	text = ReadMapper(platform="hifi", min_depth=7).settings()
	assert text == "platform=hifi min_depth=7 min_allele_fraction=0.8 min_mapq=20 min_baseq=20"


def test_unknown_platform():
	with pytest.raises(ValueError, match="platform"):
		ReadMapper(platform="sanger")


def test_available_needs_both_programs(tmp_path):
	assert not ReadMapper(path=str(tmp_path)).available(), "no minimap2 in that folder"
	assert not ReadMapper(samtools=str(tmp_path / "samtools")).available(), "no samtools there"


def test_sites_bed(tmp_path):
	path = tmp_path / "s.bed"
	write_sites_bed(str(path), [("B.2", "chr", 20, "A", "G"), ("B.1", "chr", 5, "C", "T"), ("B.3", "plasmid", 1, "A", "C")])
	assert path.read_text() == "chr\t4\t5\nchr\t19\t20\nplasmid\t0\t1\n"
