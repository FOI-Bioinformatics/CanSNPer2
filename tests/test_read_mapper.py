"""Read typing: the read mapper, minimap2 with a read preset, samtools sort and samtools
mpileup at the canSNP positions, run as one command per reference (Aligners.ReadMapper)."""
import pytest

from CanSNPer2.modules.Aligners import ALIGNERS, ReadMapper, threads_per_mapping, write_sites_bed
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
	s = script(ReadMapper(threads=3, max_reads=0).command("ref dir/ref.fna", ("r 1.fq", "r2.fq"), "out/REF_s.pileup"))
	assert "minimap2 -a -x sr -t 3 'ref dir/ref.fna' 'r 1.fq' r2.fq" in s, "paths are quoted"
	assert "| samtools sort -@ 3 -o out/REF_s.pileup.bam -" in s
	assert "samtools index out/REF_s.pileup.bam" in s
	assert "samtools mpileup -l out/REF_s.pileup.sites.bed -a -B -A -d 0 -q 20 -Q 20 -o out/REF_s.pileup out/REF_s.pileup.bam" in s


@pytest.mark.parametrize("platform", ["illumina", "ont", "hifi"])
def test_the_first_300000_reads_of_each_file_are_mapped_by_default(platform):
	'''300000 reads, or 300000 pairs when there are two files: 4 lines per FASTQ record'''
	s = script(ReadMapper(platform=platform).command("ref.fna", ("r 1.fq", "r2.fq"), "o"))
	assert " ref.fna <(head -n 1200000 'r 1.fq') <(head -n 1200000 r2.fq) |" in s


def test_a_compressed_read_file_is_decompressed_before_it_is_cut(tmp_path):
	import gzip
	path = tmp_path / "r.fq.gz"
	with gzip.open(path, "wt") as out:
		out.write("@r1\nACGT\n+\nIIII\n")
	s = script(ReadMapper(max_reads=10).command("ref.fna", (str(path),), "o"))
	assert " ref.fna <(gzip -cd {p} | head -n 40) |".format(p=path) in s


def test_max_reads_zero_maps_every_read():
	s = script(ReadMapper(max_reads=0).command("ref.fna", ("r1.fq", "r2.fq"), "o"))
	assert " ref.fna r1.fq r2.fq |" in s and "head" not in s


def test_threads_of_one_mapping_can_be_given():
	'''the run divides the threads of a sample between the references mapped at the same time'''
	mapper = ReadMapper(threads=8, max_reads=0)
	s = script(mapper.command("ref.fna", ("r.fq",), "o", threads=2))
	assert " -t 2 " in s and "samtools sort -@ 2 " in s
	assert " -t 8 " in script(mapper.command("ref.fna", ("r.fq",), "o")), "all of them by default"


@pytest.mark.parametrize("threads, mappings, each", [(8, 1, 8), (8, 2, 4), (8, 3, 2), (4, 5, 1), (1, 3, 1)])
def test_threads_are_divided_between_mappings(threads, mappings, each):
	assert threads_per_mapping(threads, mappings) == each


def test_quality_thresholds_follow_the_platform_unless_given():
	assert "-q 20 -Q 10 " in script(ReadMapper(platform="ont").command("r.fna", ("x.fq",), "o"))
	assert "-q 30 -Q 7 " in script(ReadMapper(platform="ont", min_mapq=30, min_baseq=7).command("r.fna", ("x.fq",), "o"))


def test_reader_thresholds_follow_the_platform_unless_given():
	r = ReadMapper(platform="ont").reader()
	assert isinstance(r, PileupReader) and (r.min_depth, r.min_allele_fraction) == (5, 0.7)
	r = ReadMapper(platform="ont", min_depth=3, min_allele_fraction=0.9).reader()
	assert (r.min_depth, r.min_allele_fraction) == (3, 0.9)


def test_settings_text_names_everything_that_changes_a_call():
	text = ReadMapper(platform="hifi", min_depth=7).settings()
	assert text == "platform=hifi min_depth=7 min_allele_fraction=0.8 min_mapq=20 min_baseq=20 max_reads=300000"
	assert ReadMapper(max_reads=0).settings().endswith(" max_reads=0")


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


@pytest.mark.parametrize("programs, line", [
	(["minimap2", "samtools"], "read typing: usable"),
	(["minimap2"], "read typing: not usable, samtools not found"),
	(["progressiveMauve"], "read typing: not usable, minimap2 and samtools not found"),
])
def test_selftest_reports_read_typing(tmp_path, programs, line):
	import os, subprocess, sys
	bin_dir = tmp_path / "bin"
	bin_dir.mkdir()
	for program in programs:
		(bin_dir / program).write_text("#!/bin/sh\nexit 0\n")
		(bin_dir / program).chmod(0o755)
	env = dict(os.environ, PATH=str(bin_dir))
	done = subprocess.run([sys.executable, "-m", "CanSNPer2.selftest"], capture_output=True, text=True, env=env)
	assert line in done.stdout, done.stdout
	assert done.returncode == 0, "read typing is optional, an aligner is enough"
