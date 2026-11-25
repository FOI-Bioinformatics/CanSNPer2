'''
Sequencing reads as query input: telling read files from assemblies, and grouping the read
files into samples.

A read file is told from an assembly by its content (the first record starts with "@"),
gzip-compressed or not. A sample is one read file, or two for paired-end reads. Paired-end
files are found by their mate names (_R1/_R2, with or without the Illumina _001, or _1/_2);
the sample is named after the file without the mate part and the file endings. A file that
looks like one mate of a pair without the other is refused: typing half of a pair as
single-end reads would go unnoticed. Long reads (ont, hifi) are never paired.

A samplesheet gives the samples explicitly, a tab separated file with the header
sample, reads1 and optionally reads2 and platform; relative paths are read from the folder
of the sheet.
'''

import gzip
import os
import re
from collections import namedtuple

from CanSNPer2.modules.ParsePileup import PLATFORMS, DEFAULT_PLATFORM


class ReadInputError(ValueError):
	'''read input that cannot be typed, with the reason'''


ReadSample = namedtuple("ReadSample", ["name", "files", "platform"])

READ_ENDINGS = re.compile(r"(\.(fastq|fq))?(\.gz)?$", re.IGNORECASE)
MATE = re.compile(r"^(?P<stem>.+?)_(?:R(?P<r>[12])(?:_\d{3})?|(?P<n>[12]))$")


def _head(path, size=4096):
	'''the first bytes of a file, decompressed when it is gzip-compressed; b"" when the file
		cannot be read'''
	try:
		with open(path, "rb") as raw:
			magic = raw.read(2)
		opener = gzip.open if magic == b"\x1f\x8b" else open
		with opener(path, "rb") as handle:
			return handle.read(size)
	except (OSError, EOFError):
		return b""


def input_kind(path):
	'''"reads" for a FASTQ file, "fasta" for a FASTA file, None for anything else (a query
		list, an empty or unreadable file)'''
	head = _head(path)
	if head.startswith(b"\xef\xbb\xbf"):
		head = head[3:]
	head = head.lstrip()
	if head.startswith(b"@"):
		return "reads"
	if head.startswith(b">"):
		return "fasta"
	return None


def _stem(path):
	return READ_ENDINGS.sub("", os.path.basename(path))


def _mate(path):
	'''(sample name, 1 or 2) for a file named as one mate of a pair, else None'''
	match = MATE.match(_stem(path))
	if not match:
		return None
	return match.group("stem"), int(match.group("r") or match.group("n"))


def read_sample_name(path):
	'''the sample name of a read file: the file name without its endings and its mate part'''
	mate = _mate(path)
	return mate[0] if mate else _stem(path)


def _check_platform(platform, where):
	if platform not in PLATFORMS:
		raise ReadInputError("{where}: platform {p} is not known, choose one of {known}".format(
				where=where, p=platform, known=", ".join(sorted(PLATFORMS))))


def pair_read_files(paths, platform=DEFAULT_PLATFORM):
	'''the samples of the read files, in the order of their first file'''
	_check_platform(platform, "--platform")
	groups = {}
	order = []
	for path in paths:
		mate = _mate(path) if platform == "illumina" else None
		name = mate[0] if mate else _stem(path)
		if name not in groups:
			groups[name] = {}
			order.append(name)
		slot = mate[1] if mate else 0
		group = groups[name]
		if slot in group or (slot == 0 and group) or (slot != 0 and 0 in group):
			'''two files that are not the two mates of one pair'''
			other = list(group.values())[0]
			raise ReadInputError("{a} and {b} give the same sample name {name}; rename one of them or give a --samplesheet".format(
					a=other, b=path, name=name))
		groups[name][slot] = path
	samples = []
	for name in order:
		files = groups[name]
		if 0 in files:
			samples.append(ReadSample(name, (files[0],), platform))
			continue
		if set(files) != {1, 2}:
			only = list(files.values())[0]
			raise ReadInputError("{f} is named as one mate of paired-end reads but its other mate was not given; give both files, "
				"or a --samplesheet to type it as single-end reads".format(f=only))
		samples.append(ReadSample(name, (files[1], files[2]), platform))
	return samples


def read_samplesheet(path, platform=DEFAULT_PLATFORM):
	'''the samples of a samplesheet (see the module description)'''
	folder = os.path.dirname(os.path.abspath(path))
	with open(path, encoding="utf-8-sig") as fin:
		lines = [line.rstrip("\n\r") for line in fin if line.strip()]
	header = lines[0].split("\t") if lines else []
	if header[:2] != ["sample", "reads1"]:
		raise ReadInputError("samplesheet {path} must start with the header line sample<TAB>reads1, optionally followed by reads2 and platform".format(path=path))
	samples = []
	seen = set()
	for number, line in enumerate(lines[1:], 2):
		row = dict(zip(header, [field.strip() for field in line.split("\t")]))
		name = row.get("sample", "")
		where = "samplesheet {path} line {n}".format(path=path, n=number)
		if not name:
			raise ReadInputError("{where} has no sample name".format(where=where))
		if name in seen:
			raise ReadInputError("{where}: sample {name} is listed more than once".format(where=where, name=name))
		seen.add(name)
		files = [row.get(column, "") for column in ("reads1", "reads2")]
		files = [f if os.path.isabs(f) else os.path.join(folder, f) for f in files if f]
		if not files:
			raise ReadInputError("{where}: sample {name} has no read file".format(where=where, name=name))
		for f in files:
			if not os.path.isfile(f):
				raise ReadInputError("{where}: read file {f} was not found".format(where=where, f=f))
		sample_platform = row.get("platform") or platform
		_check_platform(sample_platform, where)
		samples.append(ReadSample(name, tuple(files), sample_platform))
	if not samples:
		raise ReadInputError("samplesheet {path} lists no samples".format(path=path))
	return samples
