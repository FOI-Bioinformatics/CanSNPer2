"""
Which sequencing platform the reads of an input file look like (#119).

The platform sets the minimap2 preset and the thresholds of the pileup (ParsePileup.PLATFORMS),
and --platform was not compared with the reads: ONT reads typed with the default illumina
settings are mapped with the short read preset and filtered at base quality 20, which most of
their bases do not reach, and the run then reports few sites with data or none without saying
why. The platform is therefore read from the input as well: always, to warn when it disagrees
with the platform in effect, and as the platform itself with --platform auto.

The signals are the ones the first few hundred records of a file carry:

        median read length         short reads (illumina) from long reads
        mean base quality          clr (no usable qualities) from ont, and ont from hifi
        read name form             the instrument that wrote the file, when the names are still there
        @RG PL: of a BAM or CRAM   the platform the instrument recorded

The read names are checked first, because they name the instrument: a nanopore run of short
fragments is still ONT, so a long read platform read from the names stands whatever the read
lengths are. Only a name form that says illumina is dropped on long reads, where something has
rewritten the names. When the names are gone, as SRA downloads leave them, the lengths and the
qualities decide, and they separate the platforms only roughly: ONT reads now reach the base
quality of HiFi reads, so a file of long reads of high quality is reported as hifi, and a
nanopore run of short fragments as illumina. The platform in effect is not changed by any of
this unless --platform auto was given.
"""

import gzip
import logging
import re
import statistics
import subprocess
from collections import namedtuple

logger = logging.getLogger(__name__)

"""the records read from the head of a file, and an upper bound on the bytes they may take"""
RECORDS = 400
MAX_BYTES = 8 * 1024 * 1024

"""a median read length below this is short reads; the mean base quality from which long reads
	are taken for hifi, and the one below which the qualities carry no information (clr)"""
LONG_READ = 500
HIFI_BASEQ = 25
NO_BASEQ = 3

"""the platform of a read file: platform is None when the file does not say, source is
	read_group, read_names or reads, and note holds what else is worth reporting"""
Detection = namedtuple("Detection", ["platform", "source", "reads", "median_length", "mean_baseq", "note"])

UNKNOWN = Detection(None, "none", 0, None, None, "")

"""the read name forms of the instruments, in the order they are tried: PacBio
	(movie/zmw/ccs for HiFi, movie/zmw/start_end for CLR), ONT (a UUID, or the tags its
	basecallers write) and Illumina (instrument:run:flowcell:lane:tile:x:y)"""
NAME_FORMS = (
    ("hifi", re.compile(r"^m[0-9A-Za-z_]+/\d+/ccs")),
    ("clr", re.compile(r"^m[0-9A-Za-z_]+/\d+/\d+_\d+")),
    ("ont", re.compile(r"^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}")),
    ("ont", re.compile(r"(?:^| )(?:runid|sampleid|flow_cell_id|basecall_model_version_id)=")),
    ("illumina", re.compile(r"^[^:\s]+:\d+:[^:\s]+:\d+:\d+:\d+:\d+")),
)

"""the @RG PL: values of the platforms that have a preset; PACBIO needs the records or the
	program line to tell HiFi from CLR"""
READ_GROUP_PLATFORMS = {
    "ILLUMINA": "illumina",
    "DNBSEQ": "illumina",
    "ELEMENT": "illumina",
    "ONT": "ont",
    "OXFORD_NANOPORE": "ont",
    "NANOPORE": "ont",
}


def _open(path):
    """the file as text, decompressed when it is gzip-compressed"""
    with open(path, "rb") as raw:
        compressed = raw.read(2) == b"\x1f\x8b"
    opener = gzip.open if compressed else open
    return opener(path, "rt", errors="replace")


def _mean_baseq(quality):
    """the mean Phred+33 base quality of a quality line, None when there is none"""
    quality = quality.strip()
    if not quality or quality == "*":
        return None
    return sum(ord(c) - 33 for c in quality) / len(quality)


def _fastq_records(path):
    """the names, the read lengths and the mean base qualities of the first records of a FASTQ
    file; ([], [], []) when the file cannot be read or holds no reads"""
    names, lengths, quals = [], [], []
    read = 0
    try:
        with _open(path) as handle:
            while len(lengths) < RECORDS and read < MAX_BYTES:
                head = handle.readline()
                sequence = handle.readline()
                handle.readline()
                quality = handle.readline()
                if not quality:
                    break
                read += len(head) + len(sequence) + len(quality)
                if not head.startswith("@"):
                    """not a FASTQ file after all"""
                    return [], [], []
                sequence = sequence.strip()
                if not sequence:
                    continue
                names.append(head[1:].strip())
                lengths.append(len(sequence))
                mean = _mean_baseq(quality)
                if mean is not None:
                    quals.append(mean)
    except (OSError, EOFError, UnicodeError) as e:
        logger.debug("the platform of {path} was not read: {e}".format(path=path, e=e))
        return [], [], []
    return names, lengths, quals


def _by_name(names):
    """the platform of the read names, when at least half of them are of one instrument's form"""
    if not names:
        return None
    for platform, form in NAME_FORMS:
        if 2 * sum(1 for name in names if form.search(name)) >= len(names):
            return platform
    return None


def _detection(names, lengths, quals, note=""):
    """the platform of the records of a file: from the read names when they name an instrument,
    else from the read lengths and the base qualities. A long read platform is taken from the
    names whatever the lengths are, since a run of short fragments is still that platform; a
    name form that says illumina is not trusted on long reads, where something has rewritten
    the names"""
    if not lengths:
        return UNKNOWN._replace(note=note)
    median = statistics.median(lengths)
    mean_baseq = statistics.fmean(quals) if quals else None
    long_reads = median >= LONG_READ
    if mean_baseq is None or mean_baseq < NO_BASEQ:
        note = note or "no base qualities"
    found = Detection(None, "reads", len(lengths), median, mean_baseq, note)
    by_name = _by_name(names)
    if by_name is not None and not (by_name == "illumina" and long_reads):
        return found._replace(platform=by_name, source="read_names")
    if not long_reads:
        platform = "illumina"
    elif mean_baseq is None or mean_baseq < NO_BASEQ:
        platform = "clr"
    elif mean_baseq >= HIFI_BASEQ:
        platform = "hifi"
    else:
        platform = "ont"
    return found._replace(platform=platform)


def _header_lines(path, samtools):
    """the header lines of a BAM or CRAM file, [] when samtools cannot read it"""
    try:
        done = subprocess.run([samtools, "view", "-H", path], capture_output=True, text=True, timeout=120)
    except (OSError, subprocess.SubprocessError) as e:
        logger.debug("the header of {path} was not read: {e}".format(path=path, e=e))
        return []
    if done.returncode != 0:
        logger.debug("the header of {path} was not read: {e}".format(path=path, e=done.stderr.strip()))
        return []
    return done.stdout.splitlines()


def _alignment_records(path, samtools):
    """the names, the read lengths and the mean base qualities of the first records of a BAM or
    CRAM file"""
    names, lengths, quals = [], [], []
    try:
        process = subprocess.Popen(
            [samtools, "view", path], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True
        )
    except (OSError, subprocess.SubprocessError) as e:
        logger.debug("the reads of {path} were not read: {e}".format(path=path, e=e))
        return [], [], []
    try:
        for line in process.stdout:
            fields = line.split("\t")
            if len(fields) < 11 or fields[9] == "*":
                continue
            names.append(fields[0])
            lengths.append(len(fields[9]))
            mean = _mean_baseq(fields[10])
            if mean is not None:
                quals.append(mean)
            if len(lengths) >= RECORDS:
                break
    finally:
        process.stdout.close()
        process.terminate()
        process.wait()
    return names, lengths, quals


def _alignment_detection(path, samtools):
    """the platform of a BAM or CRAM file: the @RG PL: tag when it names one, else its records.
    A PacBio file is HiFi when a program line names ccs, and is told from the records otherwise"""
    header = _header_lines(path, samtools)
    tags = {
        field[3:].upper()
        for line in header
        if line.startswith("@RG")
        for field in line.rstrip("\n").split("\t")
        if field.startswith("PL:")
    }
    note = ""
    if len(tags) == 1:
        tag = tags.pop()
        note = "the read group says PL:{tag}".format(tag=tag)
        if tag in READ_GROUP_PLATFORMS:
            return Detection(READ_GROUP_PLATFORMS[tag], "read_group", 0, None, None, note)
        if tag == "PACBIO" and any(line.startswith("@PG") and "ccs" in line.lower() for line in header):
            return Detection("hifi", "read_group", 0, None, None, note)
    elif len(tags) > 1:
        note = "the read groups name more than one platform"
        return UNKNOWN._replace(note=note)
    return _detection(*_alignment_records(path, samtools), note=note)


def detect_platform(path, kind="reads", samtools="samtools"):
    """the platform the reads of a file look like: a Detection whose platform is None when the
    file does not say. kind is the input_kind of the file, "reads" or "alignment" """
    if kind == "alignment":
        return _alignment_detection(path, samtools or "samtools")
    if kind != "reads":
        return UNKNOWN._replace(note="not a read file")
    return _detection(*_fastq_records(path))


def describe(detection):
    """the evidence of a detection in words, for a warning"""
    if detection.source == "read_group":
        return detection.note or "the read group of the file"
    words = []
    if detection.median_length is not None:
        words.append("median read length {n:g} bp".format(n=detection.median_length))
    if detection.mean_baseq is None:
        words.append("no base qualities")
    elif detection.mean_baseq < NO_BASEQ:
        words.append("no base qualities, mean {q:.1f}".format(q=detection.mean_baseq))
    else:
        words.append("mean base quality {q:.1f}".format(q=detection.mean_baseq))
    if detection.source == "read_names":
        words.append("the read names")
    return ", ".join(words)


def reads_line_fields(detection):
    """the detection as the key=value words of the "Reads:" line of a result file; empty when
    there is no detection to report"""
    if detection is None or detection.platform is None:
        return ""
    text = " detected={p} detected_from={s}".format(p=detection.platform, s=detection.source)
    if detection.median_length is not None:
        text += " detected_median_length={n:g}".format(n=detection.median_length)
    if detection.mean_baseq is not None:
        text += " detected_mean_baseq={q:.1f}".format(q=detection.mean_baseq)
    return text
