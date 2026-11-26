"""
ParsePAF reads the PAF alignments of minimap2 (run with -c --cs=long; the short cs form is read too) and reports the
query base at every SNP position of the reference. The classification of the bases, the
--mask filter, the edge flag and the result rows are shared with the other alignment
readers, see AlignmentReader.

The database stores each SNP as a reference sequence and a position on it, which are the
coordinates PAF gives for the target (columns 6 and 8-9, the start 0-based). The alignment
is therefore read without the reference file or any other table, also when it is read
again later (--skip_alignment). A target sequence without registered SNPs is passed over.

Which alignments are used: only primary alignments (tp:A:P) with a mapping quality of at
least min_mapq. A secondary alignment or a primary one with a low mapping quality is a
repeat or paralog, whose base can belong to another copy of the region. A reference
position covered by two primary alignments with different query bases is reported as N
(AlignmentReader.observe).
"""

import re
import logging

from CanSNPer2.modules.AlignmentReader import AlignmentReader

logger = logging.getLogger(__name__)

"""the long cs tag: "=ACGT" identical bases, "*ag" reference a / query g, "+acg" bases only
	in the query, "-acg" bases only in the reference, "~gt12ag" an intron (not produced by the
	assembly presets, skipped as reference without query). The short form (--cs=short) writes
	":n" for n identical bases instead of "=ACGT"; those bases are the reference bases, which the
	database holds for the SNP positions, so a SNP inside a ":n" run is read as that base."""
CS_OP = re.compile(r"=([A-Za-z]+)|\*([a-z])([a-z])|\+([a-z]+)|-([a-z]+)|~[a-z]{2}(\d+)[a-z]{2}|:(\d+)")


class ParsePAF(AlignmentReader):
    """Reader of minimap2 PAF files with the long cs tag; the reference is the target"""

    def __init__(self, verbose=False, min_mapq=20, **kwargs):
        super(ParsePAF, self).__init__(verbose=verbose, **kwargs)
        self.min_mapq = int(min_mapq)

    def read_alignment(self, path):
        used = skipped = 0
        with open(path) as fin:
            for line in fin:
                if not line.strip():
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 12:
                    skipped += 1
                    continue
                tags = fields[12:]
                if "tp:A:P" not in tags or int(fields[11]) < self.min_mapq:
                    skipped += 1
                    continue
                cs = next((tag[5:] for tag in tags if tag.startswith("cs:Z:")), None)
                if cs is None:
                    raise ValueError(
                        "{path} has no cs tag, minimap2 must be run with -c --cs=long".format(path=path)
                    )
                used += 1
                self.read_record(fields, cs)
        if skipped and not used:
            logger.warning(
                "{path}: {skipped} record(s) skipped and none used: a record is used only when it has 12 columns, the tag tp:A:P and a mapping quality of at least {q}; check that minimap2 was run with -c --cs=long".format(
                    path=path, skipped=skipped, q=self.min_mapq
                )
            )
        logger.debug(
            "{path}: {used} alignments used, {skipped} secondary or below mapping quality {q} skipped".format(
                path=path, used=used, skipped=skipped, q=self.min_mapq
            )
        )

    def read_record(self, fields, cs):
        """report the query base at the SNP positions one alignment covers. The cs tag is on
        the forward strand of the reference, with the query bases complemented for a "-"
        alignment, so its bases need no further complement"""
        qname, qstart, qend, strand = fields[0], int(fields[2]), int(fields[3]), fields[4]
        sequence = fields[5]
        start = int(fields[7]) + 1  ## 1-based position on the reference sequence
        end = int(fields[8])
        wanted = self.snps_in(sequence, start, end)
        if not wanted:
            return
        block_length = end - start + 1
        next_snp = 0
        at = start  ## reference position of the next reference base in the cs string
        qdone = 0  ## query bases consumed so far

        def query_pos(consumed):
            """position in the query sequence of the query base after "consumed" bases"""
            return str(qstart + consumed + 1) if strand == "+" else str(qend - consumed)

        for match in CS_OP.finditer(cs):
            same, rbase, qbase, insertion, deletion, intron, short = match.groups()
            if insertion:
                qdone += len(insertion)
                continue
            if same:
                span, query_bases = len(same), same
            elif short:
                span, query_bases = int(short), None
            elif rbase:
                span, query_bases = 1, qbase
            elif deletion:
                span, query_bases = len(deletion), None
            else:
                span, query_bases = int(intron), None
            identical = bool(short)
            while next_snp < len(wanted) and wanted[next_snp] < at + span:
                pos = wanted[next_snp]
                k = pos - at
                if identical:
                    base = self.snplist[(sequence, pos)][0]  ## identical bases: the reference base
                    self.observe(
                        sequence,
                        pos,
                        base,
                        pos - start + 1,
                        block_length,
                        query_name=qname,
                        query_pos=query_pos(qdone + k),
                        query_strand=strand,
                    )
                elif query_bases is None:
                    self.observe(
                        sequence,
                        pos,
                        "-",
                        pos - start + 1,
                        block_length,
                        query_name=qname,
                        query_pos="N/A",
                        query_strand=strand,
                    )
                else:
                    self.observe(
                        sequence,
                        pos,
                        query_bases[k],
                        pos - start + 1,
                        block_length,
                        query_name=qname,
                        query_pos=query_pos(qdone + k),
                        query_strand=strand,
                    )
                next_snp += 1
            at += span
            if query_bases is not None or identical:
                qdone += span
            if next_snp == len(wanted):
                break
