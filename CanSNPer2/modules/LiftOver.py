'''
Lifting canSNP positions from a reference genome of the scheme to another genome of the same
organism, a local reference for read typing (see LocalReferences).

The local genome is aligned to the scheme reference with the minimap2 assembly aligner that
assembly typing uses (Aligners.Minimap2, asm20), and ParsePAF reports for every canSNP the
contig, position and strand of the local genome and the local base, on the strand of the
scheme reference. A canSNP is lifted only when the local genome has its ancestral or derived
base there: a gap, another base, or a position two alignments disagree on is not lifted, and
read typing then uses the scheme reference for it.

The lift is cached in a file per pair of genomes, made again when either genome file, the
canSNPs or the minimap2 version differ from those it was made with.
'''

import hashlib
import os
import subprocess
from collections import namedtuple

from CanSNPer2.modules.Aligners import Minimap2

Lifted = namedtuple("Lifted", ["contig", "position", "strand", "base"])

LIFT_SUFFIX = ".lift.tsv"


def _sha256(path):
	digest = hashlib.sha256()
	with open(path, "rb") as fin:
		for block in iter(lambda: fin.read(1 << 20), b""):
			digest.update(block)
	return digest.hexdigest()


def _name(path):
	return os.path.basename(path).rsplit(".", 1)[0]


def lift_snps(reference_fasta, local_fasta, snps, tmpdir, aligner=None):
	'''{snp_id: Lifted} of the canSNPs (snp_id, sequence, position, ancestral, derived) of the
		reference that can be lifted to the local genome'''
	aligner = aligner or Minimap2()
	paf = os.path.join(tmpdir, "lift_{local}__{ref}.paf".format(local=_name(local_fasta), ref=_name(reference_fasta)))
	run = subprocess.run(aligner.command(reference_fasta, local_fasta, paf), capture_output=True, text=True)
	if run.returncode != 0:
		raise RuntimeError("minimap2 could not align {local} to {ref}: {err}".format(
				local=local_fasta, ref=reference_fasta, err=run.stderr.strip()[-300:]))
	reader = aligner.reader()
	reader._reset()
	reader.reference = _name(reference_fasta)
	reader.snplist = {(seq, int(pos)): (anc, der, snp_id) for snp_id, seq, pos, anc, der in snps}
	positions = {}
	for seq, pos in reader.snplist:
		positions.setdefault(seq, []).append(pos)
	reader.snp_positions = {seq: sorted(p) for seq, p in positions.items()}
	reader.read_alignment(paf)
	os.remove(paf)
	lifted = {}
	for key, obs in reader.observations.items():
		anc, der, snp_id = reader.snplist[key]
		if obs["base"] in (anc, der) and obs["query_pos"] != "N/A":
			lifted[snp_id] = Lifted(obs["query_name"], int(obs["query_pos"]), obs["query_strand"], obs["base"])
	return lifted


def _signature(reference_fasta, local_fasta, snps, aligner):
	version = subprocess.run([aligner.executable(), "--version"], capture_output=True, text=True).stdout.strip()
	snp_text = "\n".join("\t".join(str(v) for v in row) for row in sorted(snps))
	return "reference_sha256={r} local_sha256={l} snps_sha256={s} minimap2={v} options={o}".format(
			r=_sha256(reference_fasta), l=_sha256(local_fasta), s=hashlib.sha256(snp_text.encode()).hexdigest(),
			v=version, o=",".join([aligner.preset] + list(aligner.options)))


def lift(reference_fasta, local_fasta, snps, cache_dir, tmpdir, aligner=None):
	'''lift_snps through the cache file <cache_dir>/<local>__<reference>.lift.tsv'''
	aligner = aligner or Minimap2()
	os.makedirs(cache_dir, exist_ok=True)
	path = os.path.join(cache_dir, "{local}__{ref}{suffix}".format(local=_name(local_fasta), ref=_name(reference_fasta), suffix=LIFT_SUFFIX))
	signature = _signature(reference_fasta, local_fasta, snps, aligner)
	if os.path.isfile(path):
		with open(path) as fin:
			if fin.readline().rstrip("\n") == "# " + signature:
				fin.readline()  ## column header
				lifted = {}
				for line in fin:
					snp_id, contig, position, strand, base = line.rstrip("\n").split("\t")
					lifted[snp_id] = Lifted(contig, int(position), strand, base)
				return lifted
	lifted = lift_snps(reference_fasta, local_fasta, snps, tmpdir, aligner)
	partial = path + ".part"
	with open(partial, "w") as out:
		out.write("# " + signature + "\n")
		out.write("snp_id\tcontig\tposition\tstrand\tbase\n")
		for snp_id in sorted(lifted):
			lf = lifted[snp_id]
			out.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(snp_id, lf.contig, lf.position, lf.strand, lf.base))
	os.replace(partial, path)
	return lifted
