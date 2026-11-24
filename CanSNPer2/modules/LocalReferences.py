'''
Read typing with a local reference: after a first round against the scheme reference(s), the
reads are mapped to the genome of the scheme closest to the sample and the canSNPs are read
at their positions lifted to that genome (LiftOver). This module merges the counts of the two
rounds and turns them into the states the tree walk reads.

Merging, per canSNP: the counts of round 2 are used when they reach min_depth; otherwise those
of round 1. When both rounds reach min_depth and decide the canSNP differently (derived,
ancestral or another base), the canSNP is a conflict and counts as no data: the two genomes
the reads were mapped to give different answers, and neither is preferred.
'''

from CanSNPer2.modules.AlignmentReader import DERIVED
from CanSNPer2.modules.ParsePileup import BASES, STATES

DECIDED = ("derived", "ancestral", "other_base")


def _depth(row):
	return int(row["depth"]) if row else 0


def merge_counts(round1, round2, min_depth):
	'''the site count rows of round 1 and round 2 (dicts by SITE_COUNT_COLUMNS) merged, in the
		order of round 1'''
	second = {row["Name"]: row for row in round2}
	merged = []
	for first in round1:
		other = second.get(first["Name"])
		use2 = other is not None and other["state"] != "not_lifted" and _depth(other) >= min_depth
		if use2 and _depth(first) >= min_depth and first["state"] in DECIDED and other["state"] in DECIDED \
				and first["state"] != other["state"]:
			row = dict(first)
			row["state"] = "conflict"
			row["source"] = "{a},{b}".format(a=first["source"], b=other["source"])
			row["lifted_position"] = other.get("lifted_position", "")
			merged.append(row)
		elif use2:
			merged.append(dict(other))
		else:
			merged.append(dict(first))
	return merged


def states_from_rows(rows):
	'''(states {snp_id: state}, result rows in the columns of SNP_FILE_COLUMNS, derived snp ids)
		of site count rows; sites without reads are left out, as unaligned positions are'''
	states, info, called = {}, [], []
	for row in rows:
		depth = _depth(row)
		if depth == 0:
			continue
		label = row["state"]
		state = STATES[label]
		if label == "derived":
			base = row["Derived base"]
		elif label == "ancestral":
			base = row["Ancestral base"]
		elif label == "other_base":
			base = max(BASES, key=lambda b: int(row[b]))
		else:
			base = "N"
		states[row["Name"]] = state
		info.append([row["Name"], row["Reference"], row["Sequence"], row["Pos"], "reads", str(depth),
				row["Ancestral base"], row["Derived base"], base, ""])
		if state == DERIVED:
			called.append(row["Name"])
	return states, info, called
