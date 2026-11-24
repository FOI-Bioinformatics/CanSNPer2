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

import gzip
import os
import sqlite3
import subprocess
from collections import namedtuple

from CanSNPer2.modules.AlignmentReader import DERIVED
from CanSNPer2.modules.ParsePileup import BASES, STATES
from CanSNPer2.modules.SchemeGenomes import read_genomes

'''separates the genome name from the contig name in the combined fasta of the screen'''
SEP = "|"

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


Candidate = namedtuple("Candidate", ["genome", "node", "path"])

LOCAL_FOLDER = "local"


def local_path(refdir, genome):
	'''the file of a local reference genome: <refdir>/local/<genome>.fna'''
	return os.path.join(refdir, LOCAL_FOLDER, genome + ".fna")


def _subtree(conn, node):
	'''{node name: name of the child of node it descends from (node itself for node)} of the
		subtree of node; the whole tree when node is not in it (no call)'''
	names = dict(conn.execute("SELECT id, name FROM nodes"))
	children = {}
	for parent, child in conn.execute("SELECT parent, child FROM tree WHERE parent != child"):
		children.setdefault(parent, []).append(child)
	ids = {name: i for i, name in names.items()}
	if node in ids:
		start = ids[node]
	else:
		parents = {child for kids in children.values() for child in kids}
		roots = [i for i in names if i not in parents]
		start = roots[0] if roots else None
	branch = {names[start]: names[start]} if start is not None else {}
	for top in children.get(start, []):
		stack = [top]
		while stack:
			i = stack.pop()
			branch[names[i]] = names[top]
			stack.extend(children.get(i, []))
	return branch


def candidates(database, node, refdir, limit=50):
	'''(candidates, missing genomes): the representative genomes of the scheme below node (all
		of them when node is not in the tree, as for no call), as Candidate(genome, node, path),
		at most limit of them, taken in turn from the subtrees of the children of node; and the
		representatives whose file is not in <refdir>/local'''
	conn = sqlite3.connect(str(database))
	try:
		branch = _subtree(conn, node)
		rows = [row for row in read_genomes(conn) if row[2] == "representative" and row[1] in branch]
	finally:
		conn.close()
	found, missing = [], []
	groups = {}
	for genome, at, _role in rows:
		path = local_path(refdir, genome)
		if not os.path.isfile(path):
			missing.append(genome)
			continue
		groups.setdefault(branch[at], []).append(Candidate(genome, at, path))
	queues = [list(group) for _key, group in groups.items()]
	while queues and len(found) < limit:
		for queue in list(queues):
			if len(found) >= limit:
				break
			found.append(queue.pop(0))
			if not queue:
				queues.remove(queue)
	order = {genome: i for i, (genome, _n, _r) in enumerate(rows)}
	found.sort(key=lambda c: order[c.genome])
	return found, missing


def _head_reads(path, n, out):
	'''copy the first n FASTQ records of path (gzip or not) to the open file out'''
	opener = gzip.open if open(path, "rb").read(2) == b"\x1f\x8b" else open
	with opener(path, "rt") as fin:
		for i, line in enumerate(fin):
			if i >= 4 * n:
				break
			out.write(line)


def screen(read_files, genomes, preset, n_reads, tmpdir, threads=2, minimap2="minimap2"):
	'''[(genome, share, matched bases)] of the candidate genomes {name: fasta}, best first: the
		first n_reads reads of each read file are mapped to all genomes at once, and each genome
		is scored by the matched bases of the primary alignments on it; share is its part of
		the total. Mapping qualities are not used: between close candidates they are low by
		design, and the matched bases still favour the closest genome'''
	combined = os.path.join(tmpdir, "screen_genomes.fna")
	with open(combined, "w") as out:
		for name, path in genomes.items():
			with open(path) as fin:
				for line in fin:
					out.write(">" + name + SEP + line[1:].split()[0] + "\n" if line.startswith(">") else line)
	reads = os.path.join(tmpdir, "screen_reads.fq")
	with open(reads, "w") as out:
		for path in read_files:
			_head_reads(path, n_reads, out)
	run = subprocess.run([minimap2, "-x", preset, "-t", str(threads), "--secondary=no", combined, reads],
			capture_output=True, text=True)
	if run.returncode != 0:
		raise RuntimeError("minimap2 failed in the screen of the local references: " + run.stderr.strip()[-300:])
	score = dict.fromkeys(genomes, 0)
	for line in run.stdout.splitlines():
		fields = line.split("\t")
		if len(fields) > 12 and "tp:A:P" in fields[12:]:
			score[fields[5].split(SEP, 1)[0]] += int(fields[9])
	os.remove(combined)
	os.remove(reads)
	total = sum(score.values()) or 1
	return sorted(((g, s / total, s) for g, s in score.items()), key=lambda t: -t[2])


def node_path(database, node):
	'''the names of the nodes from the root down to node, node included; [] for a node not in
		the tree'''
	conn = sqlite3.connect(str(database))
	try:
		names = dict(conn.execute("SELECT id, name FROM nodes"))
		parent = {child: par for par, child in conn.execute("SELECT parent, child FROM tree WHERE parent != child")}
	finally:
		conn.close()
	ids = {name: i for i, name in names.items()}
	if node not in ids:
		return []
	path, i = [], ids[node]
	while i is not None:
		path.append(names[i])
		i = parent.get(i)
	return path[::-1]


def descend(database, call, refdir, site_rows, screen_fn):
	'''Choose the local reference by walking down the tree from call (the root without a call).
		At each node, the child subtrees that hold representative genomes with a file are the
		options; a subtree whose own canSNP is ancestral in site_rows is left out. When the
		derived canSNPs of site_rows lie in exactly one of them, the walk goes there without a
		screen; otherwise one representative of each subtree is screened (screen_fn: {genome:
		fasta} -> [(genome, share, matched bases)], best first) and the walk follows the best.
		It ends at a node without such subtrees, whose representative is chosen. Reads that
		match none of the screened genomes end the walk without a choice when nothing was chosen
		before. Returns (Candidate or None, steps), steps being {"node", "how", "share",
		"runner_up"} for each step down.'''
	conn = sqlite3.connect(str(database))
	try:
		names = dict(conn.execute("SELECT id, name FROM nodes"))
		children = {}
		for parent, child in conn.execute("SELECT parent, child FROM tree WHERE parent != child"):
			children.setdefault(parent, []).append(child)
		reps = {}
		ids = {name: i for i, name in names.items()}
		for genome, node, role in read_genomes(conn):
			path = local_path(refdir, genome)
			if role == "representative" and node in ids and os.path.isfile(path):
				reps.setdefault(ids[node], []).append(Candidate(genome, node, path))
	finally:
		conn.close()
	states = {row["Name"]: row["state"] for row in site_rows}
	below = {}
	def subtree(i):
		'''node ids of the subtree of i, i first, in depth-first order'''
		if i not in below:
			out = [i]
			for c in children.get(i, []):
				out.extend(subtree(c))
			below[i] = out
		return below[i]
	def genomes_in(i):
		return [c for n in subtree(i) for c in reps.get(n, [])]
	if call in ids:
		node = ids[call]
	else:
		has_parent = {c for kids in children.values() for c in kids}
		roots = [i for i in names if i not in has_parent]
		if not roots:
			return None, []
		node = roots[0]
	steps = []
	while True:
		options = [c for c in children.get(node, []) if genomes_in(c) and states.get(names[c]) != "ancestral"]
		if not options:
			break
		supported = [c for c in options if any(states.get(names[n]) == "derived" for n in subtree(c))]
		if len(supported) == 1 or len(options) == 1:
			node = supported[0] if len(supported) == 1 else options[0]
			steps.append({"node": names[node], "how": "counts", "share": None, "runner_up": None})
			continue
		pool = supported or options
		first = {genomes_in(c)[0].genome: c for c in pool}
		ranked = screen_fn({g: genomes_in(c)[0].path for g, c in first.items()})
		if not ranked or ranked[0][2] == 0:
			if not steps:
				return None, steps
			break
		node = first[ranked[0][0]]
		steps.append({"node": names[node], "how": "screen", "share": ranked[0][1],
				"runner_up": ranked[1][1] if len(ranked) > 1 else None})
	found = reps.get(node) or genomes_in(node)
	return (found[0] if found else None), steps
