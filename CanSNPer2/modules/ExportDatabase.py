"""
Export a CanSNPer2 database to the source files that CanSNPer2-database --create reads:

        <name>_tree.txt        one line per node, the tab separated path from the root to the node
        <name>_snp.txt         the SNP annotation, one row per SNP
        <name>_references.txt  the reference genomes

These are the files kept in CanSNPer2-data/source. Creating a database from an export gives
the same tree (node names and links), SNPs and references. The children of a node are written in
the order of their node ids. The node ids themselves are not kept: a database built from the
files numbers its nodes in the order of the tree file, which is the same as the original only
when the original numbered its nodes depth first.

The files are never overwritten. They are written to temporary files in the target folder and
renamed into place when all of them are written, so a failed export leaves no partial files.
"""

import logging
import os
import sqlite3
import tempfile

from CanSNPer2.modules.DatabaseConnection import GENOME_SOURCE_COLUMNS
from CanSNPer2.modules.SchemeGenomes import GENOMES_COLUMNS, read_genomes

logger = logging.getLogger(__name__)

## the canSNP header of CanSNPer2-data, so an export can be added there as it is. Every canSNP
## of every node is written, ranked within its node, with the ones the scheme types marked. The
## date column is not written, it is not part of that header; a database built from the export
## gets the build date
CANSNP_COLUMNS = (
    "node",
    "rank",
    "reference",
    "genome",
    "sequence",
    "sequence_position",
    "derived_base",
    "ancestral_base",
    "score",
    "typed",
)
## the typed column as the canSNP file spells it
TYPED_TEXT = {1: "yes", 0: "no"}


class ExportError(Exception):
    """the database cannot be written as source files"""


def text(value):
    return "" if value is None else str(value)


def tree_paths(conn):
    """Return the path from the root to every node, in the order the tree file lists them.
    The source format writes each node with its full path, so every node needs exactly one
    parent and must be reachable from the root, the node that is linked to itself."""
    names = dict(conn.execute("SELECT id,name FROM nodes"))
    parent_of, children, roots = {}, {}, []
    for child, parent in conn.execute("SELECT child,parent FROM tree ORDER BY child"):
        if child == parent:
            roots.append(child)
            continue
        if child in parent_of:
            raise ExportError(
                "node {node} has more than one parent ({a}, {b}), the tree file can only "
                "describe a tree".format(
                    node=names.get(child, child),
                    a=names.get(parent_of[child], parent_of[child]),
                    b=names.get(parent, parent),
                )
            )
        parent_of[child] = parent
        children.setdefault(parent, []).append(child)
    if len(roots) != 1:
        raise ExportError(
            "the database has {n} root nodes (a node linked to itself), the tree file needs "
            "exactly one".format(n=len(roots))
        )
    missing = sorted((set(parent_of) | set(children) | set(roots)) - set(names))
    if missing:
        raise ExportError(
            "the tree links {n} node id(s) that are not in the nodes table: {ids}".format(
                n=len(missing), ids=", ".join(map(str, missing[:10]))
            )
        )
    paths = []
    stack = [(roots[0], [names[roots[0]]])]
    seen = set()
    while stack:
        node, path = stack.pop()
        if node in seen:
            raise ExportError("the tree has a cycle through node {node}".format(node=names[node]))
        seen.add(node)
        paths.append(path)
        for child in reversed(children.get(node, [])):
            stack.append((child, path + [names[child]]))
    unreachable = sorted(names[n] for n in set(parent_of) - seen)
    if unreachable:
        raise ExportError(
            "{n} node(s) are not connected to the root {root}: {names}".format(
                n=len(unreachable), root=names[roots[0]], names=", ".join(unreachable[:10])
            )
        )
    return paths


def write_outputs(contents):
    """Write {path: text} so that either all files exist afterwards or none does. An existing
    file is never replaced. The text goes to temporary files in the folder of each file and
    is renamed into place when all are written. Raises ExportError."""
    for path in contents:
        if os.path.lexists(path):
            raise ExportError(
                "{path} already exists, it is not overwritten; choose another folder with "
                "-o or remove the file".format(path=path)
            )
    temporary = {}
    placed = []
    try:
        for path, content in contents.items():
            folder = os.path.dirname(path) or "."
            os.makedirs(folder, exist_ok=True)
            fd, partial = tempfile.mkstemp(
                prefix=".{name}.".format(name=os.path.basename(path)), suffix=".part", dir=folder
            )
            temporary[path] = partial
            with os.fdopen(fd, "w") as out:
                out.write(content)
            os.chmod(partial, 0o666 & ~_umask())
        for path, partial in temporary.items():
            os.replace(partial, path)
            placed.append(path)
    except OSError as e:
        for path in placed:
            try:
                os.remove(path)
            except OSError:
                pass
        raise ExportError(
            "the files could not be written to {folder}: {error}; nothing was exported".format(
                folder=os.path.dirname(next(iter(contents))) or ".", error=e.strerror or e
            )
        )
    finally:
        for partial in temporary.values():
            if os.path.exists(partial):
                try:
                    os.remove(partial)
                except OSError:
                    pass


def _umask():
    mask = os.umask(0)
    os.umask(mask)
    return mask


def typed_text(value):
    """the typed column of one canSNP row; anything that is not 0 is a canSNP the scheme types,
    so that a database written by hand does not export a canSNP no node would call"""
    return TYPED_TEXT.get(value, "no" if value in (0, "0") else "yes")


def export_source(database, outdir, name=None):
    """Write the tree, canSNP and reference source files of database into outdir and return
    their paths. name is the file prefix, by default the database file name."""
    if name is None:
        name = os.path.basename(database).rsplit(".", 1)[0]
    conn = sqlite3.connect(database)
    try:
        paths = tree_paths(conn)
        references = conn.execute(
            "SELECT {cols} FROM snp_references ORDER BY id".format(cols=",".join(GENOME_SOURCE_COLUMNS))
        ).fetchall()
        ## every canSNP, the alternatives of a node included, so that the export and a database
        ## built from it hold the same canSNPs
        snps = conn.execute("""SELECT a.snp_id, a.rank, a.reference, g.genome, a.sequence, a.sequence_position,
				a.derived_base, a.ancestral_base, a.score, a.typed
			FROM snp_annotation a LEFT JOIN snp_references g ON g.id = a.genome_i
			ORDER BY a.node_id, a.rank""").fetchall()
    finally:
        conn.close()
    scheme_genomes = read_genomes(database)
    without_genome = [row[0] for row in snps if row[3] is None]
    if without_genome:
        raise ExportError(
            "{n} canSNP(s) point to a reference genome that is not in the snp_references table: "
            "{snps}".format(n=len(without_genome), snps=", ".join(map(str, without_genome[:10])))
        )

    files = {
        key: os.path.join(
            outdir, "{name}_{file}.txt".format(name=name, file="canSNPs" if key == "snp" else key)
        )
        for key in ("tree", "snp", "references")
    }
    lines = lambda rows: "".join("\t".join(row) + "\n" for row in rows)
    rows = [[text(value) for value in row[:-1]] + [typed_text(row[-1])] for row in snps]
    contents = {
        files["tree"]: lines(paths),
        files["snp"]: lines([CANSNP_COLUMNS] + rows),
        files["references"]: lines(
            [GENOME_SOURCE_COLUMNS] + [[text(value) for value in row] for row in references]
        ),
    }
    if scheme_genomes:
        """the genomes file is written only for a database that holds genomes"""
        files["genomes"] = os.path.join(outdir, "{name}_genomes.txt".format(name=name))
        contents[files["genomes"]] = lines([GENOMES_COLUMNS] + scheme_genomes)
    write_outputs(contents)
    typed = sum(1 for row in snps if typed_text(row[-1]) == "yes")
    logger.info(
        "Exported {nodes} nodes, {snps} canSNPs ({typed} typed) and {refs} reference genomes to {outdir}".format(
            nodes=len(paths), snps=len(snps), typed=typed, refs=len(references), outdir=outdir
        )
    )
    return files
