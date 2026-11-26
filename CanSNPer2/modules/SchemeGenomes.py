"""
The genomes of a scheme: which genome sits at which node of the tree, and in what role
(<scheme>_genomes.txt of CanSNPer2-data, columns genome, node, role; further columns such as
basis and source are not stored). The terminal nodes list their representative genomes and
the members folded into them; an outgroup was aligned with the scheme but sits at no node
(its node is empty), as CanSNPer2-data validate.py requires. It is stored with the node id
NO_NODE.

The rows are kept in the genomes table, whose id is the node id; the role is a column added
to it. Read typing with --local_references takes its candidate genomes from the
representatives. A database without the table, without the role column or without rows
holds no genomes, which is the case for every database built before 2.6.0.
"""

import sqlite3

from CanSNPer2.modules.DatabaseConnection import SourceFileError, open_source

GENOMES_COLUMNS = ("genome", "node", "role")
ROLES = ("representative", "member", "outgroup")
NO_NODE = -1  ## node id of a genome that sits at no node (an outgroup)


def read_genomes_file(path):
    """[(genome, node, role)] of a genomes file, in file order; SourceFileError for a file
    without the header or with a row that has no genome, node or role"""
    rows = []
    with open_source(path) as fin:
        header = fin.readline().lstrip("#").strip().split("\t")
        if tuple(header[:3]) != GENOMES_COLUMNS:
            raise SourceFileError(
                "The genomes file {path} must start with the header genome<TAB>node<TAB>role, "
                "it starts with: {h}".format(path=path, h="\t".join(header))
            )
        for number, line in enumerate(fin, 2):
            if not line.strip():
                continue
            fields = [f.strip() for f in line.rstrip("\n\r").split("\t")]
            fields += [""] * (3 - len(fields))
            if not fields[0] or not fields[2]:
                raise SourceFileError(
                    "The genomes file {path} line {n} needs a genome and a role: {line}".format(
                        path=path, n=number, line=line.strip()
                    )
                )
            rows.append(tuple(fields[:3]))
    return rows


def check_genomes(rows, nodes, path="the genomes file"):
    """the problems of genome rows against the node names of the tree"""
    problems = []
    seen = set()
    for genome, node, role in rows:
        if role == "outgroup" and node:
            problems.append(
                "{path}: genome {g} is an outgroup, which sits at no node, but is placed at {n}".format(
                    path=path, g=genome, n=node
                )
            )
        elif role != "outgroup" and not node:
            problems.append(
                "{path}: genome {g} sits at no node; only an outgroup does".format(path=path, g=genome)
            )
        elif node and node not in nodes:
            problems.append(
                "{path}: genome {g} is placed at node {n}, which is not in the tree".format(
                    path=path, g=genome, n=node
                )
            )
        if role not in ROLES:
            problems.append(
                "{path}: genome {g} has the role {r}, which is not one of {known}".format(
                    path=path, g=genome, r=role, known=", ".join(ROLES)
                )
            )
        if genome in seen:
            problems.append("{path}: genome {g} is listed more than once".format(path=path, g=genome))
        seen.add(genome)
    return problems


def load_genomes(conn, path):
    """check the genomes file against the tree of the database and store it; raises
    SourceFileError with every problem"""
    rows = read_genomes_file(path)
    node_ids = dict((name, i) for i, name in conn.execute("SELECT id, name FROM nodes"))
    problems = check_genomes(rows, node_ids, path)
    if problems:
        raise SourceFileError(
            "The genomes file cannot be loaded, {n} problem(s):\n  {list}".format(
                n=len(problems), list="\n  ".join(problems)
            )
        )
    conn.execute("""CREATE TABLE IF NOT EXISTS genomes (id integer NOT NULL, genome VARCHAR(255) NOT NULL UNIQUE,
		reference VARCHAR(255), FOREIGN KEY (id) REFERENCES nodes (id))""")
    columns = {row[1] for row in conn.execute("PRAGMA table_info(genomes)")}
    if "role" not in columns:
        conn.execute("ALTER TABLE genomes ADD COLUMN role TEXT")
    conn.executemany(
        "INSERT INTO genomes (id, genome, role) VALUES (?, ?, ?)",
        [(node_ids[node] if node else NO_NODE, genome, role) for genome, node, role in rows],
    )
    conn.commit()
    return len(rows)


def read_genomes(database):
    """[(genome, node, role)] of a database, in the order they were loaded; [] when it holds
    none. database is a path or an open sqlite3 connection"""
    conn = database if isinstance(database, sqlite3.Connection) else sqlite3.connect(str(database))
    try:
        tables = {row[0] for row in conn.execute("SELECT name FROM sqlite_master WHERE type = 'table'")}
        if "genomes" not in tables:
            return []
        columns = {row[1] for row in conn.execute("PRAGMA table_info(genomes)")}
        if "role" not in columns:
            return []
        return [
            tuple(row)
            for row in conn.execute(
                "SELECT g.genome, coalesce(n.name, ''), g.role FROM genomes g LEFT JOIN nodes n ON n.id = g.id ORDER BY g.rowid"
            )
        ]
    finally:
        if conn is not database:
            conn.close()
