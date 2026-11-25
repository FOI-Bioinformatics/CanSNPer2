"""cansnper2 ref get --local-references (CanSNPer2-download --local_references): the representative genomes of the genomes table of
a scheme are fetched into <outdir>/local/<accession>.fna. The genomes file names only the
accession, so the assembly name of the NCBI path is read from the folder listing."""
import gzip
import os
import sqlite3

import helpers
from CanSNPer2.DownloadGenomes import DownloadGenomes, assembly_from_listing
from CanSNPer2.modules.SchemeGenomes import load_genomes

LISTING = '''<html><body><a href="../">../</a>
<a href="GCA_000000011.1_ASM11v2/">GCA_000000011.1_ASM11v2/</a> 2020-01-01
<a href="GCA_000000011.2_ASM11v3/">GCA_000000011.2_ASM11v3/</a> 2021-01-01
</body></html>'''
GENOME = gzip.compress(b">c1\n" + b"ACGT" * 100 + b"\n")


def test_assembly_name_from_the_folder_listing():
	assert assembly_from_listing(LISTING, "GCA_000000011.1") == "ASM11v2"
	assert assembly_from_listing(LISTING, "GCA_000000011.2") == "ASM11v3"
	assert assembly_from_listing(LISTING, "GCA_000000099.1") is None


def database(tmp_path):
	db = helpers.make_db(tmp_path / "s.db", [("A.1", "chr", 5, "A", "G"), ("B.1", "chr", 9, "C", "T")])
	conn = sqlite3.connect(db)
	conn.execute("INSERT INTO tree (parent, child, rank_i) VALUES (1, 1, 1), (1, 2, 1)")
	(tmp_path / "g.txt").write_text("genome\tnode\trole\nGCA_000000011.1\tA.1\trepresentative\n"
		"GCA_000000012.1\tA.1\tmember\nGCF_000000013.1\tB.1\trepresentative\n")
	load_genomes(conn, str(tmp_path / "g.txt"))
	conn.close()
	return db


def downloader(tmp_path, db):
	obj = DownloadGenomes.__new__(DownloadGenomes)
	obj.source, obj.directory, obj.database_path = "genbank", str(tmp_path / "refs"), db
	obj.links = []
	def retrieve(link, target, timeout=None):
		obj.links.append(link)
		data = GENOME if link.endswith(".gz") else (LISTING.replace("GCA_000000011", "GCF_000000013") if "GCF" in link else LISTING).encode()
		with open(target, "wb") as out:
			out.write(data)
	obj.retrieve = retrieve
	obj.rename_to_genbank = lambda link, path: None
	return obj


def test_representatives_are_fetched_into_the_local_folder(tmp_path):
	db = database(tmp_path)
	obj = downloader(tmp_path, db)
	assert obj.run_local() == 0
	local = tmp_path / "refs" / "local"
	assert sorted(p.name for p in local.glob("*.fna")) == ["GCA_000000011.1.fna", "GCF_000000013.1.fna"], "members are not fetched"
	assert open(local / "GCA_000000011.1.fna").read().startswith(">c1")
	gz = [l for l in obj.links if l.endswith(".gz")]
	assert "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/000/011/GCA_000000011.1_ASM11v2/GCA_000000011.1_ASM11v2_genomic.fna.gz" in gz
	assert any(l.startswith("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/000/013/") for l in gz), "GCF accessions from the RefSeq tree"


def test_a_scheme_without_genomes_has_nothing_to_fetch(tmp_path):
	db = helpers.make_db(tmp_path / "plain.db", [("A.1", "chr", 5, "A", "G")])
	assert downloader(tmp_path, db).run_local() == 1


def test_an_accession_without_a_listed_assembly_is_a_failure(tmp_path):
	db = database(tmp_path)
	obj = downloader(tmp_path, db)
	plain = obj.retrieve
	def no_listing(link, target, timeout=None):
		if "GCF" in link and not link.endswith(".gz"):
			with open(target, "w") as out:
				out.write("<html></html>")
			return
		plain(link, target, timeout)
	obj.retrieve = no_listing
	assert obj.run_local() == 1
	assert (tmp_path / "refs" / "local" / "GCA_000000011.1.fna").exists(), "the others are still fetched"


# fetched while typing (CanSNPer2 --fetch_local)

from CanSNPer2.DownloadGenomes import LocalFetcher


def fetcher(tmp_path, listing=LISTING):
	links = []
	def retrieve(self, link, target, timeout=None):
		links.append(link)
		data = GENOME if link.endswith(".gz") else listing.encode()
		with open(target, "wb") as out:
			out.write(data)
	fetch = LocalFetcher(str(tmp_path / "refs"))
	fetch.retrieve = retrieve
	return fetch, links


def test_the_fetcher_downloads_the_genomes_asked_for(tmp_path, monkeypatch):
	fetch, links = fetcher(tmp_path)
	monkeypatch.setattr(DownloadGenomes, "retrieve", fetch.retrieve)
	monkeypatch.setattr(DownloadGenomes, "rename_to_genbank", lambda self, link, path: None)
	assert fetch(["GCA_000000011.1"]) == {"GCA_000000011.1"}
	assert open(tmp_path / "refs" / "local" / "GCA_000000011.1.fna").read().startswith(">c1")
	assert fetch.fetched == ["GCA_000000011.1"]
	n = len(links)
	assert fetch(["GCA_000000011.1"]) == {"GCA_000000011.1"} and len(links) == n, "a genome on disk is not fetched again"


def test_a_genome_that_failed_is_not_tried_again_in_the_run(tmp_path, monkeypatch):
	fetch, links = fetcher(tmp_path, listing="<html></html>")
	monkeypatch.setattr(DownloadGenomes, "retrieve", fetch.retrieve)
	assert fetch(["GCA_000000011.1"]) == set()
	assert "GCA_000000011.1" in fetch.failed
	n = len(links)
	assert fetch(["GCA_000000011.1"]) == set() and len(links) == n
	assert not (tmp_path / "refs" / "local" / "GCA_000000011.1.fna").exists()
