"""Fetching a scheme from CanSNPer2-data and building its database (CanSNPer2-download --scheme).

The repository is replaced by a dictionary of files, so no test touches the network."""
import hashlib
import json
import os
import sqlite3

import pytest

from CanSNPer2.modules import FetchScheme
from CanSNPer2.modules.FetchScheme import (API_URL, DATA_REPO, RAW_URL, NotFound, SchemeError,
                                           GitHubSchemeSource, LocalSchemeSource, build_scheme, provenance_path)

COMMIT = "a" * 40

SNP = ("snp_id\tstrain\treference\tgenome\tposition\tderived_base\tancestral_base\n"
       "A.1\ttest\tTest2026\tREF1\t100\tT\tC\n"
       "B.1\ttest\tTest2026\tREF1\t200\tG\tA\n")
TREE = "R\nR\tA.1\nR\tA.1\tB.1\n"
REFERENCES = ("genome\tstrain\tgenbank_id\trefseq_id\tassembly_name\n"
              "REF1\tref one\tGCA_000000001.1\tGCF_000000001.1\tASM1v1\n")


def sha(text):
	return hashlib.sha256(text.encode()).hexdigest()


def scheme_repo(name="test_scheme", status="verified", manifest=True, snp=SNP):
	folder = "schemes/{n}/{n}".format(n=name)
	files = {
		folder + "_snp.txt": snp,
		folder + "_tree.txt": TREE,
		folder + "_references.txt": REFERENCES,
		"schemes/{n}/README.md".format(n=name): "# Test\n\n| Field | Value |\n|---|---|\n| Status | {s} |\n".format(s=status),
		"schemes.toml": '[taxo]\nkind = "taxonomy"\n',
	}
	if manifest:
		files["database/{n}.manifest.json".format(n=name)] = json.dumps({"sources": {
			folder + "_snp.txt": sha(SNP),
			folder + "_tree.txt": sha(TREE),
			folder + "_references.txt": sha(REFERENCES)}})
	return files


class FakeRepo(object):
	'''answers the GitHub URLs GitHubSchemeSource requests from a dictionary of files'''
	def __init__(self, files, commit=COMMIT, branches=("master", "dev")):
		self.files, self.commit, self.branches = files, commit, branches
		self.requested = []

	def __call__(self, url, accept=None, timeout=None):
		self.requested.append(url)
		for branch in self.branches:
			if url == API_URL.format(repo=DATA_REPO, endpoint="commits/" + branch):
				return self.commit.encode()
		if url.startswith(API_URL.format(repo=DATA_REPO, endpoint="git/trees/")):
			tree = [{"path": p, "type": "blob"} for p in self.files]
			return json.dumps({"tree": tree}).encode()
		prefix = RAW_URL.format(repo=DATA_REPO, ref=self.commit, path="")
		if url.startswith(prefix) and url[len(prefix):] in self.files:
			return self.files[url[len(prefix):]].encode()
		raise NotFound(url)


def source(files, **kwargs):
	return GitHubSchemeSource(ref="master", fetcher=FakeRepo(files, **kwargs))


def test_build_scheme_makes_a_database_and_records_where_it_came_from(tmp_path):
	db = str(tmp_path / "out" / "test_scheme.db")
	fake = FakeRepo(scheme_repo())
	provenance = build_scheme("test_scheme", db, source=GitHubSchemeSource(ref="master", fetcher=fake))
	conn = sqlite3.connect(db)
	assert conn.execute("SELECT snp_id,position FROM snp_annotation ORDER BY snp_id").fetchall() == [("A.1", 100), ("B.1", 200)]
	assert conn.execute("SELECT genome,genbank_id FROM snp_references").fetchall() == [("REF1", "GCA_000000001.1")]
	conn.close()
	recorded = json.load(open(provenance_path(db)))
	assert recorded == provenance
	assert recorded["commit"] == COMMIT and recorded["status"] == "verified"
	assert recorded["sources"]["schemes/test_scheme/test_scheme_snp.txt"] == sha(SNP)
	# every file is read at the resolved commit, never from the moving branch
	raw = [u for u in fake.requested if u.startswith("https://raw.")]
	assert raw and all("/" + COMMIT + "/" in u for u in raw)
	assert sorted(os.listdir(tmp_path / "out")) == ["test_scheme.db", "test_scheme.source.json"]


def test_a_file_that_does_not_match_the_manifest_builds_nothing(tmp_path):
	db = str(tmp_path / "test_scheme.db")
	damaged = scheme_repo(snp=SNP.replace("200", "201"))
	with pytest.raises(SchemeError, match="does not match the SHA-256"):
		build_scheme("test_scheme", db, source=source(damaged))
	assert os.listdir(tmp_path) == []


def test_a_scheme_without_a_manifest_is_built_with_a_warning(tmp_path, caplog):
	db = str(tmp_path / "test_scheme.db")
	build_scheme("test_scheme", db, source=source(scheme_repo(manifest=False)))
	assert os.path.exists(db)
	assert "cannot be checked" in caplog.text


def test_an_existing_database_is_kept_unless_overwrite_is_given(tmp_path):
	db = tmp_path / "test_scheme.db"
	db.write_text("keep me")
	with pytest.raises(SchemeError, match="--overwrite"):
		build_scheme("test_scheme", str(db), source=source(scheme_repo()))
	assert db.read_text() == "keep me"
	build_scheme("test_scheme", str(db), overwrite=True, source=source(scheme_repo()))
	assert sqlite3.connect(str(db)).execute("SELECT count(*) FROM snp_annotation").fetchone()[0] == 2


def test_a_draft_scheme_is_built_with_a_warning(tmp_path, caplog):
	db = str(tmp_path / "test_scheme.db")
	provenance = build_scheme("test_scheme", db, source=source(scheme_repo(status="draft")))
	assert provenance["status"] == "draft"
	assert "is a draft" in caplog.text


@pytest.mark.parametrize("name, message", [
	("missing", "There is no scheme missing"),
	("taxo", "taxonomy scheme"),
	("../escape", "not a valid scheme name"),
])
def test_unusable_scheme_names_are_refused(tmp_path, name, message):
	with pytest.raises(SchemeError, match=message):
		build_scheme(name, str(tmp_path / "x.db"), source=source(scheme_repo()))
	assert os.listdir(tmp_path) == []


def test_list_schemes_gives_snp_schemes_and_their_status():
	files = scheme_repo()
	files.update(scheme_repo(name="other", status="draft"))
	files["schemes/taxo/taxo_tree.txt"] = TREE
	assert source(files).list_schemes() == [("other", "draft"), ("test_scheme", "verified")]


def test_an_unknown_branch_is_reported():
	with pytest.raises(SchemeError, match="is not a branch, tag or commit"):
		GitHubSchemeSource(ref="nosuch", fetcher=FakeRepo(scheme_repo()))


def test_without_the_github_api_the_branch_is_read_directly(tmp_path, caplog):
	files = scheme_repo()
	def fetcher(url, accept=None, timeout=None):
		if url.startswith("https://api.github.com/"):
			raise SchemeError("rate limit exceeded")
		prefix = RAW_URL.format(repo=DATA_REPO, ref="dev", path="")
		if url.startswith(prefix) and url[len(prefix):] in files:
			return files[url[len(prefix):]].encode()
		raise NotFound(url)
	provenance = build_scheme("test_scheme", str(tmp_path / "t.db"), source=GitHubSchemeSource(ref="dev", fetcher=fetcher))
	assert provenance["commit"] is None and provenance["ref"] == "dev"
	assert "commit is not recorded" in caplog.text


# --repo: another GitHub repository with the CanSNPer2-data layout

def test_another_repository_is_read_at_its_own_urls(tmp_path):
	files = scheme_repo()
	def fetcher(url, accept=None, timeout=None):
		if url == API_URL.format(repo="someone/my-schemes", endpoint="commits/main"):
			return COMMIT.encode()
		prefix = RAW_URL.format(repo="someone/my-schemes", ref=COMMIT, path="")
		if url.startswith(prefix) and url[len(prefix):] in files:
			return files[url[len(prefix):]].encode()
		raise NotFound(url)
	src = GitHubSchemeSource(ref="main", repo="someone/my-schemes", fetcher=fetcher)
	provenance = build_scheme("test_scheme", str(tmp_path / "t.db"), source=src)
	assert provenance["repository"] == "https://github.com/someone/my-schemes"
	assert provenance["commit"] == COMMIT and provenance["ref"] == "main"


@pytest.mark.parametrize("repo", ["my-schemes", "https://github.com/a/b", "a/b/c"])
def test_a_repository_that_is_not_owner_slash_name_is_refused(repo):
	with pytest.raises(SchemeError, match="owner/name"):
		GitHubSchemeSource(repo=repo, fetcher=FakeRepo({}))


# --data_dir: a local directory with the CanSNPer2-data layout

def write_tree(root, files):
	for path, content in files.items():
		full = root / path
		full.parent.mkdir(parents=True, exist_ok=True)
		full.write_text(content)
	return root


def test_a_local_scheme_needs_only_its_three_files(tmp_path):
	data = write_tree(tmp_path / "data", {
		"schemes/custom/custom_snp.txt": SNP,
		"schemes/custom/custom_tree.txt": TREE,
		"schemes/custom/custom_references.txt": REFERENCES,
	})
	src = LocalSchemeSource(str(data))
	assert src.list_schemes() == [("custom", None)]
	db = str(tmp_path / "custom.db")
	provenance = build_scheme("custom", db, source=src)
	assert sqlite3.connect(db).execute("SELECT count(*) FROM snp_annotation").fetchone()[0] == 2
	assert provenance["repository"] == str(data) and provenance["commit"] is None
	assert provenance["sources"]["schemes/custom/custom_tree.txt"] == sha(TREE)


def test_a_local_file_changed_after_its_manifest_is_built_with_a_warning(tmp_path, caplog):
	files = scheme_repo()
	files["schemes/test_scheme/test_scheme_snp.txt"] = SNP.replace("200", "201")
	src = LocalSchemeSource(str(write_tree(tmp_path / "data", files)))
	build_scheme("test_scheme", str(tmp_path / "t.db"), source=src)
	assert "changed after the manifest was written" in caplog.text
	assert sqlite3.connect(str(tmp_path / "t.db")).execute(
		"SELECT position FROM snp_annotation WHERE snp_id='B.1'").fetchone() == (201,)


def test_a_local_git_checkout_records_its_commit_and_uncommitted_files(tmp_path, caplog):
	git = pytest.importorskip("shutil").which("git")
	if not git:
		pytest.skip("git is not installed")
	import subprocess
	data = write_tree(tmp_path / "data", scheme_repo())
	run = lambda *a: subprocess.run(["git", "-C", str(data)] + list(a), check=True, capture_output=True, text=True)
	run("init", "-q")
	run("add", ".")
	run("-c", "user.name=t", "-c", "user.email=t@t", "commit", "-q", "-m", "schemes")
	head = run("rev-parse", "HEAD").stdout.strip()
	provenance = build_scheme("test_scheme", str(tmp_path / "a.db"), source=LocalSchemeSource(str(data)))
	assert provenance["commit"] == head and provenance["uncommitted_sources"] == []
	(data / "schemes/test_scheme/test_scheme_tree.txt").write_text(TREE + "R\tA.1\tC.1\n")
	provenance = build_scheme("test_scheme", str(tmp_path / "b.db"), source=LocalSchemeSource(str(data)))
	assert provenance["uncommitted_sources"] == ["schemes/test_scheme/test_scheme_tree.txt"]
	assert "uncommitted changes" in caplog.text


def test_a_missing_local_scheme_names_the_list_command(tmp_path):
	src = LocalSchemeSource(str(write_tree(tmp_path / "data", scheme_repo())))
	with pytest.raises(SchemeError, match="--list_schemes --data_dir"):
		build_scheme("missing", str(tmp_path / "m.db"), source=src)


def test_data_dir_must_be_a_directory(tmp_path):
	with pytest.raises(SchemeError, match="is not a directory"):
		LocalSchemeSource(str(tmp_path / "nowhere"))


def test_the_command_line_builds_a_local_scheme_and_refuses_mixed_options(tmp_path):
	import subprocess, sys
	data = write_tree(tmp_path / "data", scheme_repo())
	env = dict(os.environ, PYTHONPATH=os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
	command = [sys.executable, "-m", "CanSNPer2.DownloadGenomes", "--logs", str(tmp_path / "logs")]
	listed = subprocess.run(command + ["--list_schemes", "--data_dir", str(data)], capture_output=True, text=True, env=env, cwd=str(tmp_path))
	assert listed.returncode == 0, listed.stderr
	assert "test_scheme" in listed.stdout and "verified" in listed.stdout
	built = subprocess.run(command + ["--scheme", "test_scheme", "--data_dir", str(data), "--no_references"],
		capture_output=True, text=True, env=env, cwd=str(tmp_path))
	assert built.returncode == 0, built.stderr
	assert (tmp_path / "test_scheme.db").exists() and (tmp_path / "test_scheme.source.json").exists()
	mixed = subprocess.run(command + ["--scheme", "test_scheme", "--data_dir", str(data), "--branch", "dev"],
		capture_output=True, text=True, env=env, cwd=str(tmp_path))
	assert mixed.returncode != 0 and "cannot be combined" in mixed.stderr


def test_a_scheme_folder_inside_a_larger_git_repository_reports_its_own_paths(tmp_path):
	import shutil, subprocess
	if not shutil.which("git"):
		pytest.skip("git is not installed")
	project = tmp_path / "project"
	data = write_tree(project / "sub" / "data", scheme_repo())
	run = lambda *a: subprocess.run(["git", "-C", str(project)] + list(a), check=True, capture_output=True, text=True)
	run("init", "-q")
	run("add", "sub/data/schemes/test_scheme/test_scheme_snp.txt")
	run("-c", "user.name=t", "-c", "user.email=t@t", "commit", "-q", "-m", "snp only")
	(data / "schemes/test_scheme/test_scheme_snp.txt").write_text(SNP + "C.1\ttest\tTest2026\tREF1\t300\tA\tG\n")
	provenance = build_scheme("test_scheme", str(tmp_path / "t.db"), source=LocalSchemeSource(str(data)))
	assert provenance["uncommitted_sources"] == [
		"schemes/test_scheme/test_scheme_references.txt",
		"schemes/test_scheme/test_scheme_snp.txt",
		"schemes/test_scheme/test_scheme_tree.txt"]


def test_without_a_branch_the_default_branch_of_the_repository_is_used():
	class MainRepo(FakeRepo):
		def __call__(self, url, accept=None, timeout=None):
			if url == "https://api.github.com/repos/someone/my-schemes":
				return json.dumps({"default_branch": "main"}).encode()
			if url == API_URL.format(repo="someone/my-schemes", endpoint="commits/main"):
				return COMMIT.encode()
			raise NotFound(url)
	src = GitHubSchemeSource(ref=None, repo="someone/my-schemes", fetcher=MainRepo({}))
	assert src.ref == "main" and src.commit == COMMIT
	with pytest.raises(SchemeError, match="was not found on GitHub"):
		GitHubSchemeSource(ref=None, repo="someone/missing", fetcher=MainRepo({}))
