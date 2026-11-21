'''
Fetch a SNP scheme from the CanSNPer2-data repository and build a database from it.

CanSNPer2-data keeps each scheme as three text files, schemes/<name>/<name>_snp.txt,
<name>_tree.txt and <name>_references.txt, and records the SHA-256 of each in
database/<name>.manifest.json. The scheme files are fetched rather than the committed
database, so the database is built by the CanSNPer2 version that will read it, and the
fetched files are checked against the hashes in the manifest before anything is built.

The branch is resolved to a commit first and every file is read at that commit, so a push
to the branch during the download cannot mix two versions of a scheme. The commit is
written to a provenance file beside the database.

master holds verified schemes only; dev also holds draft schemes (see the Status section of
schemes/README.md in CanSNPer2-data).

The same layout is read from another GitHub repository (a fork, or a repository of custom
schemes) with GitHubSchemeSource(repo="owner/name"), and from a local directory, such as a
clone or a folder of custom schemes, with LocalSchemeSource. Only schemes/<name>/ with its
three files is required; schemes.toml, the README with its Status and the manifest are used
when they are present.
'''

import datetime
import hashlib
import json
import logging
import os
import re
import subprocess
import tempfile
import tomllib
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

from CanSNPer2 import __version__

logger = logging.getLogger(__name__)

DATA_REPO = "FOI-Bioinformatics/CanSNPer2-data"
RAW_URL = "https://raw.githubusercontent.com/{repo}/{ref}/{path}"
API_URL = "https://api.github.com/repos/{repo}/{endpoint}"
FETCH_TIMEOUT = 60 ## seconds per connection
GITHUB_REPO = re.compile(r"^[A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+$")
SCHEME_NAME = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")
STATUS_LINE = re.compile(r"^\|\s*Status\s*\|\s*([A-Za-z]+)\s*\|", re.MULTILINE)


class SchemeError(Exception):
	'''A scheme that cannot be fetched, verified or built'''
	pass


class NotFound(SchemeError):
	'''The requested file does not exist at the requested commit'''
	pass


def fetch(url, accept=None, timeout=FETCH_TIMEOUT):
	'''Return the body of url as bytes. A 404 raises NotFound, any other failure SchemeError'''
	request = Request(url, headers={"User-Agent": "CanSNPer2/{v}".format(v=__version__)})
	if accept:
		request.add_header("Accept", accept)
	token = os.environ.get("GITHUB_TOKEN") or os.environ.get("GH_TOKEN")
	if token and url.startswith("https://api.github.com/"):
		## anonymous API use is limited to 60 requests an hour per address, which shared
		## machines such as CI runners reach; a token raises it. Only sent to the API
		request.add_header("Authorization", "Bearer {token}".format(token=token))
	try:
		with urlopen(request, timeout=timeout) as response:
			return response.read()
	except HTTPError as e:
		if e.code in (404, 422): ## the GitHub API answers 422 for an unknown branch
			raise NotFound("{url} was not found".format(url=url))
		if e.code in (403, 429) and url.startswith("https://api.github.com/"):
			raise SchemeError("Could not fetch {url}: {e}. The GitHub API allows 60 anonymous requests an hour, set GITHUB_TOKEN to a GitHub token to raise the limit".format(url=url, e=e))
		raise SchemeError("Could not fetch {url}: {e}".format(url=url, e=e))
	except (URLError, OSError) as e:
		raise SchemeError("Could not fetch {url}: {e}".format(url=url, e=e))


def sha256(data):
	return hashlib.sha256(data).hexdigest()


class SchemeSource(object):
	'''A collection of schemes laid out as CanSNPer2-data: schemes/<name>/<name>_snp.txt and so on,
		an optional schemes.toml and optional manifests under database/. Subclasses say where the
		files are read from (raw, paths) and how that place is recorded (location, provenance_fields)'''
	manifest_mismatch_is_error = True

	def __init__(self):
		self._config = None

	def raw(self, path):
		raise NotImplementedError

	def paths(self):
		'''every file path of the collection, relative to its top'''
		raise NotImplementedError

	def location(self):
		'''the collection, as it is named in messages'''
		raise NotImplementedError

	def list_hint(self):
		'''the command that lists the schemes of this collection'''
		raise NotImplementedError

	def provenance_fields(self, hashes):
		raise NotImplementedError

	def source_name(self):
		'''the collection as the scheme_info table records it: owner/name or a folder'''
		raise NotImplementedError

	def known_issues(self, name):
		'''The snp_ids that schemes.toml lists under known_issues for the scheme, as a
			comma-separated string. A schemes.toml that cannot be read, or that has no entry
			for the scheme, gives an empty string'''
		try:
			entry = self.config().get(name, {})
		except SchemeError as e:
			logger.warning("The known issues of {name} could not be read: {e}".format(name=name, e=e))
			return ""
		issues = entry.get("known_issues", []) if isinstance(entry, dict) else []
		if isinstance(issues, str):
			issues = issues.split(",")
		return ",".join(str(snp).strip() for snp in issues if str(snp).strip())

	def config(self):
		'''The entries of schemes.toml; a scheme without an entry follows the naming convention'''
		if self._config is None:
			try:
				self._config = tomllib.loads(self.raw("schemes.toml").decode())
			except NotFound:
				self._config = {}
			except (tomllib.TOMLDecodeError, UnicodeDecodeError) as e:
				raise SchemeError("schemes.toml of {loc} could not be read: {e}".format(loc=self.location(), e=e))
		return self._config

	def scheme_files(self, name):
		'''Paths of the files of an SNP scheme, with the overrides of schemes.toml'''
		if not SCHEME_NAME.match(name):
			raise SchemeError("{name} is not a valid scheme name".format(name=name))
		entry = self.config().get(name, {})
		kind = entry.get("kind", "snp")
		if kind != "snp":
			raise SchemeError("{name} is a {kind} scheme, CanSNPer2 can only use SNP schemes".format(name=name, kind=kind))
		folder = "schemes/{name}".format(name=name)
		database_dir = entry.get("database_dir", "database")
		return {
			"snp": entry.get("snp", "{folder}/{name}_snp.txt".format(folder=folder, name=name)),
			"tree": entry.get("tree", "{folder}/{name}_tree.txt".format(folder=folder, name=name)),
			"references": "{folder}/{name}_references.txt".format(folder=folder, name=name),
			"readme": entry.get("readme", "{folder}/README.md".format(folder=folder)),
			"manifest": "{dir}/{name}.manifest.json".format(dir=database_dir, name=name),
			"source_type": entry.get("source_type", "CanSNPer"),
		}

	def status(self, readme_path):
		'''draft or verified, as the scheme README gives it; None when it gives none'''
		try:
			text = self.raw(readme_path).decode(errors="replace")
		except NotFound:
			return None
		found = STATUS_LINE.search(text)
		return found.group(1).lower() if found else None

	def list_schemes(self):
		'''The SNP schemes of the collection, as (name, status) pairs'''
		paths = set(self.paths())
		names = set()
		for path in paths:
			parts = path.split("/")
			if len(parts) == 3 and parts[0] == "schemes" and parts[2] == "{n}_snp.txt".format(n=parts[1]):
				names.add(parts[1])
		for name, entry in self.config().items():
			if isinstance(entry, dict) and entry.get("kind", "snp") == "snp" and entry.get("snp") in paths:
				names.add(name)
		schemes = []
		for name in sorted(names):
			schemes.append((name, self.status(self.scheme_files(name)["readme"])))
		return schemes

	def download(self, name, directory):
		'''Copy the scheme files of name into directory and check them against the manifest.
			Returns the local paths and the provenance record'''
		files = self.scheme_files(name)
		local = {}
		hashes = {}
		for role in ("snp", "tree", "references"):
			path = files[role]
			try:
				data = self.raw(path)
			except NotFound:
				if role == "snp":
					raise SchemeError("There is no scheme {name} in {loc}, {hint} lists those there are".format(name=name, loc=self.location(), hint=self.list_hint()))
				raise SchemeError("The scheme {name} has no {path} in {loc}".format(name=name, path=path, loc=self.location()))
			local[role] = os.path.join(directory, os.path.basename(path))
			with open(local[role], "wb") as out:
				out.write(data)
			hashes[path] = sha256(data)
		self.verify(name, files["manifest"], hashes)
		status = self.status(files["readme"])
		if status == "draft":
			logger.warning("The scheme {name} is a draft: it is in {loc} but has not been verified yet".format(name=name, loc=self.location()))
		provenance = {
			"scheme": name,
			"status": status,
			"sources": hashes,
			"source_type": files["source_type"],
			"built": datetime.date.today().isoformat(),
			"cansnper2": __version__,
		}
		provenance.update(self.provenance_fields(hashes))
		return local, provenance

	def verify(self, name, manifest_path, hashes):
		'''Compare the files with the SHA-256 the manifest records for them. Fetched from GitHub,
			a difference means a damaged transfer or a scheme changed without a rebuild, and the
			files are not used. In a local directory it means the files were edited after the
			manifest was written, which is reported'''
		try:
			manifest = json.loads(self.raw(manifest_path).decode())
		except NotFound:
			logger.warning("{path} was not found, the files of {name} cannot be checked".format(path=manifest_path, name=name))
			return
		except ValueError as e:
			raise SchemeError("{path} could not be read: {e}".format(path=manifest_path, e=e))
		recorded = manifest.get("sources", {})
		matched = True
		for path, digest in sorted(hashes.items()):
			if path not in recorded:
				logger.warning("{manifest} records no hash for {path}, it is not checked".format(manifest=manifest_path, path=path))
			elif recorded[path] != digest:
				if self.manifest_mismatch_is_error:
					raise SchemeError("{path} does not match the SHA-256 in {manifest}: the download is damaged or the scheme was changed without being rebuilt".format(path=path, manifest=manifest_path))
				logger.warning("{path} does not match the SHA-256 in {manifest}, it was changed after the manifest was written".format(path=path, manifest=manifest_path))
				matched = False
		if matched:
			logger.info("The files of {name} match {manifest}".format(name=name, manifest=manifest_path))


class GitHubSchemeSource(SchemeSource):
	'''One commit of a GitHub repository laid out as CanSNPer2-data.

		ref: a branch (master or dev), tag or commit of the repository; None gives the
			default branch of the repository
		repo: owner/name on GitHub, CanSNPer2-data by default
	'''
	def __init__(self, ref="master", repo=DATA_REPO, fetcher=fetch):
		super().__init__()
		if not GITHUB_REPO.match(repo):
			raise SchemeError("{repo} is not a GitHub repository, give it as owner/name".format(repo=repo))
		self.repo = repo
		self.fetcher = fetcher
		self.ref = ref or self.default_branch()
		self.commit = self.resolve_commit()

	def default_branch(self):
		'''the default branch of the repository (master for CanSNPer2-data, often main elsewhere)'''
		try:
			return json.loads(self.fetcher(API_URL.format(repo=self.repo, endpoint="").rstrip("/")).decode())["default_branch"]
		except NotFound:
			raise SchemeError("The repository {repo} was not found on GitHub".format(repo=self.repo))
		except (SchemeError, ValueError, KeyError) as e:
			logger.warning("Could not look up the default branch of {repo}, master is used: {e}".format(repo=self.repo, e=e))
			return "master"

	def resolve_commit(self):
		'''The commit ref points to now. Without it (the GitHub API is rate limited for
			anonymous use) the files are read from the branch itself, which is reported'''
		url = API_URL.format(repo=self.repo, endpoint="commits/{ref}".format(ref=self.ref))
		try:
			commit = self.fetcher(url, accept="application/vnd.github.sha").decode().strip()
		except NotFound:
			raise SchemeError("{ref} is not a branch, tag or commit of {repo}, or the repository does not exist".format(ref=self.ref, repo=self.repo))
		except SchemeError as e:
			logger.warning("Could not resolve {ref} to a commit, the files are read from {ref} directly and the commit is not recorded: {e}".format(ref=self.ref, e=e))
			return None
		if not re.match(r"^[0-9a-f]{40}$", commit):
			logger.warning("Unexpected answer when resolving {ref} to a commit, the files are read from {ref} directly".format(ref=self.ref))
			return None
		logger.info("{repo} {ref} is at commit {commit}".format(repo=self.repo, ref=self.ref, commit=commit))
		return commit

	def raw(self, path):
		return self.fetcher(RAW_URL.format(repo=self.repo, ref=self.commit or self.ref, path=path))

	def paths(self):
		url = API_URL.format(repo=self.repo, endpoint="git/trees/{ref}?recursive=1".format(ref=self.commit or self.ref))
		try:
			tree = json.loads(self.fetcher(url).decode())
		except (SchemeError, ValueError) as e:
			raise SchemeError("Could not list the files of {loc}: {e}".format(loc=self.location(), e=e))
		return [item["path"] for item in tree.get("tree", []) if item.get("type") == "blob"]

	def location(self):
		return "{ref} of {repo}".format(ref=self.ref, repo=self.repo)

	def list_hint(self):
		repo = "" if self.repo == DATA_REPO else " --repo {repo}".format(repo=self.repo)
		return "CanSNPer2-download --list_schemes{repo} --branch {ref}".format(repo=repo, ref=self.ref)

	def source_name(self):
		return self.repo

	def provenance_fields(self, hashes):
		return {
			"repository": "https://github.com/{repo}".format(repo=self.repo),
			"ref": self.ref,
			"commit": self.commit,
		}


class LocalSchemeSource(SchemeSource):
	'''A local directory laid out as CanSNPer2-data, such as a clone of it or a directory of
		custom schemes. Files are read as they are on disk; when the directory is a git
		checkout, its commit and any scheme file with uncommitted changes are recorded'''
	manifest_mismatch_is_error = False

	def __init__(self, directory):
		super().__init__()
		if not os.path.isdir(directory):
			raise SchemeError("{d} is not a directory".format(d=directory))
		self.directory = os.path.abspath(directory)

	def raw(self, path):
		full = os.path.join(self.directory, path)
		try:
			with open(full, "rb") as handle:
				return handle.read()
		except FileNotFoundError:
			raise NotFound("{full} was not found".format(full=full))
		except OSError as e:
			raise SchemeError("Could not read {full}: {e}".format(full=full, e=e))

	def paths(self):
		found = []
		for root, dirs, files in os.walk(self.directory):
			dirs[:] = [d for d in dirs if not d.startswith(".")]
			for name in files:
				found.append(os.path.relpath(os.path.join(root, name), self.directory).replace(os.sep, "/"))
		return found

	def location(self):
		return self.directory

	def list_hint(self):
		return "CanSNPer2-download --list_schemes --data_dir {d}".format(d=self.directory)

	def git(self, *arguments):
		'''output of a git command in the directory, None when it is not a git checkout'''
		try:
			result = subprocess.run(["git", "-C", self.directory] + list(arguments), capture_output=True, text=True)
		except OSError:
			return None
		return result.stdout if result.returncode == 0 else None

	def source_name(self):
		return self.directory

	def provenance_fields(self, hashes):
		fields = {"repository": self.directory, "ref": None, "commit": None}
		commit = self.git("rev-parse", "HEAD")
		if commit:
			fields["commit"] = commit.strip()
			## changed since the commit (staged or not) and untracked, with paths relative to the
			## directory, which may be a subfolder of the repository
			changed = (self.git("diff", "--name-only", "--relative", "HEAD", "--", *sorted(hashes)) or "") + \
						(self.git("ls-files", "--others", "--", *sorted(hashes)) or "")
			fields["uncommitted_sources"] = sorted(set(line for line in changed.splitlines() if line.strip()))
			if fields["uncommitted_sources"]:
				logger.warning("{files} in {d} have uncommitted changes, the recorded commit does not fully describe the database".format(
						files=", ".join(fields["uncommitted_sources"]), d=self.directory))
		return fields


def provenance_path(database):
	'''<database without .db>.source.json beside the database'''
	base = database[:-3] if database.endswith(".db") else database
	return base + ".source.json"


def build_scheme(name, database, ref="master", overwrite=False, source=None):
	'''Fetch scheme name from CanSNPer2-data at ref and build it into database.
		The database is built under a temporary name and moved into place only when complete,
		so a failed build leaves no partial database behind'''
	from CanSNPer2.modules.InitiateDatabase import CanSNPDatabase
	if os.path.exists(database) and not overwrite:
		raise SchemeError("{db} already exists, use --overwrite to replace it".format(db=database))
	if source is None:
		source = GitHubSchemeSource(ref=ref)
	target_dir = os.path.dirname(os.path.abspath(database))
	os.makedirs(target_dir, exist_ok=True)
	with tempfile.TemporaryDirectory(prefix="cansnper2_scheme_") as tmpdir:
		local, provenance = source.download(name, tmpdir)
		scheme_info = {
			"scheme": name,
			"status": provenance["status"] or "unknown",
			"source": source.source_name(),
			"ref": provenance.get("ref"),
			"commit": provenance.get("commit") or "unknown",
			"known_issues": source.known_issues(name),
		}
		fd, partial = tempfile.mkstemp(suffix=".db.part", dir=target_dir)
		os.close(fd)
		os.remove(partial) ## the database loader creates the file itself
		try:
			CanSNPDatabase(database=partial, create=True, tree=local["tree"], annotation=local["snp"],
							references=local["references"], source_type=provenance["source_type"], verbose=False,
							scheme_info=scheme_info)
			os.replace(partial, database)
		finally:
			if os.path.exists(partial):
				os.remove(partial)
	with open(provenance_path(database), "w") as out:
		json.dump(provenance, out, indent=2, sort_keys=True)
		out.write("\n")
	logger.info("Built {db} from {name} in {loc} ({commit})".format(db=database, name=name,
				loc=source.location(), commit=provenance["commit"] or "commit unknown"))
	return provenance
