'''
The scheme_info table of a CanSNPer2 database: which scheme the database holds, the status of
that scheme, where its source files came from and which CanSNPer2 version built it.

	scheme_info(key TEXT PRIMARY KEY, value TEXT)

--create and --import always write format, built_with and built. The other keys are written
when they are known:

	scheme        the name of the scheme
	status        verified, draft or unknown
	source        the repository (owner/name) or folder the scheme files were read from
	ref           the branch, tag or commit asked for
	commit        the commit the files were read at, unknown when it could not be resolved
	known_issues  comma-separated snp_ids that schemes.toml lists as known issues

The table is optional to read: a database built before 2.3.0 or by hand has none, and is
used with a warning.
'''

import datetime
import logging
import os
import sqlite3
from urllib.request import pathname2url

from CanSNPer2 import __version__

FORMAT = "2.3"
STATUSES = ("verified", "draft", "unknown")
## the order in which --info lists the keys; keys not listed here follow in name order
KEY_ORDER = ("scheme", "status", "source", "ref", "commit", "known_issues", "format", "built_with", "built")
MISSING_TEXT = ("the database records no scheme information (built with CanSNPer2 before 2.3.0 or by hand); "
	"check where it came from")

TABLE_SQL = "CREATE TABLE IF NOT EXISTS scheme_info (key TEXT PRIMARY KEY, value TEXT)"


def build_info(**known):
	'''The keys --create and --import always write, together with the known keys that are
		not None'''
	info = {"format": FORMAT, "built_with": __version__, "built": datetime.date.today().isoformat()}
	info.update({key: str(value) for key, value in known.items() if value is not None})
	return info


def write_scheme_info(conn, info):
	'''Create the table when it is missing and write the keys of info into it, replacing the
		values a key already has. A value of None is not written.'''
	conn.execute(TABLE_SQL)
	conn.executemany("INSERT OR REPLACE INTO scheme_info (key, value) VALUES (?, ?)",
		[(key, str(value)) for key, value in info.items() if value is not None])
	conn.commit()


def read_scheme_info(path):
	'''The keys and values of the scheme_info table of the database at path, {} when the
		database has no such table. The database is opened read-only, so nothing is created
		or changed.'''
	if not os.path.exists(path):
		return {}
	## the path is quoted, so a name with ? or # is not read as part of the URI
	conn = sqlite3.connect("file:{path}?mode=ro".format(path=pathname2url(os.path.abspath(path))), uri=True)
	try:
		found = conn.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='scheme_info'").fetchone()
		if not found:
			return {}
		return {key: value for key, value in conn.execute("SELECT key, value FROM scheme_info")}
	finally:
		conn.close()


def ordered_keys(info):
	'''the keys of info in the order --info lists them'''
	return [key for key in KEY_ORDER if key in info] + sorted(key for key in info if key not in KEY_ORDER)


def known_issues(info):
	'''the snp_ids recorded as known issues, as a list'''
	return [snp for snp in (info.get("known_issues") or "").split(",") if snp.strip()]


def describe(info):
	'''One line on the scheme, for example
		brucella_suis (draft), FOI-Bioinformatics/CanSNPer2-data dev@abc1234, built 2026-11-21 with CanSNPer2 2.3.0'''
	parts = ["{scheme} ({status})".format(scheme=info.get("scheme") or "unnamed scheme",
		status=info.get("status") or "status not recorded")]
	origin = info.get("source") or ""
	ref, commit = info.get("ref"), info.get("commit")
	if ref:
		origin = "{origin} {ref}".format(origin=origin, ref=ref).strip()
	if commit and commit != "unknown":
		origin = "{origin}@{commit}".format(origin=origin, commit=commit[:7])
	elif commit == "unknown" and origin:
		origin = "{origin} (commit unknown)".format(origin=origin)
	if origin:
		parts.append(origin)
	built = "built"
	if info.get("built"):
		built += " " + info["built"]
	if info.get("built_with"):
		built += " with CanSNPer2 " + info["built_with"]
	parts.append(built)
	line = ", ".join(parts)
	issues = known_issues(info)
	if issues:
		line += "; known issues: {snps}".format(snps=", ".join(issues))
	return line


def scheme_notice(path):
	'''(level, text) of the notice logged when a database is opened. The level is
		logging.WARNING when the database has no scheme_info table, when the status is draft,
		unknown or not recorded, or when known issues are recorded; otherwise logging.INFO.'''
	info = read_scheme_info(path)
	if not info:
		return logging.WARNING, MISSING_TEXT
	text = "Scheme of {db}: {line}".format(db=os.path.basename(path), line=describe(info))
	status = info.get("status") or "unknown"
	level = logging.INFO
	if status == "draft":
		text += ". The scheme is a draft and has not been verified"
		level = logging.WARNING
	elif status != "verified":
		text += ". The status of the scheme is unknown; it has not been verified"
		level = logging.WARNING
	if known_issues(info):
		text += ". Calls at the SNPs listed as known issues may be unreliable"
		level = logging.WARNING
	return level, text


def log_scheme_notice(path, logger):
	'''Log the notice of scheme_notice once, at its level. A database that cannot be read
		here is left to the checks of the command.'''
	try:
		level, text = scheme_notice(path)
	except sqlite3.Error as e:
		logger.warning("The scheme information of {db} could not be read: {e}".format(db=path, e=e))
		return
	logger.log(level, text)
