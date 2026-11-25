"""Messages and help name options as cansnper2 spells them (2.6.0).

Since 2.6.0 only the cansnper2 command is installed. A message that tells the user to give
--allow_unverified or --mod_file names an option the verbs do not have, so every option named
in a user-facing string of the package must exist in the cansnper2 parsers."""

import argparse
import ast
import pathlib
import re

from CanSNPer2.cli import build_parsers

PACKAGE = pathlib.Path(__file__).resolve().parent.parent / "CanSNPer2"
OPTION = re.compile(r"(?<![\w-])--[a-z][a-z0-9]*(?:[-_][a-z0-9]+)*")
'''modes of the database command before 2.6.0; they are named only in messages that the
	verbs cannot reach (a verb sets exactly one mode)'''
OLD_MODES = {"--create", "--import", "--export", "--info", "--convert", "--remove"}
'''options of the external programs that messages quote (minimap2 -c --cs=long)'''
EXTERNAL = {"--cs"}
'''parsers of the commands from before 2.6.0, no longer reachable from an installed command'''
OLD_PARSERS = {("SNPDatabase.py", "build_parser"), ("CanSNPerTree.py", "build_parser"),
	("DownloadGenomes.py", "build_parser"), ("selftest.py", "command")}


MESSAGE_CALLS = {"error", "warning", "info", "critical", "exception", "print", "exit"}
HELP_KEYWORDS = {"help", "description", "epilog"}


def user_strings():
	'''(file, line, text) of every string constant the user can see: the arguments of log
		calls, print, exit, parser.error and raised exceptions, and the help, description and
		epilog texts of the cansnper2 parsers; OLD_PARSERS are left out'''
	for path in sorted(PACKAGE.rglob("*.py")):
		tree = ast.parse(path.read_text())
		old_parser = set()
		for node in ast.walk(tree):
			if isinstance(node, ast.FunctionDef) and (path.name, node.name) in OLD_PARSERS:
				old_parser.update(id(sub) for sub in ast.walk(node))
		seen = set()
		for node in ast.walk(tree):
			if id(node) in old_parser:
				continue
			roots = []
			if isinstance(node, ast.Raise) and node.exc is not None:
				roots.append(node.exc)
			elif isinstance(node, ast.Call):
				name = getattr(node.func, "attr", None) or getattr(node.func, "id", None)
				if name in MESSAGE_CALLS:
					roots.extend(node.args)
				roots.extend(k.value for k in node.keywords if k.arg in HELP_KEYWORDS)
			for root in roots:
				for sub in ast.walk(root):
					if isinstance(sub, ast.Constant) and isinstance(sub.value, str) and id(sub) not in seen:
						seen.add(id(sub))
						yield path.relative_to(PACKAGE.parent), sub.lineno, sub.value


def verb_options():
	_parser, parsers = build_parsers()
	options = {"--help", "--version"}
	for sub in parsers.values():
		for action in sub._actions:
			options.update(action.option_strings)
	return options


def test_no_option_with_an_underscore_is_named():
	found = ["{0}:{1}: {2}".format(path, line, name) for path, line, text in user_strings()
		for name in OPTION.findall(text) if "_" in name]
	assert not found, "\n".join(found)


def test_every_option_named_exists_in_cansnper2():
	known = verb_options() | OLD_MODES | EXTERNAL
	found = ["{0}:{1}: {2}".format(path, line, name) for path, line, text in user_strings()
		for name in OPTION.findall(text) if name not in known]
	assert not found, "\n".join(found)


def test_no_removed_command_is_named_outside_the_upgrade_notes():
	'''the old command names may appear in the top-level help, which says they are no longer
		installed'''
	pattern = re.compile(r"CanSNPer2-(?:database|download|test)\b")
	found = ["{0}:{1}".format(path, line) for path, line, text in user_strings()
		if pattern.search(text) and "no longer installed" not in text]
	assert not found, "\n".join(found)


def test_log_files_are_named_after_the_verb(tmp_path):
	'''cansnper2 db info writes cansnper2-db-info-<date>.log; the name has no colon, which
		file names on Windows cannot contain'''
	import datetime
	from CanSNPer2.modules.LogSetup import get_logpath, log_tag
	assert log_tag("cansnper2 scheme get") == "cansnper2-scheme-get"
	now = datetime.datetime(2026, 11, 25, 13, 5)
	first = get_logpath(str(tmp_path), "cansnper2-db-info", now=now)
	assert pathlib.Path(first).name == "cansnper2-db-info-Nov-25-2026.log"
	open(first, "w").close()
	second = get_logpath(str(tmp_path), "cansnper2-db-info", now=now)
	assert ":" not in pathlib.Path(second).name and second != first


def test_scheme_and_ref_verbs_take_the_shared_log_levels():
	'''--quiet had no effect in scheme get and ref get, and -v printed debug output there'''
	import logging
	from CanSNPer2.DownloadGenomes import log_level
	_parser, parsers = build_parsers()
	for verb in ("scheme get", "ref get"):
		extra = ["NAME"] if verb == "scheme get" else []
		parse = lambda *options: parsers[verb].parse_args(extra + list(options))
		assert log_level(parse()) == logging.WARNING, verb
		assert log_level(parse("--quiet")) == logging.ERROR, verb
		assert log_level(parse("-v")) == logging.INFO, verb
		assert log_level(parse("--debug")) == logging.DEBUG, verb
