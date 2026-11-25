#!/usr/bin/env python3
'''Write docs/cli.md, the command line reference, from the cansnper2 parsers.

The help of argparse is the single source: run this after changing an option or a help
text, and commit the result. tests/test_cli_reference.py fails when docs/cli.md is out of
date. argparse formats the help slightly differently between Python versions, so the file
is generated, and compared exactly, with the Python version named in PYTHON.'''

import argparse
import contextlib
import io
import os
import pathlib
import sys

PYTHON = (3, 13)
WIDTH = "80"
OUTPUT = pathlib.Path(__file__).resolve().parent.parent / "docs" / "cli.md"
'''verbs in the order of the reference'''
ORDER = ["type", "selftest", "scheme list", "scheme get", "ref get", "db create", "db import",
	"db export", "db info", "db convert", "db modify", "db remove"]

HEADER = """# Command line reference

This page lists the help of every `cansnper2` command, as printed by
`cansnper2 <command> --help`. It is generated from the program by
`tools/gen_cli_docs.py`; do not edit it by hand.

Options that take a value show their default in parentheses. Paths are relative
to the current directory, except for `cansnper2 type`, where the output, log and
temporary folders are relative to `--workdir`.

"""


def help_text(parser):
	with contextlib.redirect_stdout(io.StringIO()) as out:
		parser.print_help()
	return out.getvalue().rstrip() + "\n"


def summary(top, verb):
	'''the one-line help of a verb, as the parent command lists it'''
	parser = top
	for word in verb.split():
		commands = next(action for action in parser._actions if isinstance(action, argparse._SubParsersAction))
		line = next(choice.help for choice in commands._choices_actions if choice.dest == word)
		parser = commands.choices[word]
	return line[0].upper() + line[1:]


def anchor(verb):
	return "cansnper2-" + verb.replace(" ", "-")


def render():
	os.environ["COLUMNS"] = WIDTH
	os.environ.pop("FORCE_COLOR", None)
	os.environ["NO_COLOR"] = "1"
	from CanSNPer2.cli import build_parsers
	top, parsers = build_parsers()
	commands = next(action for action in top._actions if isinstance(action, argparse._SubParsersAction))
	parsers = dict(parsers, selftest=commands.choices["selftest"])
	assert sorted(parsers) == sorted(ORDER), sorted(parsers)
	parts = [HEADER, "| Command | Purpose |\n| --- | --- |\n"]
	for verb in ORDER:
		parts.append("| [`cansnper2 {v}`](#{a}) | {h} |\n".format(v=verb, a=anchor(verb), h=summary(top, verb)))
	parts.append("\n## cansnper2\n\n```text\n" + help_text(top) + "```\n")
	for verb in ORDER:
		parts.append("\n## cansnper2 {v}\n\n```text\n{h}```\n".format(v=verb, h=help_text(parsers[verb])))
	return "".join(parts)


def main():
	if sys.version_info[:2] != PYTHON:
		sys.exit("generate docs/cli.md with Python {0}.{1}, argparse formats differently in other "
			"versions".format(*PYTHON))
	OUTPUT.parent.mkdir(exist_ok=True)
	OUTPUT.write_text(render())
	print("wrote {0}".format(OUTPUT))


if __name__ == "__main__":
	main()
