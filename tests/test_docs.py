"""The documentation stays consistent with the program: every cansnper2 command line in a
code block of README.md and docs/*.md parses with the current parsers, every relative link
and anchor resolves, and the commands before 2.6.0 and their option spellings appear only in
the pages about upgrading."""

import pathlib
import re
import shlex

import pytest

from CanSNPer2.cli import build_parsers

ROOT = pathlib.Path(__file__).resolve().parent.parent
PAGES = [ROOT / "README.md", ROOT / "CHANGELOG.md"] + sorted((ROOT / "docs").glob("*.md"))
"""pages that describe the commands before 2.6.0 on purpose"""
HISTORY = {"upgrading.md", "CHANGELOG.md", "cli.md"}
"""a command line followed by this comment is an example of an error and is not parsed"""
FAILS = "# fails"


def code_blocks(text):
    """(language, body) of every fenced code block"""
    return re.findall(r"^```(\w*)\n(.*?)^```", text, re.M | re.S)


def command_lines(text):
    """every cansnper2 command of the shell code blocks, with continuation lines joined"""
    for language, body in code_blocks(text):
        if language not in ("sh", "bash", "console", ""):
            continue
        body = re.sub(r"\\\n\s*", " ", body)
        for line in body.splitlines():
            line = line.strip()
            if line.startswith("$ "):
                line = line[2:]
            if line.startswith("cansnper2 ") or line == "cansnper2":
                yield line


def examples():
    for page in PAGES:
        for line in command_lines(page.read_text()):
            yield pytest.param(line, id="{0}: {1}".format(page.name, line[:60]))


@pytest.mark.parametrize("line", list(examples()))
def test_command_lines_in_the_documentation_parse(line):
    if FAILS in line:
        pytest.skip("an example of an error")
    words = shlex.split(line, comments=True)
    for redirect in (">", "|", "2>", "&&"):
        if redirect in words:
            words = words[: words.index(redirect)]
    parser, _parsers = build_parsers()
    try:
        parser.parse_args(words[1:])
    except SystemExit as stop:
        assert stop.code == 0, "does not parse: {0}".format(line)


def slug(heading):
    """the anchor GitHub gives a heading"""
    text = re.sub(r"[`*_]", "", heading.strip().lower())
    text = re.sub(r"[^\w\- ]", "", text)
    return text.replace(" ", "-")


def anchors(page):
    text = re.sub(r"^```.*?^```", "", page.read_text(), flags=re.M | re.S)
    return {slug(match) for match in re.findall(r"^#+ (.+)$", text, re.M)}


def links():
    for page in PAGES:
        text = re.sub(r"^```.*?^```", "", page.read_text(), flags=re.M | re.S)
        for target in re.findall(r"\]\(([^)\s]+)\)", text):
            if re.match(r"[a-z]+:", target):
                continue
            yield pytest.param(page, target, id="{0} -> {1}".format(page.name, target))


@pytest.mark.parametrize("page,target", list(links()))
def test_relative_links_resolve(page, target):
    path, _, anchor = target.partition("#")
    destination = (page.parent / path).resolve() if path else page
    assert destination.exists(), "{0} links to a missing file {1}".format(page.name, path)
    if anchor and destination.suffix == ".md":
        assert anchor in anchors(destination), "{0}: no heading for #{1} in {2}".format(
            page.name, anchor, destination.name
        )


def test_old_commands_and_spellings_only_in_the_upgrade_pages():
    old = re.compile(r"CanSNPer2-(?:database|download|test)\b|(?<![\w-])--[a-z]+_[a-z_]+")
    found = []
    for page in PAGES:
        if page.name in HISTORY:
            continue
        for number, line in enumerate(page.read_text().splitlines(), 1):
            if old.search(line):
                found.append("{0}:{1}: {2}".format(page.name, number, line.strip()[:100]))
    assert not found, "\n".join(found)


def test_us_spelling():
    british = re.compile(
        r"\b(colour|behaviour|recognis|organis(e|ed|ing|ation)|analys(e|ed|ing)\b|licence\b|favour)", re.I
    )
    found = [
        "{0}:{1}".format(page.name, number)
        for page in PAGES
        for number, line in enumerate(page.read_text().splitlines(), 1)
        if british.search(line)
    ]
    assert not found, "\n".join(found)
