"""docs/cli.md, the command line reference, is generated from the cansnper2 parsers by
tools/gen_cli_docs.py and has to stay current with them."""

import importlib.util
import pathlib
import sys

import pytest

ROOT = pathlib.Path(__file__).resolve().parent.parent
REFERENCE = ROOT / "docs" / "cli.md"


def generator():
    spec = importlib.util.spec_from_file_location("gen_cli_docs", ROOT / "tools" / "gen_cli_docs.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_every_option_of_every_verb_is_in_the_reference():
    from CanSNPer2.cli import build_parsers

    text = REFERENCE.read_text()
    _top, parsers = build_parsers()
    for verb, parser in parsers.items():
        assert "## cansnper2 " + verb in text, verb
        for action in parser._actions:
            for option in action.option_strings:
                assert option in text, "{0} {1}".format(verb, option)


def test_the_reference_is_current(monkeypatch):
    gen = generator()
    if sys.version_info[:2] != gen.PYTHON:
        pytest.skip("docs/cli.md is generated with Python {0}.{1}".format(*gen.PYTHON))
    monkeypatch.setenv("COLUMNS", gen.WIDTH)
    assert REFERENCE.read_text() == gen.render(), (
        "docs/cli.md is out of date, run: python tools/gen_cli_docs.py"
    )
