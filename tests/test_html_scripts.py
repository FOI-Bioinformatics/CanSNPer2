"""The tree page and the results dashboard are single html files whose javascript is written
from Python string templates. A quote that Python unescapes, or a missing bracket, made the
browser reject the whole script, and both pages stayed empty; the content tests did not see it.
The scripts are checked with node --check when node is installed."""

import os
import re
import shutil
import subprocess

import pytest

from CanSNPer2.modules.NewickTree import NewickTree
from test_dashboard import REPORTS, make_db

NODE = shutil.which("node")
pytestmark = pytest.mark.skipif(NODE is None, reason="node is not installed")


def check_scripts(page, tmp_path):
    scripts = re.findall(r"<script[^>]*>(.*?)</script>", open(page).read(), re.S)
    assert scripts, "the page has a script"
    js = tmp_path / (os.path.basename(page) + ".js")
    js.write_text("\n".join(scripts))
    result = subprocess.run([NODE, "--check", str(js)], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr


def test_the_dashboard_script_parses(tmp_path):
    db = str(tmp_path / "t.db")
    make_db(db)
    nt = NewickTree(db, name="t", outdir=str(tmp_path))
    check_scripts(nt.write_dashboard(REPORTS, str(tmp_path / "summary.html")), tmp_path)


@pytest.mark.parametrize("summary", [False, True])
def test_the_tree_page_script_parses(tmp_path, summary):
    db = str(tmp_path / "t.db")
    make_db(db)
    nt = NewickTree(db, name="t", outdir=str(tmp_path))
    nt.draw_tree({"A.1": 1}, ["A.1"], True, summary=summary)
    check_scripts(nt.tree_file, tmp_path)
