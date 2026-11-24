"""Since 2.6.0 one command is installed: cansnper2, and its help names the verbs.

The four separate commands made the same jobs look like four tools with four
spellings of the same options. They are no longer installed; the old modules keep
their code for the transition (the test suite and the release gate still drive
them) and their --help points at the verbs."""
import os
import re
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
VERBS = ["type", "selftest", "db", "scheme", "ref"]
OLD_MODULES = ["CanSNPer2.CanSNPerTree", "CanSNPer2.SNPDatabase",
	"CanSNPer2.DownloadGenomes", "CanSNPer2.selftest"]


def help_output(*args):
	env = dict(os.environ, PYTHONPATH=ROOT)
	run = subprocess.run([sys.executable, "-m", "CanSNPer2.cli"] + list(args) + ["--help"],
		capture_output=True, text=True, cwd=ROOT, env=env)
	assert run.returncode == 0, run.stderr
	return run.stdout


def test_only_cansnper2_is_installed():
	'''the console scripts of setup.py: the old commands were replaced in 2.6.0'''
	with open(os.path.join(ROOT, "setup.py")) as f:
		entry_points = re.search(r"console_scripts['\"]?\s*:\s*\[(.*?)\]", f.read(), re.S).group(1)
	scripts = re.findall(r"['\"]([A-Za-z0-9_-]+)=", entry_points)
	assert scripts == ["cansnper2"], "setup.py installs {0}".format(scripts)


def test_the_help_names_every_verb():
	out = help_output()
	for verb in VERBS:
		assert verb in out, "cansnper2 --help does not name {0}".format(verb)


def test_the_old_modules_point_at_the_verbs():
	'''their code still runs during the transition, and their help has to send
		whoever finds it to cansnper2'''
	env = dict(os.environ, PYTHONPATH=ROOT)
	for module in OLD_MODULES:
		run = subprocess.run([sys.executable, "-m", module, "--help"],
			capture_output=True, text=True, cwd=ROOT, env=env)
		assert run.returncode == 0, run.stderr
		out = run.stdout + run.stderr
		for line in COMMAND_MENTIONS:
			assert line in out, "{module} --help does not mention {line}".format(
				module=module, line=line)


COMMAND_MENTIONS = ["cansnper2 type", "cansnper2 db", "cansnper2 scheme", "cansnper2 ref"]
