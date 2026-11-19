"""Every CanSNPer2 command must name the whole command family in its --help.

CanSNPer2-database and CanSNPer2-download were invisible to users who only ran
CanSNPer2 --help: the sibling commands were mentioned nowhere but the README,
so nothing pointed a new user at database creation or reference download."""
import subprocess
import sys

COMMANDS = ["CanSNPer2", "CanSNPer2-database", "CanSNPer2-download", "CanSNPer2-test"]
MODULES = ["CanSNPer2.CanSNPerTree", "CanSNPer2.SNPDatabase",
	"CanSNPer2.DownloadGenomes", "CanSNPer2.selftest"]


def help_output(module):
	run = subprocess.run([sys.executable, "-m", module, "--help"],
		capture_output=True, text=True)
	assert run.returncode == 0, run.stderr
	return run.stdout


def test_every_help_lists_every_command():
	for module in MODULES:
		out = help_output(module)
		for command in COMMANDS:
			assert command in out, "{module} --help does not mention {command}".format(
				module=module, command=command)
