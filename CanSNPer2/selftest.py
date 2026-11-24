#!/usr/bin/env python3
'''
CanSNPer self test for conda
Check that CanSNPer2 can be imported and which aligner programs are installed
'''

import shutil
import sys
import argparse

from CanSNPer2.modules.TopLevel import run_main
from CanSNPer2 import COMMANDS_HELP

def run_command():
	'''the check itself, without a parser: cansnper2 selftest calls this directly (2.6.0)'''
	'''CanSNPer2 needs no third party python packages any more: the ete3/PyQt tree drawing
		was replaced by the html viewer in modules/NewickTree.py, so ete3 could be dropped
		(it imports cgi, removed in python 3.13, which capped the package below 3.13, issue #184).
		Import the package itself to catch syntax/import breakage of the installed copy.'''
	try:
		import CanSNPer2.modules.NewickTree
	except ImportError as e:
		print("ERROR: CanSNPer2 could not be imported: {e}; reinstall CanSNPer2".format(e=e))
		sys.exit(1)

	'''each aligner is reported; at least one has to be usable (found, with delta-filter for nucmer)'''
	from CanSNPer2.modules.Aligners import ALIGNERS, DEFAULT_ALIGNER
	usable = []
	for name, cls in sorted(ALIGNERS.items()):
		aligner = cls()
		found = shutil.which(aligner.program)
		print("{program}: {state} (--aligner {name}{default})".format(program=aligner.program,
			state=found or "not found", name=name, default=", default" if name == DEFAULT_ALIGNER else ""))
		ready = found is not None
		if hasattr(aligner, "filter_program"):
			filter_found = shutil.which(aligner.filter_program())
			print("{program}: {state} (needed by {name})".format(program=aligner.filter_program(),
				state=filter_found or "not found", name=aligner.program))
			ready = ready and filter_found is not None
		if ready:
			usable.append(name)

	if not usable:
		print("ERROR: no aligner was found in PATH; install progressiveMauve, minimap2 or MUMmer4 (nucmer and delta-filter)")
		sys.exit(1)
	print("usable aligners: {names}".format(names=", ".join(usable)))


def command():
	'''only --debug today, but the parser makes --help work and show the command list'''
	parser = argparse.ArgumentParser(
		description='CanSNPer2-test - check the installation: python imports and aligner programs in PATH',
		epilog=COMMANDS_HELP,
		formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument('--debug', action='store_true', help="print the traceback of an unexpected error")
	parser.parse_args()
	run_command()


def main():
	'''the entry point of the command: command() with the handler of unexpected errors'''
	return run_main(command)


if __name__ == "__main__":
	'''without this guard "python -m CanSNPer2.selftest" did nothing at all, the same missing guard
		as issue #128
	'''
	main()
