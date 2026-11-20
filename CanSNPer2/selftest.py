#!/usr/bin/env python3
'''
CanSNPer self test for conda
Check so that all required packages are installed
'''

import shutil
import sys
import argparse

from CanSNPer2 import COMMANDS_HELP

def main():
	'''no options today, but the parser makes --help work and show the command list'''
	parser = argparse.ArgumentParser(
		description='CanSNPer2-test - check the installation: python imports and aligner programs in PATH',
		epilog=COMMANDS_HELP,
		formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.parse_args()

	failed = []
	'''CanSNPer2 needs no third party python packages any more: the ete3/PyQt tree drawing
		was replaced by the html viewer in modules/NewickTree.py, so ete3 could be dropped
		(it imports cgi, removed in python 3.13, which capped the package below 3.13, issue #184).
		Import the package itself to catch syntax/import breakage of the installed copy.'''
	try:
		import CanSNPer2.modules.NewickTree
	except ImportError as e:
		print("ERROR: CanSNPer2 could not be imported: {e}; install the requirements: pip install -r requirements.txt".format(e=e))
		sys.exit(1)

	'''the default aligner must be installed; the others are optional (--aligner)'''
	from CanSNPer2.modules.Aligners import ALIGNERS, DEFAULT_ALIGNER
	for name, aligner in sorted(ALIGNERS.items()):
		found = shutil.which(aligner.program)
		if name == DEFAULT_ALIGNER and not found:
			failed.append("{program} was not found in PATH, it is required to align the genomes".format(program=aligner.program))
		elif name != DEFAULT_ALIGNER:
			print("{program}: {state} (optional, --aligner {name})".format(program=aligner.program,state=found or "not found",name=name))

	if failed:
		for error in failed:
			print("ERROR: "+error)
		sys.exit(1)
	print("required python packages installed!")

if __name__ == "__main__":
	'''without this guard "python -m CanSNPer2.selftest" did nothing at all, the same missing guard
		as issue #128
	'''
	main()
