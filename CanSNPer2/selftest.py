#!/usr/bin/env python3
'''
CanSNPer self test for conda
Check so that all required packages are installed
'''

import shutil
import sys

def main():
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

	if not shutil.which("progressiveMauve"):
		failed.append("progressiveMauve was not found in PATH, it is required to align the genomes")

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
