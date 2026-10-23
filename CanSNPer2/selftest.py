#!/usr/bin/env python3
'''
CanSNPer self test for conda
Check so that all required packages are installed
'''

import shutil
import sys

def main():
	failed = []
	try:
		import ete3
	except ImportError:
		raise ImportError("ete3 could not be found!")

	'''ete3 imports its drawing classes from the treeview subpackage, which silently fails
		when no Qt binding is installed, and CanSNPer2 cannot draw any tree without them'''
	for name in ("NodeStyle","TreeStyle","AttrFace","faces"):
		if not hasattr(ete3,name):
			failed.append("ete3 is missing the drawing class {name}, install PyQt5 (pip install PyQt5)".format(name=name))

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
