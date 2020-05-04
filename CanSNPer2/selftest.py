#!/usr/bin/env python3
'''
CanSNPer self test for conda
Check so that all required packages are installed
'''

try:
	import ete3
except ImportError:
	raise ImportError("ete3 could not be found!")

try:
	import flextaxd
except ImportError:
	raise ImportError("flextaxd could not be found!")


def main():
	print("required python packages installed!")
