#!/usr/bin/env python3 -c
import sys, os
from subprocess import Popen,PIPE,STDOUT

import logging
logger = logging.getLogger(__name__)

class AlignmentError(Exception):
	"""docstring for MauveE"""
	def __init__(self, value):
		self.value = value
	def __str__(self):
		return repr(self.value)

class RunAlignment(object):
	"""docstring for RunAlignment"""
	def __init__(self, database,organism,query **kwargs):
		super(RunAlignment, self).__init__()

	def get_references(self,selected=False):
		'''references must be available in the refdir'''
		return [ref for ref in os.listdir(self.refdir) if ref.endswith(".fa")]

	def run(self,database,organism):
		'''Run RunAlignment'''
		print("Run alignments to references using progressiveMauve")

		'''Read query input file if a txt file is supplied insead of fasta files'''
		if self.query[0].endswith(".txt"):
			self.query=self.read_query_input(self.query)

		'''Main function of RunAlignment
				1. Align reads with BWA/Minimap2
				2. Call SNPs
				3. Return results
		'''
