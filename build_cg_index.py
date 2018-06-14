#!/usr/bin/env python
# -*- coding: utf-8 -*-

def build_cg_index(filename_input, filename_output) :
	# load sequence file

	chrname = ''
	seq = ''
	with open(filename_input, 'r') as seq_file :
		
		for line in seq_file:
			if(line[0] == '>'):
				chrname = line[1:].strip()
			else:
				seq += line.strip().upper()

	# build cg indexes
	
	indexes = [m.start() for m in re.finditer("CG", seq)]
	indexes = np.array(indexes)
	if filename_output == None :
		filename_output = chrname + '.cgidx'
	np.savez(filename_output, chrsize = len(seq), indexes = indexes)

	return filename_output

if __name__ == "__main__":
	import argparse
	import os
	import numpy as np
	import re

	# command line arguments

	parser = argparse.ArgumentParser(description = "Build CpG indexes from reference genomes")
	parser.add_argument("-i", "--input", dest = "input", required = True,
		help = "reference genome file by chromsome chr.fa", metavar = "FILE")
	parser.add_argument("-o", "--output", dest = "output", required = False,
		help = "output cg index file name chr.cgidx", metavar = "FILE")

	args = parser.parse_args()
	filename_input = args.input
	filename_output = args.output

	# parameter check

	if not os.path.isfile(filename_input) :
		print("reference genome .fa file \"" + filename_input + "\" does not exist!")
		exit(-1)

	if filename_output != None and os.path.isfile(filename_output) :
		os.unlink(filename_output)

	filename_output = build_cg_index(filename_input, filename_output)
	print('cg index file "' + filename_output + '" write complete.')