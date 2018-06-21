#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import os
import numpy as np
import re

def build_cg_index(filename_ref_seq, filename_cg_index) :
	# load sequence file
	print('[*] loading reference sequence')
	list_ref_seq = []
	with open(filename_ref_seq, 'r') as ref_seq_file :
		chrname = ''
		seq = ''
		for line in ref_seq_file:
			if(line[0] == '>'):
				
				# save current seq for current chr

				if(chrname != ''):
					list_ref_seq += [(chrname, seq)]
				
				# new chrname & seq

				chrname = line[1:].strip()
				seq = ''
				print('    loading reference sequence: {0}'.format(chrname))
			else:
				seq += line.strip().upper()
		
		# write the last chr

		if(chrname != ''):
			list_ref_seq += [(chrname, seq)]
	ref_seq_file.close()

	# build cg indexes

	print('[*] building cg indexes')
	list_cg_indexes = []
	offset = 0
	for ref_seq in list_ref_seq:
		chrname = ref_seq[0]
		seq = ref_seq[1]
		print('    building cg indexes for reference sequence: {0}'.format(chrname))
		chrlen = len(seq)
		offset += chrlen
		indexes = [m.start() for m in re.finditer("CG", seq)]
		list_cg_indexes += [(offset, chrlen, np.array(indexes))]
	
	if filename_cg_index == None :
		filename_cg_index = os.path.splitext(os.path.basename(filename_ref_seq))[0] + '.cgidx'
	np.savez(filename_cg_index, list_cg_indexes)

	return filename_cg_index

if __name__ == "__main__":

	# command line arguments

	parser = argparse.ArgumentParser(description = "Build CpG indexes from reference genomes")
	parser.add_argument("-i", "--input", dest = "input", required = True,
		help = "reference genome file", metavar = "FILE")
	parser.add_argument("-o", "--output", dest = "output", required = False,
		help = "output cg index file name .cgidx", metavar = "FILE")

	args = parser.parse_args()
	filename_input = args.input
	filename_output = args.output

	# parameter check

	if not os.path.isfile(filename_input) :
		print("[~] reference genome .fa file \"" + filename_input + "\" does not exist!")
		exit(-1)

	if filename_output != None and os.path.isfile(filename_output) :
		os.unlink(filename_output)

	filename_output = build_cg_index(filename_input, filename_output)
	print('[*] cg index file "' + filename_output + '" write complete.')