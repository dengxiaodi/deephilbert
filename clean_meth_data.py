#!/usr/bin/env python
# -*- coding: utf-8 -*-

import csv
import argparse
import os
import shutil
import numpy as np
from numpy import genfromtxt

def clean_meth_data(filename_input, filename_output) :
	methy_data = genfromtxt(filename_input, delimiter = '\t', skip_header = 1,
		usecols = (0, 1, 2, 3),
		dtype = ('S12', '<f8', 'S8', '<u8'), 
		names = ('ref', 'beta', 'chr', 'pos'))

	# filter out all unknown chrs

	methy_data = methy_data[methy_data['chr'] != '*']

	# filter out all nan beta values

	methy_data = methy_data[~np.isnan(methy_data['beta'])]
	if filename_output == None :
		filename_output = os.path.splitext(os.path.basename(filename_input))[0] + '.clean.csv'
	np.savetxt(filename_output, 
		methy_data, 
		fmt=('%s, %f, %s, %ld'), 
		delimiter = ',')

	return filename_output

if __name__ == "__main__":

	# command line arguments

	parser = argparse.ArgumentParser(description = "clean methylation data")
	parser.add_argument("-i", "--input", dest = "input", required = True,
		help = "methylation.txt with betavalue downloaded from TCGA data repository", metavar = "FILE")
	parser.add_argument("-o", "--output", dest = "output", required = False,
		help = "output filename", metavar = "FILE")

	args = parser.parse_args()
	filename_input = args.input
	filename_output = args.output

	# parameter check

	if not os.path.isfile(filename_input) :
		print("[~] TCGA methylation.text file \"{}\" does not exist!".format(filename_input))
		exit(-1)

	# read & parse input file

	print('[*] parsing methylation data file ...')
	filename_output = clean_meth_data(filename_input, filename_output)
	print('[*] {0} written complete'.format(filename_output))