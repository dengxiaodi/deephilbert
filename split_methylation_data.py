#!/usr/bin/env python
# -*- coding: utf-8 -*-

import csv
import argparse
import os
import shutil
import numpy as np
from numpy import genfromtxt

def split_methylation_file(filename, output_folder) :
	methy_data = genfromtxt(filename, delimiter = '\t', skip_header = 1,
		usecols = (0, 1, 2, 3),
		dtype = ('S12', '<f8', 'S8', '<u8'), 
		names = ('ref', 'beta', 'chr', 'pos'))

	# filter out all unknown chrs

	methy_data = methy_data[methy_data['chr'] != '*']

	# filter out all nan beta values

	methy_data = methy_data[~np.isnan(methy_data['beta'])]

	chrs = np.unique(methy_data['chr'])
	for chr in chrs:
		np.savetxt(output_folder + '/' + chr + '.csv', 
			methy_data[methy_data['chr'] == chr], 
			fmt=('%s, %f, %s, %ld'),
			delimiter = ',')
		print('** {chr}.csv written complete'.format(chr = chr))

if __name__ == "__main__":

	# command line arguments

	parser = argparse.ArgumentParser(description = "split methylation data by chr")
	parser.add_argument("-i", "--input", dest = "input", required = True,
		help = "methylation.txt with betavalue downloaded from TCGA data repository", metavar = "FILE")
	parser.add_argument("-o", "--output", dest = "output", required = True,
		help = "output foldername", metavar = "FOLDER")

	args = parser.parse_args()
	filename_input = args.input
	folder_output = args.output

	# parameter check

	if not os.path.isfile(filename_input) :
		print("TCGA methylation.text file \"" + filename_input + "\" does not exist!")
		exit(-1)

	if os.path.isdir(folder_output) :
		shutil.rmtree(folder_output, ignore_errors = True)
	os.makedirs(folder_output)

	# read & parse input file

	print('[*] parsing methylation data file ...')
	split_methylation_file(filename_input, folder_output)
	print('[*] splite complete')