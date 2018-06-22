#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import os
import shutil
import numpy as np
from numpy import genfromtxt

if __name__ == "__main__":

	# command line arguments

	parser = argparse.ArgumentParser(description = "batch generating hilbert map")
	parser.add_argument("-i", "--input", dest = "input", required = True,
		help = "input data folder", metavar = "FOLDER")
	parser.add_argument("-m", "--meta", dest = "meta", required = True,
		help = "cleaned meta filename", metavar = "FILE")
	parser.add_argument("-o", "--output", dest = "output", required = False,
		help = "output folder", metavar = "FOLDER")

	args = parser.parse_args()
	folder_input = args.input
	filename_meta = args.meta
	folder_output = args.output

	if not os.path.isdir(folder_input) :
		print('[~] "{0}" input does not exist!'.format(folder_input))
		exit(-1)

	if not os.path.isfile(filename_meta) :
		print('[~] meta file "{0}" does not exist!'.format(filename_meta))
		exit(-1)

	if folder_output == None :
		folder_output = os.path.basename(os.path.normpath(folder_input)) + '_hb'

	if os.path.isdir(folder_output) :
		shutil.rmtree(folder_output, ignore_errors = True)
	os.makedirs(folder_output)

	# read & parse meta file

	print('[*] parsing meta file ...')
	meta_data = genfromtxt(filename_meta, delimiter = ',',
		autostrip = True,
		dtype = ('|S64', '|S128', '|S32', '|S64', 'b', '<f8', '|S16', 'b', '<i4', '|S32', 'b'),
		names = ('file_id', 'file_name', 'primary_site', 'disease_type', 'gender', 'age', 'race', 'is_tumor', 'sample_type_id', 'sample_type', 'is_alive'))
	
	# processing meta data

	for meta in meta_data :
		print(meta)

	print('[*] complete')