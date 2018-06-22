#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import os
import shutil
import numpy as np
from numpy import genfromtxt
from multiprocessing import Process
from multiprocessing import Pool
from multiprocessing import cpu_count

if __name__ == "__main__":

	# command line arguments

	parser = argparse.ArgumentParser(description = "batch generating hilbert map")
	parser.add_argument("-i", "--input", dest = "input", required = True,
		help = "input data folder", metavar = "FOLDER")
	parser.add_argument("-m", "--meta", dest = "meta", required = True,
		help = "cleaned meta filename", metavar = "FILE")
	parser.add_argument("-o", "--output", dest = "output", required = False,
		help = "output folder", metavar = "FOLDER")
	parser.add_argument("-p", "--process", dest = "process", required = False,
		help = "number of parallel processes", metavar = "INTEGER")

	args = parser.parse_args()
	folder_input = args.input
	filename_meta = args.meta
	folder_output = args.output
	n_p = args.process

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

	if n_p == None :
		n_p = cpu_count()
		print('[!] missing parallel parameters, using CPU count {0}'.format(n_p))
		if n_p <= 0 :
			n_p = 1
			print('[!] undetected CPU count, using 1')

	# read & parse meta file

	print('[*] parsing meta file ...')
	meta_data = genfromtxt(filename_meta, delimiter = ',',
		autostrip = True,
		dtype = ('|S64', '|S128', '|S32', '|S64', 'b1', '<f8', '|S16', 'b1', '<i4', '|S32', 'b1'),
		names = ('file_id', 'file_name', 'primary_site', 'disease_type', 'gender', 'age', 'race', 'is_tumor', 'sample_type_id', 'sample_type', 'is_alive'))
	
	# processing meta data

	for i in xrange(0, len(meta_data), n_p):
        p_data =  meta_data[i:i + n_p]
        print("i:\n")
        print(p_data)

	print('[*] complete')