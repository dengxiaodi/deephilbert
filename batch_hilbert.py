#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import os
import sys
import shutil
import numpy as np
from numpy import genfromtxt
from multiprocessing import Process
from multiprocessing import Pool
from multiprocessing import cpu_count, current_process

script_path = os.path.abspath(os.path.dirname(sys.argv[0]))
clean_meth_script = os.path.join(script_path, 'clean_meth_data.py')
generate_cg_hilbert_script = os.path.join(script_path, 'generate_cg_hilbert.py')

def batch_hilbert(meta_chunk, filename_index, folder_input, folder_output) :
	proc_name = current_process().name
	print('[*] processing {0} records on process {1}'.format(len(meta_chunk), proc_name))
	for meta in meta_chunk :
		
		# mkdir & clean meth data
		
		folder_meta = os.path.join(folder_output, meta['file_id'])
		try:
			os.stat(folder_meta)
		except:
			os.mkdir(folder_meta)

		# read input file

		filename_meth = os.path.join(folder_input, meta['file_id'], meta['file_name'])
		if not os.path.isfile(filename_meth) :
			print('[!] methylation file "{0}" does not exist, skip'.format(filename_meth))
			pass
		filename_meth_clean = os.path.join(folder_meta, os.path.splitext(meta['file_name'])[0] + '.clean.csv')
		os.system(clean_meth_script + ' -i ' + filename_meth + ' -o ' + filename_meth_clean)

		# generate hilbert map for cleaned meth data

		filename_hilbert_map = os.path.join(folder_meta, os.path.splitext(meta['file_name'])[0])
		os.system(generate_cg_hilbert_script + ' -i ' + filename_meth_clean + ' -d '  + filename_index + ' -o ' + filename_hilbert_map)


if __name__ == "__main__":

	# command line arguments

	parser = argparse.ArgumentParser(description = "batch generating hilbert map")
	parser.add_argument("-i", "--input", dest = "input", required = True,
		help = "input data folder", metavar = "FOLDER")
	parser.add_argument("-d", "--index", dest = "index", required = True,
		help = "cg index cgidx.npz file", metavar = "FILE")
	parser.add_argument("-m", "--meta", dest = "meta", required = True,
		help = "cleaned meta filename", metavar = "FILE")
	parser.add_argument("-o", "--output", dest = "output", required = False,
		help = "output folder", metavar = "FOLDER")
	parser.add_argument("-p", "--process", dest = "process", required = False,
		help = "number of parallel processes", metavar = "INTEGER")

	args = parser.parse_args()
	folder_input = args.input
	filename_meta = args.meta
	filename_index = args.index
	folder_output = args.output
	n_p = args.process

	if not os.path.isdir(folder_input) :
		print('[~] "{0}" input does not exist!'.format(folder_input))
		exit(-1)

	if not os.path.isfile(filename_meta) :
		print('[~] meta file "{0}" does not exist!'.format(filename_meta))
		exit(-1)

	if not os.path.isfile(filename_index) :
		print("[~] CG index .cgidx.npz file \"{0}\" does not exist!".format(filename_index))
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

	p = Pool(n_p)
	meta_data_chunks = np.array_split(meta_data, n_p)
	for meta_chunk in meta_data_chunks:
		p.apply_async(batch_hilbert, (meta_chunk, filename_index, folder_input, folder_output, ))

	p.close()
	p.join()

	print('[*] complete')