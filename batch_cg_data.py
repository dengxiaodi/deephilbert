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

def load_cg_indexes(filename_index) :
	cg_indexes_file = np.load(filename_index)
	cg_indexes = cg_indexes_file['indexes']
	dict_cg_indexes = cg_indexes.item()

	return dict_cg_indexes

def load_meta_data(filename_meta) :
	meta_data = genfromtxt(filename_meta, delimiter = ',',
		autostrip = True,
		dtype = ('|S64', '|S128', '|S32', '|S64', 'b1', '<f8', '|S16', 'b1', '<i4', '|S32', 'b1'),
		names = ('file_id', 'file_name', 'primary_site', 'disease_type', 'gender', 'age', 'race', 'is_tumor', 'sample_type_id', 'sample_type', 'is_alive'))
	return meta_data

def load_meth_data(filename_meth) :
	meth_data = genfromtxt(filename_meth, delimiter = '\t', skip_header = 1,
		usecols = (0, 1, 2, 3), autostrip = True,
		dtype = ('|S12', '<f8', '|S8', '<u8'), 
		names = ('ref', 'beta', 'chr', 'pos'))

	# filter out all unknown chrs

	meth_data = meth_data[meth_data['chr'] != '*']

	# filter out all nan beta values

	meth_data = meth_data[~np.isnan(meth_data['beta'])]

	return meth_data

def write_meth_data(meth_data, filename_output) :
	np.savetxt(filename_output, 
		meth_data, 
		fmt=('%s, %f, %s, %ld'), 
		delimiter = ',')

def generate_cg_data(meth_data, dict_cg_indexes, filename_output) :
	cg_meth = np.array([], dtype = '<f8')
	genome_meth = np.array([], dtype = '<f8')
	for chrname in dict_cg_indexes :
		chrdata = dict_cg_indexes[chrname]
		chr_size = chrdata[1]
		chr_cg_index = chrdata[2]

		meth = np.zeros(chr_size)
		chr_meth_data = meth_data[meth_data['chr'] == chrname]
		meth[chr_meth_data['pos'] - 1] = chr_meth_data['beta']
		cg = meth[chr_cg_index - 1]

		cg_meth = np.append(cg_meth, cg)
		genome_meth = np.append(genome_meth, meth)
	
	filename_cg_meth = os.path.splitext(os.path.basename(filename_output))[0] + '.cg.meth'
	filename_genome_meth = os.path.splitext(os.path.basename(filename_output))[0] + '.genome.meth'
	np.savez(filename_cg_meth, meth = cg_meth)
	np.savez(filename_genome_meth, meth = genome_meth)

def process_meta_data(meta, dict_cg_indexes, folder_input, folder_output, total_count, current_count) :
	proc_name = current_process().name
	print('[*] processing meta {0} on process {1} [{2} / {3}]'.format(meta['file_id'], proc_name, current_count, total_count))

	# mkdir & clean meth data
		
	folder_meta = os.path.join(folder_output, meta['file_id'])
	try:
		os.stat(folder_meta)
	except:
		os.mkdir(folder_meta)

	# load meth data

	filename_meth = os.path.join(folder_input, meta['file_id'], meta['file_name'])
	if not os.path.isfile(filename_meth) :
		print('[!] methylation file "{0}" does not exist, skip'.format(filename_meth))
		pass
	filename_meth_clean = os.path.join(folder_meta, os.path.splitext(meta['file_name'])[0] + '.clean.csv')

	meth_data = load_meth_data(filename_meth)
	write_meth_data(meth_data, filename_meth_clean)

	# generate cg data

	filename_cg_meth = os.path.join(folder_meta, os.path.splitext(meta['file_name'])[0])
	generate_cg_data(meth_data, dict_cg_indexes, filename_cg_meth)

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
	parser.add_argument("-p", "--process", dest = "process", required = False, type = int,
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

	print('[*] parsing meta file')
	meta_data = load_meta_data(filename_meta)

	print('[*] loading cg indexes')
	dict_cg_indexes = load_cg_indexes(filename_index)

	# processing meta data

	p = Pool(n_p)
	total_count = len(meta_data)
	current_count = 0
	for meta in meta_data:
		current_count += 1
		p.apply_async(process_meta_data, (meta, dict_cg_indexes, folder_input, folder_output, total_count, current_count, ))

	p.close()
	p.join()

	print('[*] complete')