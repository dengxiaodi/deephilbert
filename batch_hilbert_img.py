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

def load_meta_data(filename_meta) :
	meta_data = genfromtxt(filename_meta, delimiter = ',',
		autostrip = True,
		dtype = ('|S64', '|S128', '|S32', '|S64', 'b1', '<f8', '|S16', 'b1', '<i4', '|S32', 'b1'),
		names = ('file_id', 'file_name', 'primary_site', 'disease_type', 'gender', 'age', 'race', 'is_tumor', 'sample_type_id', 'sample_type', 'is_alive'))
	return meta_data

def load_cg_meth(filename_cg_meth) :
	cg_meth = np.load(filename_cg_meth)
	return cg_meth['meth']

def generate_hilbert_map(v, r, N, filename) :
	plt_data = np.array(hb.hilbert_plot(v, r))
	x = plt_data[:, 0]
	y = plt_data[:, 1]
	z = plt_data[:, 2]

	dmax = 2**r
	xs, ys = np.mgrid[0:dmax:N, 0:dmax:N]
	zs = griddata((x, y), z, (xs, ys))
	plt.imsave(filename, zs.T, cmap = cm.gray_r, vmin = 0, vmax = 1)

def process_hilbert(meta, r, N, folder_input, folder_output, total_count, current_count) :
	proc_name = current_process().name
	print('[*] processing meta {0} on process {1} ({2}/{3})'.format(meta['file_id'], proc_name, current_count, total_count))

	# mkdir & clean meth data
		
	folder_meta = os.path.join(folder_output, meta['file_id'])
	try:
		os.stat(folder_meta)
	except:
		os.mkdir(folder_meta)

	# load meth data

	filename_cg_meth = os.path.join(folder_meta, os.path.splitext(meta['file_name'])[0] + '.cg.meth.npz')
	if not os.path.isfile(filename_cg_meth) :
		print('[!] methylation data file "{0}" does not exist, skip'.format(filename_cg_meth))
		pass

	cg_meth_data = load_cg_meth(filename_cg_meth)

	# generate hilbert map

	filename_hilbert = os.path.join(folder_meta, os.path.splitext(meta['file_name'])[0] + '.png')
	generate_hilbert_map(cg_meth_data, r, N, filename_hilbert)
	print('[*] meta {0} on process {1} complete (#{2})'.format(meta['file_id'], proc_name, current_count))
	print('[*] {0} written'.format(filename_hilbert))

if __name__ == "__main__":

	# command line arguments

	parser = argparse.ArgumentParser(description = "batch generating hilbert image")
	parser.add_argument("-i", "--input", dest = "input", required = True,
		help = "input data folder", metavar = "FOLDER")
	parser.add_argument("-m", "--meta", dest = "meta", required = True,
		help = "cleaned meta filename", metavar = "FILE")
	parser.add_argument("-o", "--output", dest = "output", required = False,
		help = "output folder", metavar = "FOLDER")
	parser.add_argument("-p", "--process", dest = "process", required = False, type = int,
		help = "number of parallel processes", metavar = "INTEGER")
	parser.add_argument("-r", "--resolution", dest = "resolution", required = False, type = int, default = 10,
		help = "hilbert resolution", metavar = "INTEGER")
	parser.add_argument("-N", "--image-point", dest = "N", required = False, type = int, default = 300,
		help = "image point in imaginary form", metavar = "INTEGER")

	args = parser.parse_args()
	folder_input = args.input
	filename_meta = args.meta
	folder_output = args.output
	n_p = args.process
	image_size = complex(args.N)
	resolution = args.resolution

	if not os.path.isdir(folder_input) :
		print('[~] "{0}" input does not exist!'.format(folder_input))
		exit(-1)

	if not os.path.isfile(filename_meta) :
		print('[~] meta file "{0}" does not exist!'.format(filename_meta))
		exit(-1)

	if folder_output == None :
		folder_output = os.path.basename(os.path.normpath(folder_input)) + '_hilbert'

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

	# processing meta data

	p = Pool(n_p)
	total_count = len(meta_data)
	current_count = 0
	for meta in meta_data:
		current_count += 1
		p.apply_async(process_hilbert, (meta, resolution, image_size, folder_input, folder_output, total_count, current_count, ))

	p.close()
	p.join()

	print('[*] complete')