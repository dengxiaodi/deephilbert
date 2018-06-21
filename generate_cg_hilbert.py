#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import os
import shutil
import numpy as np
import matplotlib as mpl
if os.environ.get('DISPLAY','') == '':
    print('[!] No display found, using non-interactive Agg backend')
    mpl.use('Agg')
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import matplotlib.cm as cm

import hilbert as hb

def generate_hilbert_map(v, r) :
	plt_data = np.array(hb.hilbert_plot(v, r))
	x = plt_data[:, 0]
	y = plt_data[:, 1]
	z = plt_data[:, 2]
	N = 300j

	print(max(x))
	print(max(y))
	dmax = 2**r
	xs, ys = np.mgrid[0:dmax:N, 0:dmax:N]
	zs = griddata((x, y), z, (xs, ys))
	plt.imsave("t.png", zs.T, cmap = cm.gray_r, vmin = 0, vmax = 1)

if __name__ == "__main__":

	# command line arguments

	parser = argparse.ArgumentParser(description = "Generate CG hilbert maps")
	parser.add_argument("-i", "--input", dest = "input", required = True,
		help = "methylation.txt with betavalue downloaded from TCGA data repository", metavar = "FILE")
	parser.add_argument("-d", "--index", dest = "index", required = True,
		help = "cg index cgidx.npz file", metavar = "FILE")
	parser.add_argument("-o", "--output", dest = "output", required = True,
		help = "output foldername", metavar = "FILE")

	args = parser.parse_args()
	filename_input = args.input
	filename_index = args.index
	filename_output = args.output

	# parameter check

	if not os.path.isfile(filename_input) :
		print("[~] cleaned methylation file \"{0}\" does not exist!".format(filename_input))
		exit(-1)

	if not os.path.isfile(filename_index) :
		print("[~] CG index .cgidx.npz file \"{0}\" does not exist!".format(filename_index))
		exit(-1)

	if os.path.isfile(filename_output) :
		os.unlink(filename_output)
	
	# load indexes

	print('[*] loading cg indexes')
	cg_indexes_file = np.load(filename_index)
	cg_indexes = cg_indexes_file['indexes']
	dict_cg_indexes = cg_indexes.item()

	# load methylation data

	print('[*] loading methyaltion data')
	meth_data = np.genfromtxt(filename_input, delimiter = ',', 
		dtype = ('S12', '<f8', 'S8', '<u8'),
		names = ('ref', 'beta', 'chr', 'pos'))

	print(meth_data)

	print('[*] building methylation data list')
	cg_meth = np.array([], dtype = '<f8')
	genome_meth = np.array([], dtype = '<f8')
	for chrname in dict_cg_indexes :
		print('    parsing: {0}'.format(chrname))
		chrdata = dict_cg_indexes[chrname]
		chr_size = chrdata[1]
		chr_cg_index = chrdata[2]

		meth = np.zeros(chr_size)
		chr_meth_data = meth_data[meth_data['chr'] == chrname]
		meth[chr_meth_data['pos'] - 1] = chr_meth_data['beta']
		cg = meth[chr_cg_index - 1]

		print(np.sort(chr_meth_data['pos'] - 1))
		print(np.sort(chr_cg_index - 1))

		cg_meth = np.append(cg_meth, cg)
		genome_meth = np.append(genome_meth, meth)

	print('meth_data size: {0}'.format(len(meth_data)))
	print('hitted meth_data size: {0}'.format(len(cg_meth[cg_meth != 0])))

	# generate hilbert curve

	# generate_hilbert_map(cg, 10)

	print('[*] complete')