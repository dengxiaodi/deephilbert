#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import os
import shutil
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

import hilbert as hb

def generate_hilbert_map(v, r) :
	plt_data = np.array(hb.hilbert_plot(v, r))
	x = plt_data[:, 0]
	y = plt_data[:, 1]
	z = plt_data[:, 2]
	N = 300j
	xmin = int(np.floor(np.min(x)))
	xmax = int(np.ceil(np.max(x)))
	ymin = int(np.floor(np.min(y)))
	ymax = int(np.ceil(np.max(y)))
	extent = (xmin, xmax, ymin, ymax)
	xs, ys = np.mgrid[xmin:xmax:N, ymin:ymax:N]
	zs = griddata((x, y), z, (xs, ys))

	plt.imshow(rdata.T, extent = extent)
	plt.savefig("test.png")

if __name__ == "__main__":

	# command line arguments

	parser = argparse.ArgumentParser(description = "Generate CG hilbert maps")
	parser.add_argument("-i", "--input", dest = "input", required = True,
		help = "methylation.txt with betavalue downloaded from TCGA data repository", metavar = "FILE")
	parser.add_argument("-d", "--index", dest = "index", required = True,
		help = "cg index cgidx.npz file", metavar = "FILE")
	parser.add_argument("-o", "--output", dest = "output", required = True,
		help = "output foldername", metavar = "FOLDER")
	parser.add_argument("-c", "--chr", dest = "chr", required = False,
		help = "output only chr", metavar = "CHROMESOME")

	args = parser.parse_args()
	filename_input = args.input
	filename_index = args.index
	folder_output = args.output

	# parameter check

	if not os.path.isfile(filename_input) :
		print("TCGA methylation.text file \"" + filename_input + "\" does not exist!")
		exit(-1)

	if not os.path.isfile(filename_index) :
		print(".cgidx.npz file \"" + filename_index + "\" does not exist!")
		exit(-1)

	if os.path.isdir(folder_output) :
		shutil.rmtree(folder_output, ignore_errors = True)
	os.makedirs(folder_output)

	# load indexes

	cg_index = np.load(filename_index)
	chrsize = cg_index['chrsize']
	cg_indexes = cg_index['indexes']

	print('* load cg indexes complete')

	# load methylation data

	methy_data = np.genfromtxt(filename_input, delimiter = ',', 
		dtype = ('S12', '<f8', 'S8', '<u8'),
		names = ('ref', 'beta', 'chr', 'pos'))
	methy_genome = np.zeros(chrsize)
	methy_genome[methy_data['pos']] = methy_data['beta']
	methy_cg = methy_genome[cg_indexes]

	# generate hilbert curve

	print(methy_genome)
	print(methy_cg)

	v = np.random.normal(0, 1, size = 5000)

	generate_hilbert_map(v, 5)

	print('* complete')