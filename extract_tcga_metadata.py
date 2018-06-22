#!/usr/bin/env python
# -*- coding: utf-8 -*-

def extract_tcga_metadata(metadata, output_file) :
	total_count = 0
	success_count = 0
	error_count = 0
	metadata_writer = csv.writer(output_file, delimiter = ',')
	for metadata_item in metadata : 
		file_id = metadata_item['file_id']
		file_name = metadata_item['file_name']
		total_count += 1 

		# check file existence

		# if not os.path.isfile('./' + file_id + '/' + file_name) :
		# 	print("warning: file {0} does not exist, please download again".format(file_id))
		# 	error_count += 1
		# 	continue
		
		# case info 

		cases = metadata_item['cases']
		if len(cases) :
			case_info = cases[0]
		else :
			print("[!]: cases info of file {0} is missing, skipped".format(file_id))
			error_count += 1
			continue

		primary_site = case_info['primary_site']
		disease_type = case_info['disease_type']

		# demographic info

		if 'demographic' in case_info:
			demographic_info = case_info['demographic']
			gender = (demographic_info['gender'] == 'male')
			race = demographic_info['race']
		else :
			print("[!]: demographic info of file {0} is missing, skipped".format(file_id))
			error_count += 1
			continue

		# diagnosis info

		diagnosis_info = case_info['diagnoses'][0]
		age = diagnosis_info['age_at_diagnosis']

		# sample info 

		sample_info = case_info['samples'][0]
		is_tumor = (sample_info['sample_type_id'] > 0 and sample_info['sample_type_id'] < 10)
		sample_type_id = sample_info['sample_type_id']
		sample_type = sample_info['sample_type']
		is_alive = (sample_info['state'] == 'live')

		metadata_writer.writerow([
			file_id,
			file_name,
			primary_site,
			disease_type,
			gender,
			age,
			race,
			is_tumor,
			sample_type_id,
			sample_type,
			is_alive])
		success_count += 1

	return total_count, success_count, error_count

if __name__ == "__main__":
	import argparse
	import os
	import json
	import csv

	# command line arguments

	parser = argparse.ArgumentParser(description = "Extract TCGA metadata")
	parser.add_argument("-m", "--metadata", dest = "metadata", required = True,
		help = "metadata.json downloaded from TCGA data repository", metavar = "FILE")
	parser.add_argument("-o", "--output", dest = "output", required = True,
		help = "output csv filename", metavar = "FILE")

	args = parser.parse_args()
	filename_metadata = args.metadata
	filename_output = args.output

	# check files existence

	if not os.path.isfile(filename_metadata) :
		print("[*] TCGA metadata.json file \"{0}\" does not exist!".format(filename_metadata))
		exit(-1)

	# open files

	with open(filename_metadata, "r") as file_metadata, \
		 open(filename_output, "wb") as file_output :
		metadata = json.load(file_metadata)
		total_count, success_count, error_count = extract_tcga_metadata(metadata, file_output)

	print("[*] complete: {0} items scanned, {1} items extracted, {2} item failed.".format(total_count, success_count, error_count))






