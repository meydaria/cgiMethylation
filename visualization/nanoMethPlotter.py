#! /usr/bin/env python3

"""
This script extract methylation information from Nanopolish data.
Nanopolish data file names need to be given in a fofn with corresponding sample name. 
Position need to be given in bed file format.

@author: Daria Meyer
"""

import argparse
import os
import sys
import pandas as pd

parser = argparse.ArgumentParser(description='Visualize DNA methylation per CpG dinucleotide for a given region.')

parser.add_argument('-b', '--bed', 
					help='bed file with regions for which methylations will be generated', 
					type = str,
					required=True)

parser.add_argument('-i', '--infile', 
					help='fofn (file of filenames) containing tab-separated the sample name and the input file names of mehtylation data files', 
					type = str,
					required=True)

parser.add_argument('-o', '--outdir', help='Output directory name where methylation subfiles will be stored')

args = parser.parse_args()


# read in the bed file  --> nested list
bed = []
with open(args.bed, 'r') as fh:
	for line in fh.readlines():
		if(line.strip()):
			bed.append(line.split())

# extract sample names and file names from the fofn  --> dictionary with key = sample name, value = file name
methFiles = {}
with open(args.infile, 'r') as fh:
	for line in fh.readlines():
		# add the file only if it exists
		if(line.strip() and not line.startswith('#')):
			sample = line.split()[0]
			filename = line.split()[1]
			if(os.path.exists(filename)):
				methFiles[sample] = filename
			else:
				sys.exit(f"The file {filename} does not exist")

# set the output directory to the passed argument, or to the local directory if no argument was passed
outdir = os.getcwd()
if args.outdir is not None:
	if(os.path.exists(args.outdir)):
		outdir = args.outdir
		print(f"Output will be wirtten to {outdir}")
	else:
		print(f"Warnning: {args.outdir} does not exist! Output will be wirtten to {outdir}")




def get_regioMeth(methFiles, chro, start, stop):
	"""
	The function create_df reads the methylation data from the different samples and concatenates them into one pandas dataframe
	Currently input is methylation information per CpG
	input: 
		dictionary of file names with 
			key = sample name
			value = file name
		one region for which the methylation information will be extracted
	output: 
		pandas dataframe with methylation information for all listed samples
	"""
	li = []  # list of pandas dataframes containing regional methylation informaiton
	for sample, filename in methFiles.items():
		## extract all lines which match the given coordinates (from bed file)
		lines = [] 
		with open(filename, 'r') as fh:
			print(f"Reading {filename}...")
			# Skip first line (header)
			line = fh.readline()
			while True:
				line = fh.readline()  # read file line by line to not crash the RAM
				if not line:
					break
				splitted = line.split()
				if(splitted[0]==chro and int(splitted[1]) > int(start) and int(splitted[1]) < int(stop)):  # check if in region of interest
					lines.append(splitted)	

		## create a pandas dataframe of the resulting lines	
		tmp = pd.DataFrame(lines, columns =['chromosome', 'position', 'strand', 'methylation', 'read_id'])
		tmp = tmp.astype({'position': float, 'methylation':float})
		
		## for R: add a column containing the sample name, needed for R
		tmp['sample'] = sample
		
		li.append(tmp)
	## concatenate the different dataframes into one big dataframe 
	print("Concatenating dataframes...")
	df = pd.concat(li, axis=0, ignore_index=True)
	return(df)

	


## MAIN 
#******************************************************************************************************************************
for region in bed:
	## extract methylation information for one region from all samples
	df = get_regioMeth(methFiles, region[0], region[1], region[2])
	
	## check if the region of interest has a name given in the fourth column
	if(len(region)>3):
		name = region[3]
	else:
		name = region[0] + "_" + region[1] + "_" + region[2]
	
	df.to_csv(os.path.join(outdir, name+".csv"))

