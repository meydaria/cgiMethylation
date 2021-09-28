#! /bin/python3 

"""
This script returns the average methyaltion frequency per CpG island for a given meteore result file. 
The returned value is the weighted average of methylation probability, which is calculated among all reads, that cover the CpG islands. 
The weight is the coverage.

Take care: this script runs very long! (Multiple days)

Contact:
daria.meyer@uni-jena.de

Dependencies:
    docopt >= 0.6.2

Usage:
	 get_cgiMethylation_megalodon.py --sample=<sample> --name=<name> --outfile=<outfile>
	 get_cgiMethylation_megalodon.py

Required Arguments:
	 --sample <sample>                        	Path to the freq-perCG.tsv file, that contains the methylation frequency per CpG site. 
	 --name <name>								Sample name
	 --outfile <outfile>                        Outfile in which the result will be stored

Options:
	 -h, --help                          Show this help message and exit.

Example:
	get_cgiMethylation_megalodon.py --sample /meteore/megalodon_results/hnoPool_megalodon-freq-perCG.tsv --name hnoPool_megalodon --outfile averageMethylation_megalodon_perCGI.txt
"""



import os
import sys
import time
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from docopt import docopt
args = docopt(__doc__)

sample_infile = args['--sample']
sample_name = args['--name']	
outfile = args['--outfile']	



### read in the data
print("Reading in the dataframes...")

## CpG islands
cgis_infile="/home/nu36par/projects/cpgIslands/humanCpG_islands.bed"
cgis = pd.read_csv(cgis_infile, sep="\t", names=["chromosome", "start", "stop", "cgi_ID", "CpG_num"]) 
## change the format of chromosome column by deleting the leading "chr"
# cgis["chromosome"] = cgis["chromosome"].str[3:]

## Methylation data
# hnoPool_megalodon="/data/dessertlocal/daria/hncData_oncgnostcs/meteore/megalodon_results/hnoPool_megalodon-freq-perCG.tsv"
# hnoPool_nanopolish="/data/dessertlocal/daria/hncData_oncgnostcs/meteore/nanopolish_results/hnoPool_nanopolish-freq-perCG.tsv"

## format the data
## hnoPool_megalodon = pd.read_csv(hnoPool_megalodon, sep="\t", names=["chromosome", "start_position", "stop_position", "coverage", "methylation_frequency", "strand"], header=0, dtype={0: str}, nrows = 1000) 
sample_methylations = pd.read_csv(sample_infile, sep="\t", dtype={0: str}) 
sample_methylations.name = sample_name
# hnoPool_nanopolish = pd.read_csv(hnoPool_nanopolish, sep="\t", dtype={0: str}) 
# hnoPool_nanopolish.name = 'hnoPool_nanopolish'

print("sample name= " + sample_name)
print("dataframe: ")
print(sample_methylations)

print("## Finished reading dataframes. Starting the calculation now...")
# setting start time to know how long the calculations take.
total_start_time = time.time()

## initialize an empty array which in the end will contain a dictionary for each combination of METEORE sample and CpG island
d = []		
## iterate over each CpG island
for (index, cgi_chro, cgi_start, cgi_end, cgi_id, cpg_num) in cgis.itertuples():
	line = {}
	cgi_id = cgi_id.strip()
	## add this line, to get information about the CpG island
	line["cgi_id"]= cgi_id 
	line["cgi_chro"]= cgi_chro
	line["cgi_start"]= cgi_start
	line["cgi_end"]= cgi_end
	# print(line)
    ## iterate over the different meteore result files
	# for sample in [hnoPool_megalodon, hnoPool_nanopolish]:
	# loop_start_time = time.time()
	## extract only those positions that match to CGIs
	short_SN = sample_methylations.query('Chr == @cgi_chro & @cgi_start <= Pos_start & Pos_start <= @cgi_end')
	## if data exists for this METEORE sample for the given CpG island
	if not short_SN.empty:
		try:
			weighted_average = np.average(short_SN['Methylation'], weights=short_SN['Coverage'])
		except:
			# write NaN if there are no values for this sample
			result = pd.DataFrame(d)
			result.to_csv('tmp_average_methylation.txt', index = False)
			weighted_average = float("NaN")
		line[sample_methylations.name] = weighted_average
	## append the dictionary containing data for one sample (eg. hnoPool_megalodon) to d
	# print(line)
	d.append(line)

	# print(d)
    ## get a intermediate result after 20 entries
	# if(cgi_id == " cgi_1"):
	# 	result = pd.DataFrame(d)
	# 	result.to_csv('/data/dessertlocal/daria/hncData_oncgnostcs/tmp_average_methylation.txt', index = False)

	# loop_time = time.time() - loop_start_time
	# total_time = time.time() - total_start_time

	# print(f'### {cgi_id} finished. Took {total_time} so far.')
	# now = time.time() - start_time
	# print(f'finished cgi_island {cgi_id} took {now}')

result = pd.DataFrame(d)
result.to_csv(outfile, index = False, sep = "\t")

total_time = time.time() - total_start_time
print(f'### Finished. Took {total_time} in total.')



### alternative ways to get shortSN
	# short_SN = sample[(sample['Chr']==cgi_chro)&(sample['Pos'].between(cgi_start, cgi_end, inclusive=True))]
	# short_SN = sample.loc[(sample['Chr']==cgi_chro)&(sample['Pos'].between(cgi_start, cgi_end, inclusive=True))]	

	# # look only at those entries on the same chromosome
	# first = time.time() - start_time
	# short_SN = sample[sample['Chr']==cgi_chro]
	# now = time.time() - first
	# print(f'Subset dataframe to Chromosome took {now}')

	# # look only at the position of the read which matches a the
	# first = time.time() - start_time
	# short_SN = short_SN[short_SN['Pos'].between(cgi_start, cgi_end, inclusive=True)]
	# now = time.time() - first
	# print(f'Subset dataframe took Positions {now}')


# Dict = {}
# # dictionary for each cpg
# for (index, cgi_chro, cgi_start, cgi_end, cpg_num, cgi_id) in cgis.itertuples():
# 	# print("####" + cgi_id + " on chro " + cgi_chro)
# 	# create a new dictionary for this CGI
# 	Dict[cgi_id]={}
# 	# print(Dict)
# 	## hnoPool and Nanopolish
# 	for sample in [hnoPool_megalodon, hnoPool_megalodon, ahControl_nanopolish, ahControl_megalodon]:
# 		# print(sample.name)
# 		# look only at those entries on the same chromosome
# 		short_SN = sample[sample['Chr']==cgi_chro]
# 		# print(short_SN)
# 		# get all CpG positions from hnoPool nanopolish and add them
# 		short_SN = short_SN[short_SN['Pos'].between(cgi_start, cgi_end, inclusive=True)]
# 		# print(short_SN)
# 		if not short_SN.empty:
# 			Dict[cgi_id][sample.name] = {}
# 			for (idx, chro, pos, methylFreq, cov) in short_SN.itertuples():
# 				Dict[cgi_id][sample.name][pos] = [methylFreq, cov]		
# 		# print(Dict)

# # print(Dict)
# # print(Dict.keys())
# resultList = []
# for key in Dict:
# 	# if there exists data for the CGI
# 	if(Dict[key]):
# 		# print(f'{key} \t {Dict[key].keys()}')
# 		# print(f'### This is key {key}')
# 		for sample in Dict[key]:
# 			# print(f'# sample {sample}')
# 			# print(len(Dict[key][sample].keys()))
# 			# print(len(Dict[key][sample].values()))
# 			total_methylFreq = 0
# 			total_coverage = 0
# 			for row in Dict[key][sample].items():
# 				positions = row[0]
# 				methylFreq = row[1][0]
# 				coverage = row[1][1]
# 				# calculate the weighted average
# 				total_methylFreq += methylFreq*coverage
# 				total_coverage += coverage
# 			if total_coverage > 0:
# 				total_methylFreq = total_methylFreq/total_coverage
# 				# print(total_methylFreq)
# 				# print(total_coverage)
# 				resultList.append([key, sample, total_methylFreq, total_coverage])
# 
# # create a tab seperated output containing in each row the following information
# # cgi_ID  |	sample_name   |  averaged methylation for complete CGI   |    summed coverage over all positions ### TODO: maybe change this? 		
# for entry in resultList:
# 	print(f'{entry[0]} \t {entry[1]} \t {entry[2]} \t {entry[3]}')
	
# TODO: Next steps:
# visualize the data -> stacked bar chart per cgi? color code per cgi?


