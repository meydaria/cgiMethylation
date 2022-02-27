#! /bin/python3 

"""
This script returns the average methyaltion frequency per CpG island for a given meteore result file. 
The returned value is the weighted average of methylation probability, which is calculated among all reads, that cover the CpG islands. 
The weight is the coverage.

The file freq-perCG.tsv does not necessarily need to be a meteore output file. Important is only that is has the right format.
It has to be a tab-seperated file with the following columns: "Chr", "Pos", "Methyl_freq", "Cov"
Containing information for the CpGs:
1. The chromosome on which it is located (take care "1" vs "chr1" etc.).
2. The position of the C from the CpG.
3. The methylation frequency (ranging between 0 and 1).
4. The coverage at this position.

Take care: this script runs very long! (Multiple days)

TODO: implement more efficiently

Contact:
daria.meyer@uni-jena.de

Dependencies:
    docopt >= 0.6.2

Usage:
	 getCgiMethylation.py --cgis=<cgis> --sample=<sample> --name=<name> --outfile=<outfile>
	 getCgiMethylation.py

Required Arguments:
	 --cgis <cgis>						Path to the bed file, which contains the locations of the CpG islands of interest.
	 --sample <sample>					Path to the freq-perCG.tsv file, that contains the methylation frequency per CpG site. 
	 --name <name>						Sample name
	 --outfile <outfile>                Outfile in which the result will be stored

Options:
	 -h, --help                          Show this help message and exit.

Example:
	get_cgiMethylation_megalodon.py --cgis /home/nu36par/projects/cpgIslands/humanCpG_islands.bed --sample /meteore/megalodon_results/hnoPool_megalodon-freq-perCG.tsv --name hnoPool_megalodon --outfile averageMethylation_megalodon_perCGI.txt
"""


# import sys
# import matplotlib
# import matplotlib.pyplot as plt
import os
import time
import pandas as pd
import numpy as np
from docopt import docopt
args = docopt(__doc__)

cgis_infile = args['--cgis']
sample_infile = args['--sample']
sample_name = args['--name']	
outfile = args['--outfile']		


### read in the data
print(f"Reading in the CPG islands from {cgis_infile} ...")

## CpG islands
# cgis_infile="/home/nu36par/projects/cpgIslands/humanCpG_islands.bed"
cgis = pd.read_csv(cgis_infile, sep="\t", names=["chromosome", "start", "stop", "cgi_ID", "CpG_num"]) 
## change the format of chromosome column by deleting the leading "chr" ## needed for the older version of meteore
# cgis["chromosome"] = cgis["chromosome"].str[3:]

## Methylation data
print("Reading in the methylation sample data...")
sample_methylations = pd.read_csv(sample_infile, sep="\t", dtype={0: str}, skipinitialspace = True, names=["Chr", "Pos", "second_Pos", "Coverage", "Methylation"], header = 0) 
sample_methylations.name = sample_name

print("sample name = " + sample_name)
print("dataframe: ")
print(cgis.head(10))
print(sample_methylations.head(10))
print(list(sample_methylations.columns))
print(list(cgis.columns))

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
	# print(f'line: {line}')

    ## iterate over the different meteore result files
	# for sample in [hnoPool_megalodon, hnoPool_nanopolish]:
	# loop_start_time = time.time()
	
	## extract only those positions that match to CGIs
	# short_SN = sample_methylations.query('Chr == @cgi_chro & @cgi_start <= Pos_start & Pos_start <= @cgi_end')
	short_SN = sample_methylations.query('Chr == @cgi_chro & @cgi_start <= Pos & Pos <= @cgi_end')
	# print(f'sample_methylations: {sample_methylations.head(5)}')
	# print(f'Short SN: {short_SN}')
	
	
	## if data exists for this METEORE sample for the given CpG island
	if short_SN.empty:
		line[sample_methylations.name] = float("NaN")
		line["Coverage"] = 0
	else:
		try:
			# weighted_average = np.average(short_SN['Methyl_freq'], weights=short_SN['Cov'])
			weighted_average = np.average(short_SN['Methylation'], weights=short_SN['Coverage'])
			coverage = np.mean(short_SN['Coverage'])
		except:
			# write NaN if there are no values for this sample
			result = pd.DataFrame(d)
			result.to_csv('tmp_average_methylation.txt', index = False)
			weighted_average = float("NaN")
			coverage = 0
		line[sample_methylations.name] = weighted_average
		line["Coverage"] = coverage
		# line["Coverage"] = np.mean(short_SN['Cov'])
	## append the dictionary containing data for one sample (eg. hnoPool_megalodon) to d
	# print(f'line: {line}')
	d.append(line)

	# print(d)
    ## get a intermediate result after 20 entries
	# if(cgi_id == "cgi_6"):
	# 	result = pd.DataFrame(d)
	# 	result.to_csv(os.path.dirname(outfile) + "/tmp.txt", index = False, sep = "\t", na_rep = "NaN")

	# loop_time = time.time() - loop_start_time
	# total_time = time.time() - total_start_time
	# print(f'### {cgi_id} finished. Took {total_time} so far.')

result = pd.DataFrame(d)
result.to_csv(outfile, index = False, sep = "\t")

total_time = time.time() - total_start_time
print(f'### Finished. Took {total_time} in total.')





