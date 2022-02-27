#!/bin/bash

# This script does some first basic processing steps after receiving the data from ara. 
# Input: 
# 	- foldername containing: fastqs, guppy_logs, sequencing_summary.txt
# Output:
# 	- new folder containing the guppy logs
# 	- one fastq file (containing all the concatenated fastqs)
# 	- a sorted bam file contianing the alignment

coverageBool=false

while getopts hcn:f: opt
do
	case $opt in
		h )
      		echo "Usage: mapping.sh [-c] -n -f <folder> "
      		echo "     -h           Display this help message and quit"
      		echo "     -f           specify the name of the folder in which the data is located."
      		echo "     -n           specify the name of the fastq, which will be used as sample name."
      		echo "     -c           optional parameter; If -c is set, the coverage will be calculated."
      		exit 0
      		;;
		f) DIR=$OPTARG;;
		n) SAMPLE=$OPTARG;;
		c) coverageBool=true;;
    	\?) 
			echo "Invalid Option: -$OPTARG" 1>&2
      		exit 1
      		;;
   esac
done

GENOME='/data/fass1/genomes/Eukaryots/homo_sapiens_done/ucsc/hg38.fa'
cpgIslands="/data/dessertlocal/daria/cpgIslands"
scripts="/home/nu36par/scripts"


# Create a folder guppy_logs and move all guppy logs in there
if [ -d "${DIR}/guppy_logs" ]; then
  ### Warning: Folder already exists ###
  echo "Warning: ${DIR}/guppy_logs exists already. Cannot move Guppy log files."
  echo "Skip moving guppy logs."
else
  ###  Control will jump here if $DIR does NOT exists ###
  echo "Copying guppy logs to new folder ${DIR}/guppy_logs..."
  mkdir "${DIR}/guppy_logs"
  mv $DIR/guppy_basecaller_log*.log $DIR/guppy_logs/
fi

# # check whether an outfile exists already which would be overwritten
# if [[ -f "${DIR}/output.fastq" ]]; then
# 	### Error if file already exists ###
# 	echo "Warning: The file ${DIR}/output.fastq already exists already. This file will be used for the following steps!"
# 	echo "Skip concatenating fastq files!"
# else
# 	###  Concatenate all fastq files ###
# 	echo "## Concatenating fastq files"
# 	mv $DIR/${SAMPLE}.fastq > $DIR/output.fastq
# 	# rm $DIR/fastq_runid*.fastq
# fi


## mapping
if [[ -f "${DIR}/${SAMPLE}.sorted.bam" ]]; then
	### Error if file already exists ###
	echo "Warning: The file ${DIR}/${SAMPLE}.sorted.bam exists already. This file will be used for the following steps!"
	echo "Skip mapping!"
else
	if [[ -f "${DIR}/${SAMPLE}.fastq" ]]; then
		echo "## Mapping with minimap2 ..."
		nice minimap2 -a -x map-ont $GENOME "${DIR}/${SAMPLE}.fastq" | samtools sort -T tmp -o "${DIR}/${SAMPLE}.sorted.bam"
		samtools index "${DIR}/${SAMPLE}.sorted.bam"
	fi
fi

## quality control
if [[ -f "${DIR}/sequencing_summary.txt" ]]; then
	source /home/nu36par/miniconda3/etc/profile.d/conda.sh
	conda activate pycoQC
	echo "## Quality Control with pycoQC ..."
	pycoQC -f "${DIR}/sequencing_summary.txt" -a "${DIR}/${SAMPLE}.sorted.bam" -o "${DIR}/pycoQC.html"
else
	echo "Warning: The file ${DIR}/sequencing_summary.txt does not exist. Therefore no pycoQC-Report can be generated"
	echo "Skip quality contol!"
fi


# get amount of reads
echo "Amount of basecalled reads:"
awk '{s++}END{print s/4}' "${DIR}/${SAMPLE}.fastq"  # Anzahl der Zeilen /4
# get amount of bases
echo "Amount of basecalled bases:"
awk 'BEGIN{sum=0;}{if(NR%4==2){sum+=length($0);}}END{print sum;}' "${DIR}/${SAMPLE}.fastq"


## calculate coverage if wanted
if "$coverageBool"; then
	echo "## Calculating the coverage..."
	mosdepth --no-per-base --by $cpgIslands/humanCpG_islands.bed ${SAMPLE} $readData/${SAMPLE}.sorted.bam  
	# get coverage of the human genome: 
	echo "Coverage of the human genome on average:"
	tail -n2 t0004c_versuch_2022_0056.mosdepth.summary.txt | head -n 1 | cut -f4
	## Coverage per CpG island
	echo "Coverage over all CpG islands:"
	tail -n1 t0004c_versuch_2022_0056.mosdepth.summary.txt | cut -f4
fi 




# ## calculate coverage if wanted
# if "$coverageBool"; then
# 	echo "## Calculating the coverage..."
# 	# get coverage on CpG islands
# 	echo "Coverage of CpG islands on average:"
# 	samtools depth -a -b $cpgIslands/humanCpG_islands.bed $DIR/${SAMPLE}.sorted.bam | awk 'BEGIN{SUM=0}{sum+=$3}END{print sum/24200434}'  
# 	# get coverage of the human genome: 
# 	echo "Coverage of the human genome on average:"
# 	hg38_size=$(grep -v ">" $GENOME | wc | awk '{print $3-$1}')          # 3257347282
# 	samtools depth "${DIR}/${SAMPLE}.sorted.bam" | awk -v hsa_size=$hg38_size '{sum+=$3} END { print "Sum = ",sum,"Average = ",sum/hsa_size}' 
# 	## Coverage per CpG island
# 	nice bash $scripts/getCoverage_samtoolsUCSC.sh -p $cpgIslands/humanCpG_islands.bed -b "${DIR}/${SAMPLE}.sorted.bam" -o "${DIR}/cgiCoverage.txt"
# 	# get amount of reads
# 	echo "Amount of basecalled reads:"
# 	awk '{s++}END{print s/4}' "${DIR}/${SAMPLE}.fastq"  # Anzahl der Zeilen /4
# 	# get amount of bases
# 	echo "Amount of basecalled bases:"
# 	awk 'BEGIN{sum=0;}{if(NR%4==2){sum+=length($0);}}END{print sum;}' "${DIR}/${SAMPLE}.fastq"
# fi 