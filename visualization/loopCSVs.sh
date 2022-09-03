#!/bin/bash

# Helper script to loop over all created CSV files and generate the mehtylation plots for them.

folder=$1
outfolder=$2

files="${folder}/*.csv"
for f in $files
do
  echo "Processing $f file..."
  # take action on each file. $f store current file name
  Rscript /home/nu36par/scripts/ontAnalyses/visualization/plotMethylation.R -f $f -o $outfolder
done
