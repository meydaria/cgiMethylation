#! /usr/bin/env Rscript

############################################################################################
## This script compares the methylation information per CpG island for different samples ##  
## Input:                                                                                ## 
##    - tab-separated file with average methylation per CpG island for control sample    ##
##    - tab-separated file with average methylation per CpG island for tumor sample      ## 
##    - optional: output folder                                                          ##
## Output:                                                                               ##
##    - a plot comparing the average methylation per CpG island                          ##
##    - the file "differential_coordinates.csv" containing the following information     ##
##      cgi_chro    cgi_start   cgi_end    cgi_id                                        ##
############################################################################################

# get command line arguments
args = commandArgs(trailingOnly=TRUE)


# test if there is at least one argument: if not, return an error
if (length(args)<2) {
  stop("Two input files have to be given: first the control, than the tumor sample", call.=FALSE)
} else if(length(args)==2){
  # set default directory for output files
  args[3] = getwd()
}

## read in the average methylation per CpG island for both sample
control = read.table(args[1], header = TRUE, sep = "\t")
tumor = read.table(args[2], header = TRUE, sep = "\t")
outdir = args[3]

## throw away all CpG isalnds with an coverage <3
control = control[which(control$coverage>3),]
tumor = tumor[which(tumor$coverage>3),]

## remove the coverage column
tumor$coverage <- NULL
control$coverage <- NULL

## merge the two dataframes
mergedData = merge(tumor, control, by = c("cgi_id", "cgi_chro", "cgi_start", "cgi_end"))

## plot comparison of all CpG islands
png(paste(outdir, "averageMethylation_plot.png", sep = ""), width = 800, height = 600)
plot(mergedData$control, mergedData$tumor, pch = 16,
     xlab = "av. methylation control",
     ylab = "av. methylation tumor",
     main = "Average methylation per CpG island")
dev.off()

## select interesting CpG islands and write their coordinates to file
differential = mergedData[which((mergedData$control==0) & (mergedData$tumor>0.7)),]
differential_coordinates = differential[,c("cgi_chro", "cgi_start", "cgi_end", "cgi_id")]
write.table(differential_coordinates, paste(outdir, "differentialMethylation_coordinates.csv", sep = ""), sep = ";", row.names = FALSE, quote = FALSE)


