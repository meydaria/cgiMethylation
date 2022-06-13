#! /usr/bin/env Rscript

############################################################################################
## This script compares the methylation information per CpG island for different samples ##  
## Input:                                                                                ## 
##    - 2x tab-separated file with average methylation per CpG island                    ##
##    - 2x sample name                                                                   ## 
##    - optional: output folder                                                          ##
## Output:                                                                               ##
##    - a plot comparing the average methylation per CpG island                          ##
##    - the file "differential_coordinates.csv" containing the following information     ##
##      cgi_chro    cgi_start   cgi_end    cgi_id                                        ##
############################################################################################

# Rscript methylationComparison_twoSamples.R averageMethylation_ah01_nanopolish_perCGI.txt averageMethylation_t0044c_nanopolish_perCGI.txt ah01 t0044c

# get command line arguments
args = commandArgs(trailingOnly=TRUE)


# test if there is at least one argument: if not, return an error
if (length(args)<4) {
  stop("Four input files have to be given: firstly the two methylation files, than the two sample names respectively", call.=FALSE)
} else if(length(args)==4){
  # set default directory for output files
  args[5] = getwd()
}

## read in the average methylation per CpG island for both sample
control = read.table(args[1], header = TRUE, sep = "\t")
tumor = read.table(args[2], header = TRUE, sep = "\t")
control.name = args[3]
tumor.name = args[4]
outdir = args[5]

## throw away all CpG isalnds with an coverage <3
control = control[which(control$coverage>3),]
tumor = tumor[which(tumor$coverage>3),]

## remove the coverage column
tumor$coverage <- NULL
control$coverage <- NULL

colnames(control)[(!colnames(tumor)%in%c("cgi_id", "cgi_chro", "cgi_start", "cgi_end"))] = "control"
colnames(tumor)[(!colnames(tumor)%in%c("cgi_id", "cgi_chro", "cgi_start", "cgi_end"))] = "tumor"

## merge the two dataframes
mergedData = merge(control, tumor, by = c("cgi_id", "cgi_chro", "cgi_start", "cgi_end"))

print(head(mergedData))
print(head(mergedData$control.name))
print(outdir)

## plot comparison of all CpG islands
png(paste(outdir, "/averageMethylation_plot.png", sep = ""), width = 800, height = 600)
plot(mergedData$control, mergedData$tumor, pch = 16,
     xlab = paste("av. methylation ", control.name, sep = ""),
     ylab = paste("av. methylation ", tumor.name, sep = ""),
     main = "Average methylation per CpG island")
dev.off()

## select interesting CpG islands and write their coordinates to file
differential = mergedData[which(((mergedData$control==0) & (mergedData$tumor>0.7)) | ((mergedData$control>0.7) & (mergedData$tumor==0))),]
# differential = mergedData[which((mergedData$control==0) & (mergedData$tumor>0.7)),]  # only one way round
differential_coordinates = differential[,c("cgi_chro", "cgi_start", "cgi_end", "cgi_id")]
write.table(differential_coordinates, paste(outdir, "/differentialMethylation_coordinates.csv", sep = ""), sep = ";", row.names = FALSE, quote = FALSE)


