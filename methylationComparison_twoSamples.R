#! /usr/bin/env Rscript

############################################################################################
## This script compares the methylation information per CpG island for different samples ##  
## Output:                                                                               ##
## -  a plot comparing the average methylation per CpG island                            ##
## - the file "differential_coordinates.csv" containing the following information        ##
##      cgi_chro    cgi_start   cgi_end    cgi_id                                        ##
############################################################################################

setwd("/data/dessertlocal/daria/sepsisData_oncgnostic/")
setwd("/run/user/1000/gvfs/sftp:host=gsuffa.bioinf.uni-jena.de,user=nu36par")

meteore_results="/run/user/1000/gvfs/sftp:host=gsuffa.bioinf.uni-jena.de,user=nu36par/data/dessertlocal/daria/sepsisData_oncgnostic/methylation_calling/meteore_results/"
methylationFolder="/run/user/1000/gvfs/sftp:host=gsuffa.bioinf.uni-jena.de,user=nu36par/data/dessertlocal/daria/sepsisData_oncgnostic/methylation_calling/methylation_analysis/"
coverage = "/run/user/1000/gvfs/sftp:host=gsuffa.bioinf.uni-jena.de,user=nu36par/data/dessertlocal/daria/sepsisData_oncgnostic/coverage/"
mainFolder = "/run/user/1000/gvfs/sftp:host=gsuffa.bioinf.uni-jena.de,user=nu36par/data/dessertlocal/daria/"


## read in the average methylation per CpG island for both sample
t0044c_nanopolish = read.table(paste(mainFolder, "hncData_oncgnostics/t0044c/averageMethylation_t0044c_nanopolish_perCGI.txt", sep = ""), header = TRUE, sep = "\t")
ah01_nanopolish = read.table(paste(mainFolder, "sepsisData_oncgnostic/ah01/averageMethylation_ah01_nanopolish_perCGI.txt", sep = ""), header = TRUE, sep = "\t")
# ah01_nanopolish = read.table(paste(mainFolder, "sepsisData_oncgnostic/methylation_calling/methylation_analysis/averageMethylation_perCGI_control_nanopolish.txt", sep = ""), header = TRUE, sep = "\t")

## throw away all CpG isalnds with an coverage <3
ah01_nanopolish = ah01_nanopolish[which(ah01_nanopolish$Coverage>3),]
t0044c_nanopolish = t0044c_nanopolish[which(t0044c_nanopolish$Coverage>3),]
## delete the coverage column
t0044c_nanopolish$Coverage <- NULL
ah01_nanopolish$Coverage <- NULL
# ah01_nanopolish$cgi_chro <- NULL


## merge the two dataframes
t0044c_merged = merge(t0044c_nanopolish, ah01_nanopolish, by = c("cgi_id", "cgi_chro", "cgi_start", "cgi_end"))

## plot comparison of all CpG islands
png(paste(mainFolder, "hncData_oncgnostics/t0044c/plots/averageMethylation_vs_ah01.png", sep = ""), width = 800, height = 600)
plot(t0044c_merged$ah01_nanopolish, t0044c_merged$t0044c_nanopolish, pch = 16,
     xlab = "av. methylation ah01",
     ylab = "av. methylation T0044-C")
dev.off()

## select interesting CpG islands and write their coordinates to file
differential = t0044c_merged[which((t0044c_merged$ah01_nanopolish==0) & (t0044c_merged$t0044c_nanopolish>0.7)),]
differential_coordinates = differential[,c("cgi_chro", "cgi_start", "cgi_end", "cgi_id")]
write.table(differential_coordinates, paste(mainFolder, "hncData_oncgnostics/t0044c/cgiMethylations/differential_coordinates.csv", sep = ""), sep = ";", row.names = FALSE, quote = FALSE)


