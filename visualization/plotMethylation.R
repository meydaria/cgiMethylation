#!/usr/bin/env Rscript
library("optparse")
library(ggplot2)

## Process the command line arguments
################################################################################
option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-o", "--outfolder"), type="character", default="getwd()", 
              help="output folder name [default= %default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

sample = strsplit(basename(opt$file), "[.]")[[1]][1]

## Functions
################################################################################

plot_together <- function(df){
    ## Generate the dot line plot per CpG-island
    outfile = paste(opt$out,paste0(sample, ".png"), sep ="/")
    plot <- ggplot(df, aes(position, methylation, group = read_id, colour = sample)) + 
      geom_line() + 
      geom_point() +
      ggtitle(paste("Methylation Plot for ", sample, sep = "")) +
      xlab(paste("Position on ", df$chromosome[1], sep = "")) 
    ggsave(outfile, plot, "png", width = 8, height = 6)
}

plot_single <- function(df){
  ## Generate the dot line plot per CpG-island
  outfile = paste(opt$out,paste0(sample, "_single.png"), sep ="/")
  plot <- ggplot(df, aes(position, methylation, group = read_id, colour = sample)) + 
    geom_line() + 
    geom_point() +
    ggtitle(paste("Methylation Plot for ", sample, sep = "")) +
    facet_wrap("sample") + 
    xlab(paste("Position on ", df$chromosome[1], sep = "")) +
    theme(axis.text.x = element_text(angle=90, vjust=1, hjust=1))
  ggsave(outfile, plot, "png", width = 8, height = 6)
}



## program...
df = read.table(opt$file, header=TRUE, sep = ",")
print(head(df))
print(dim(df))
plot_together(df)
plot_single(df)
print(opt$out)
