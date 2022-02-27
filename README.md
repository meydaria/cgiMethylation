# cgiMethylation
This repository contains some scripts that I use to analyze ONT sequencing results, for which methylations were called using METEORE.


# getCgiMethylation.py
The script returns the weighted average methyaltion frequency per CpG island for a given meteore result file. Weighted, because the amount of reads that support the methylation information is used to weigh the methylation information.

# mapping.sh
This script does the first basic steps after basecalling
- merge all guppy output files into one folder
- mapping with minimap2 (https://github.com/lh3/minimap2)
- quality control with pycoQC (https://github.com/a-slide/pycoQC)
- coverage calculation with mosdepth (https://github.com/brentp/mosdepth)
