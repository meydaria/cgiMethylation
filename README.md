# cgiMethylation
This repository contains some scripts that I use to analyze ONT sequencing results, for which methylations were called using METEORE.


# getCgiMethylation.py
The script returns the weighted average methyaltion frequency per CpG island for a given meteore result file. Weighted, because the amount of reads that support the methylation information is used to weigh the methylation information.

