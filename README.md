# Summary
The purpose of this repo is to annotate ChIP peak regions that have already 
been called. 

The first round annotates peaks by location relative to genes. 
Priority is given to intragenic peaks and then peaks upstream of 
genes (limited by `filter_tss_upstream`) and then
peaks downstream of genes (limited by `filter_tts_downstream`). The second 
round is dependent on genes that were found to be significantly differentially 
expressed (DE). In round 2, peaks are annotated to upstream or downstream DE 
genes (limited by `outlier_filter`). 

This script then annotates ChIP peaks by overlapping them with
any of the following: 
* other ChIP peaks regions
* DNase regions
* MNase regions
* Motifs

Last, text files are output that contain peak-centric and gene-centric annotations.

# Python Library Dependencies
* pandas
* numpy
* argparse
