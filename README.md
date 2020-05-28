# Summary

At its core the purpose of this code is to annotate ChIP peak regions to 
genes. Other parameters allow for users to annotate both the peaks and the 
genes using other types of datasets which we will describe further below.

## Dependencies

* Bedtools
* Samtools
* Python
* Python Libraries
 * pandas
 * numpy
 * argparse

## Annotation Method

This script uses two functions to annotate peaks to genes which we call 
"round 1" and "round 2". 
* In round 1, this script annotates peaks by location relative to genes. 
 * First peaks within genes are annotated and then removed from further 
 annotation
 * Then peaks upstream of genes (limited by `filter_tss_upstream`) are 
 annotated and then removed from further annotation
 * Last peaks downstream of genes are annotated (limited by 
 `filter_tts_downstream`). 
* Round 2 is dependent on genes that were found to be significantly 
differentially expressed (DE). 
 * Tab-delimited differential expression files are listed in the `round2ann` 
parameter. The must contain a column titled `gene_id`.
 * Genes called differentially expressed can be limited using the 
 `rnaDESignificance` parameter
 * In round 2, peaks are annotated to upstream or downstream DE 
genes. Maximum distance that a peak can be annotated to a gene is limited 
by the `outlier_filter` parameter. 

## A note about the peakwise and genewise output files

This script outputs two main files that contain the subtring "peakwise_ann" and 
"genewise_ann". The reason for the seperate files is due to the fact that 
annotations can be either peaks-based (peakwise), gene-based (genewise), 
or both. Multiple peaks can annotate to a single gene and some peaks may annotate to multiple genes (eg. if a peak overlaps multiple genes). Therefore peak-centric 
annotations (such as the location of a peak summit) will sometimes be reported 
in semicolon-delimited lists in "genewise_ann" output files, and genewise 
columns will sometimes be output as semicolon-delimited lists in "peakwise_ann" output files. If a column in both peak and gene centric then it will sometimes
be reported in semicolon-delimited lists in both "peakwise_ann" and 
"genewise_ann" output files.

## Extra annotations

Given specific set parameters, this script outputs extra columns to 
further annotate the peak-gene relations. 

* summit annotation (both): if `narrowpeak_file` is set, the script will 
annotate the peak's summits using round1 peak-to-gene annotation and 
compare the results. For example, if a peak is intragenic but its 
summit is downstream of the gene further than the threshold set 
by `filter_tts_downstream` then the peak will annotate to the gene, 
but the "summit_ann" column will report "False". Another example is if 
a peak overlaps two genes, but the summit is only intragenic in one then
this column will report both "True" and "False" for each respective peak.


This script then annotates ChIP peaks by overlapping them with
any of the following: 
* other ChIP peaks regions
* DNase regions
* MNase regions
* Motifs

Last, text files are output that contain peak-centric and gene-centric annotations.


