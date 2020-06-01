# Summary

At its core the purpose of this code is to annotate ChIP peak regions to 
genes. Other parameters allow for users to annotate both the peaks and the 
genes using other types of datasets which we will describe further below.

## Dependencies

* Python
* Python Libraries
 * pandas
 * numpy
 * argparse
 * pybedtools

## Annotation Method

This script uses two functions to annotate peaks to genes which we call "round 1" and "round 2". 
* In round 1, this script annotates peaks by location relative to genes. 
 * First peaks within genes are annotated and then removed from further annotation
 * Then peaks upstream of genes (limited by `filter_tss_upstream`) are annotated and then removed from further annotation
 * Last peaks downstream of genes are annotated (limited by `filter_tts_downstream`). 
* Round 2 is dependent on genes that were found to be significantly differentially expressed (DE). 
 * Tab-delimited differential expression files are listed in the `round2ann` 
parameter. The must contain a column titled `gene_id`.
 * Genes called differentially expressed can be limited using the `rnaDESignificance` parameter
 * In round 2, peaks are annotated to upstream or downstream DE genes. Maximum distance that a peak can be annotated to a gene is limited by the `outlier_filter` parameter. You can set thresholds of DE significance (eg. maximum adjusted p-value, minimum fold-change) using the `rnaDESignificance` parameter

## A note about the peakwise and genewise output files

This script outputs two main files that contain the subtring "peakwise_ann" and "genewise_ann". The reason for the seperate files is due to the fact that annotations can be either peaks-based (peakwise), gene-based (genewise), or both. Multiple peaks can annotate to a single gene and some peaks may annotate to multiple genes (eg. if a peak overlaps multiple genes). Therefore peak-centric annotations (such as the location of a peak summit) will sometimes be reported in semicolon-delimited lists in "genewise_ann" output files, and genewise columns will sometimes be output as semicolon-delimited lists in "peakwise_ann" output files. If a column is both peak and gene centric then it will sometimes be reported in semicolon-delimited lists in both "peakwise_ann" and "genewise_ann" output files.

## Extra annotations

Given specific set parameters, this script outputs extra columns to further annotate the peak-gene relations. 

* **summit annotation** (both): If `narrowpeak_file` is set to a narrowPeak format file that corresponds exactly to `bed_file`, the script will annotate the peak's summits using round1 peak-to-gene annotation and compare the results. For example, if a peak is intragenic but its summit is downstream of the gene further than the threshold set by `filter_tts_downstream` then the peak will annotate to the gene, but the "summit_ann" column will report "False". Another example is if a peak overlaps two genes, but the summit is only intragenic in one then this column will report both "True" and "False" for each respective peak.
* **overlapping ChIP peak to bed features** (peakwise): If `comparePeaksToBeds` is set to bed/broadPeak/narrowPeak file(s), the script will annotate the current ChIP peaks by checking if any peaks overlap the features in `comparePeaksToBeds` files. The script will report True or False by default. Set `compBedScores` to report other characteristics about the feature that is overlapped. 
	* Overlap will be reported in a column with the respective name set in `comparePeaksToBedsNames`. The other characheteristics set in  `compBedScores` will be reported in a column named in the format [comparePeaksToBedsName]:[compBedScore].
* **overlapping ChIP summit-region to bed features** (peakwise): If `compareSumRegToBed` is set to bed/broadPeak/narrowPeak file(s), the script will annotate the current ChIP by comparing the region around the summit (specified by `summitRegion`) to the features defined in `compareSumRegToBed` files. In the peakwise ouput, the script will report the locations of the features in the given bed file(s) that overlap with the ChIP regions respectively in columns named in `compareSumRegToBedColNames`. In the genewise output, the columns will simply report True/False. 
	* Example of bed files that users can set in `compareSumRegToBed`:
		* Motif bed files (which we want to check overlaps within about +/-50bp summit region of the peak)
		* MNase bed file (which we want to check overlaps the summit)
* **overlapping ChIP peak to tsv features** (peakwise): If `comparePeaksToText` is set to text file(s), the script will annotate the peaks by checking if any peaks overlap the features from `comparePeaksToText` by performing a left-join on the columns: chrom, start, stop. By default, the script will report True or False if peaks overlap (>0bp) any row given in `comparePeaksToText`. Set `addPeaksToTextFeatures` to report other characteristics about the feature that is overlapped.
	* Overlap will be reported in a column with the respective name set in `comparePeaksToTextNames`. The other characheteristics set in  `compBedScores` will be reported in a column named in the format [comparePeaksToTextNames]:[addTextFeatures].
* **comparisons to other genewise annotations** (genewise): If `compareGenesToText` is set to text files(s), the script will annotate the genes by performing a left-join with the gene IDs annotated in the current sample to the "gene_id" column in the `compareGenesToText` files.
	* The respective column names given in `compareGenesToTextNames` will be output True if a gene was also reported in the other genewise file and false otherwise. 
	* To report other features from the external genewise file use the `addGeneToTextFeatures` parameter. These other features will be reported in a column named in the format [compareGenesToTextNames]:[addGeneToTextFeatures]
	* Examples of text files that users can set in `compareGenesToText`:
		* Tables that include common names, descriptions, or aliases for each gene
		* Tables that include gene expresion values (eg. TPM, RPKM)
		* Differential Expression Output (eg. DESeq2, EdgeR)
		* Genewise tables from previous runs of this script 
* **dnase annotation**: If `dnase_files` is set to bed file(s) then the script will report four different metrics:
	* **dnase-peak overlap** (peakwise): True/False whether the DNase regions overlap in column named by `dnase_name` 
	* **dnase-peak percent overlap** (peakwise):  The percentage of the ChIP peak that is overlapping with DNase region(s). Columns are named by "`dnase_name`:percentDNAseOverlapsPeak".
	* **dnase-gene percent overlap** (genewise):  The percentage of the gene that is overlapping with DNase region(s).

## Future Updates
In the future, we plan to add the following function/parameter(s) to this script:
* when comparing overlapping chip peak to bed features, there should be a way to report the bed features locations rather than just True when the bed feature(s) overlaps the peak
* when  comparing overlapping chip summit regions to bed features, there should be a way to report just True rather than the bed features locations when the bed feature(s) overlaps the summit region
* when  comparing overlapping chip summit regions to bed features, there should be a way to report more than the bed features locations of the bed feature(s) that overlaps the summit region (eg. the feature name given in the 4th column of the bed file)