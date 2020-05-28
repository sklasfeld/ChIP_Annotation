#!/usr/bin/env python

import os
import sys
import pandas as pd
import numpy as np
import pybedtools


def compare_bedfiles(bedfile1, bedfile2, *positional_parameters, \
	**keyword_parameters):
	"""compare two bedfiles to see if they overlap

	mandatory parameters:
	* bedfile1: bed file 
	* bedfile2: bed file
	
	optional parameters:
	* distance: Minimum number of bp must overlap (Default:0)
	* keep_tmps: keep temp files made in this
	* verbal: print number of peaks that overlap (Default: False)"""

	distance=0
	verbal=False
	keep_tmps=False

	if ('distance' in keyword_parameters):
		distance = keyword_parameters['distance']
	if ('verbal' in keyword_parameters):
		verbal = keyword_parameters['verbal']
	if ('keep_tmps' in keyword_parameters):
		keep_tmps = keyword_parameters['keep_tmps']

	if not os.path.isfile(bedfile1) or os.stat(bedfile1).st_size == 0:
		sys.exit((("\nERROR: Cannot find file: %s\n") % bedfile1))
	if not os.path.isfile(bedfile2) or os.stat(bedfile2).st_size == 0:
		sys.exit((("\nERROR: Cannot find file: %s\n") % bedfile1))

	file1_df = pd.read_csv(bedfile1,sep='\t', header=None, dtype=str)
	ncols_file1 = file1_df.shape[1]
	if ncols_file1 > 10:
		sys.exit((("\nERROR: %s is not in BED format nor NARROWPEAK format!!! " + \
			"Must be tab delimited with 10 columns maximum!!!\n") % bedfile1))

	file2_df = pd.read_csv(bedfile2,sep='\t', header=None, dtype=str)
	ncols_file2 = file2_df.shape[1]
	if ncols_file2 > 10:
		sys.exit((("\nERROR: %s is not in BED format nor NARROWPEAK format!!! " + \
			"Must be tab delimited with 10 columns maximum!!!\n") % bedfile2))
	
	chrInFile1 = list(file1_df.iloc[:,0].unique())
	chrInFile2 = list(file2_df.iloc[:,0].unique())
	inFile1Not2 = [x for x in chrInFile1 if x not in chrInFile2]
	inFile2Not1 = [x for x in chrInFile2 if x not in chrInFile1]
	if len(inFile1Not2) > 0 and len(inFile2Not1) > 0:
		chrInFile1_str = ",".join(chrInFile1)
		chrInFile2_str = ",".join(chrInFile2)
		err_msg = (("The chromosomes columns in %s do not match the " \
			+ "chromosomes in %s.\nThe chromosomes in %s are %s.\n" \
			+ "The chromosome in %s are %s") % (bedfile1, \
				bedfile2, bedfile1, chrInFile1_str, \
				bedfile2, chrInFile2_str))
		sys.exit(err_msg)

	bedTool1 = pybedtools.BedTool(bedfile1)
	bedTool2 = pybedtools.BedTool(bedfile2)
	intersect_bedTool = bedTool1.intersect(bedTool2, wo=True)

	colsList_for_file1 = ["chr", "start", "stop", "name", "signal", "strand", \
		"fold_change", "pValue", "qValue", "summit"]
	colsList_for_file2 = ["chr_b", "start_b", "stop_b", "name_b", "signal_b", \
		"strand_b", "fold_change_b", "pValue_b", "qValue_b", "summit_b"]
	final_col_list=[]
	final_col_list = \
		colsList_for_file1[0:ncols_file1] + \
		colsList_for_file2[0:ncols_file2] + ['overlap']
	overlap_df = intersect_bedTool.to_dataframe(
				names=final_col_list)
	overlap_df = overlap_df.loc[overlap_df['overlap'] >= distance, :]

	# report number of overlap
	if verbal:
		numOverlap = len(
			overlap_df.loc[:,["chr","start","stop"]].drop_duplicates())
		sys.stdout.write("NUMBER OF PEAKS THAT OVERLAP (%s, %s): %i\n" % \
			(bedfile1, bedfile2, numOverlap))
	
	return(overlap_df)



