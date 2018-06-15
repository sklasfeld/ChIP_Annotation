#!/usr/bin/env python

import os
import sys
import pandas as pd
import numpy as np
import subprocess


def cmd(cmd_str, speak):
	if speak:
		sys.stdout.write("%s\n" % cmd_str)
		sys.stdout.flush()
	os.system(cmd_str)

def peakoverlap(file_name):
	"""calculate peak overlap from bed file"""
	cmd = ("cut -f1,2,3 %s| sort -u| wc -l" % file_name)
	ps = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
		stderr=subprocess.STDOUT)
	output = ps.communicate()[0]
	break_out = output.split()
	return int(break_out[0])

def compare_bedfiles(peakfile1, peakfile2, outfile, *positional_parameters, \
	**keyword_parameters):
	"""compare two bedfiles to see if they overlap

	mandatory parameters:
	* peakfile1: bed file 
	* peakfile2: bed file
	* outfile: tab-delimited file containing regions that contain a 
	minimum $distance of overlap. The last column is the number of 
	bp that overlap.
	
	optional parameters:
	* distance: Minimum number of bp must overlap (Default:2)
	* verbal: print number of peaks that overlap (Default: False)"""

	distance=2
	verbal=False
	bedtools_path = ""

	if ('distance' in keyword_parameters):
		distance = keyword_parameters['distance']
	if ('verbal' in keyword_parameters):
		verbal = keyword_parameters['verbal']
	if ('bedtools_path' in keyword_parameters):
		bedtools_path = keyword_parameters['bedtools_path']

	file1_df = pd.read_csv(peakfile1,sep='\t', header=None, dtype=str)
	file2_df = pd.read_csv(peakfile2,sep='\t', header=None, dtype=str)
	chrInFile1 = list(file1_df.iloc[:,0].unique())
	chrInFile2 = list(file2_df.iloc[:,0].unique())
	inFile1Not2 = [x for x in chrInFile1 if x not in chrInFile2]
	if len(inFile1Not2) > 0:
		chrInFile1_str = ",".join(chrInFile1)
		chrInFile2_str = ",".join(chrInFile2)
		err_msg = (("The chromosomes columns in %s do not match the " \
			+ "chromosomes in %s.\nThe chromosomes in %s are %s.\n" \
			+ "The chromosome in %s are %s") % (peakfile1, \
				peakfile2, peakfile1, chrInFile1_str, \
				peakfile2, chrInFile2_str))
		sys.exit(err_msg)

	tmp_intersect1="intersectbed1.tmp"
	cmd1=("%sbedtools intersect -a %s -b %s -wo > %s" % \
		(bedtools_path, peakfile1, peakfile2, tmp_intersect1))
	cmd(cmd1, verbal)

	# remove lines below $distance
	cmd2=("""awk -F"\\t" '$NF>=%i' %s > %s""" % \
		(distance, tmp_intersect1, outfile))
	cmd(cmd2, verbal)

	# move tmp file
	cmd3=("rm %s" % tmp_intersect1)
	cmd(cmd3, verbal)

	if verbal:
		# report number of overlap
		sys.stdout.write("NUMBER OF PEAKS THAT OVERLAP: %s\n" % \
			peakoverlap(outfile))

	overlap_df = pd.read_csv(outfile, sep='\t', header=None, \
		names = ["chr", "start", "stop", "name", "signal", "strand", \
		"chr_b","start_b","stop_b","name_b","signal_b","strand_b", \
		"overlap"], \
		dtype={"chr" : object, "start" : np.int64, \
        "stop" : np.int64, "name" : object, "signal" : np.float64, \
        "strand":object, "chr_b" : object, "start_b" : np.int64, \
        "stop_b" : np.int64, "name_b" : object, "signal_b" : np.float64, \
        "strand_b":object,"overlap":np.int64})
	return(overlap_df)



