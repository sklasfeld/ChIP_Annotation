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

	overlap_df = pd.read_csv(outfile, sep='\t', header=None, dtype=str)
	return(overlap_df)



