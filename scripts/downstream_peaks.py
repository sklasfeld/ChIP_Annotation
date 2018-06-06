#!/usr/bin/env python

import os
import sys
import subprocess
import pandas as pd
import numpy as np

def annotate(prefix, intraButNotUpstream_df, gene_file, maximum_bp, outfile, \
	*positional_parameters, **keyword_parameters):
	"""This function annotates downstream peaks to the closest TTS 
	within maximum number of bp. Note that if a TTS is within range but on the same 
	strand as a closer TTS than we ignore the further TTS.

	mandatory parameters:
	* prefix - prefix for experiment
	* intraButNotUpstream_df - data table that has all peaks that are not within 
	genes nor annotated as upstream of a nearby gene
	* gene_file - bed file containing tts sites.
	* maximum_bp - maximum bp upstream
	* outfile - text file output
	optional parameters:
	* ignore_conv_peaks - remove convergent peak information (peaks that annotate 
	to two different peaks) and put their info in a seperate file (default:False)"""
	
	# default keyword parameters
	ignore_conv_peaks = False

	# user-set keyword parameters
	if ('ignore_conv_peaks' in keyword_parameters):
		ignore_conv_peaks = keyword_parameters['ignore_conv_peaks']

	

	# get tts_table from gene table
	tts_table = pd.read_csv(geneBed_file, sep='\t', header=None)
	tts_table.columns = ['tts_chr', 'tts_start', 'tts_stop', \
		'gene_id', 'tts_score', 'tts_strand']
	tts_table.loc[tts_table["tts_strand"] == "+", "tts_start"] = \
		tts_table.loc[tts_table["tts_strand"] == "+", "tts_stop"] - 1
	tts_table.loc[tts_table["tts_strand"] == "-", "tts_stop"] = \
		tts_table.loc[tts_table["tts_strand"] == "-", "tts_start"] + 1
	
	


	# get table with joined info on tts and peaks
	peaks2genes_df = intraButNotUpstream_df.merge(tts_table, how = "inner", right_on="tts_chr", left_on="chr")
	peaks2genes_df["start"] = peaks2genes_df["start"].astype('int64')
	peaks2genes_df["stop"] = peaks2genes_df["stop"].astype('int64')
	peaks2genes_df["tts_start"] = peaks2genes_df["tts_start"].astype('int64')
	peaks2genes_df["tts_stop"] = peaks2genes_df["tts_stop"].astype('int64')

	# Only get rows that are upstream of genes

	# negative strand downstream 
	neg_downstream_df = peaks2genes_df[((peaks2genes_df["tts_strand"]=="-") & \
		(peaks2genes_df["stop"] <= peaks2genes_df["tts_start"]))]
	# positive strand downstream
	pos_downstream_df = peaks2genes_df[((peaks2genes_df["tts_strand"]=="+") & \
		(peaks2genes_df["start"] >= peaks2genes_df["tts_stop"]))]

	# calculate distance for each group
	neg_downstream_df.loc[:,"distance_from_gene"] = (neg_downstream_df.loc[:,"tts_start"] - \
		(neg_downstream_df.loc[:,"stop"]-1))
	pos_downstream_df.loc[:,"distance_from_gene"] = (pos_downstream_df.loc[:,"start"] - \
		(pos_downstream_df.loc[:,"tts_stop"]-1))



	# Concatonate the groups together now that their distance is calculated
	downstream_df = neg_downstream_df.append(pos_downstream_df)

	downstream_df = downstream_df[downstream_df["distance_from_gene"]<maximum_bp]

	# group by peak name
	# * if the peak has only one annotation then keep the annotation
	# * if the peak has 2+ annotations in same direction then 
	# annotate the peak to the closest one
	# * if the peak has 2+ annotations in different directions than annotate
	# the peak to the 2 closest peaks if NOT ignoring convergent peaks. Otherwise,
	# do not annotate this peak
	downstream_filtered_arr = []
	convergent_peaks_arr = []
	downstream_pGroups = downstream_df.groupby(list(intraButNotUpstream_df.columns))
	for p_info, p_group in downstream_pGroups:
		if p_group.shape[0]==1:
			downstream_filtered_arr.append(p_group)
		else:
			if len(p_group["tts_strand"].unique()) == 1:
				closest_gene_series = p_group.loc[p_group["distance_from_gene"].argmin(),:]
				closest_gene_df = pd.DataFrame(columns=downstream_df.columns)
				closest_gene_df = closest_gene_df.append(closest_gene_series)
				downstream_filtered_arr.append(closest_gene_df)
			else:
				# seperate the peaks downstream of positive strand genes
				# and downstream of negative strand genes
				neg_downstream_peak_df = p_group[p_group["tts_strand"]=="-"]
				pos_downstream_peak_df = p_group[p_group["tts_strand"]=="+"]
				# get the peaks closest to an downstream gene
				closest_gene_neg_series = \
					neg_downstream_peak_df.loc[neg_downstream_peak_df["distance_from_gene"].argmin(),:]
				closest_gene_pos_series = \
					pos_downstream_peak_df.loc[pos_downstream_peak_df["distance_from_gene"].argmin(),:]
				closest_gene_df = pd.DataFrame(columns=downstream_df.columns)
				closest_gene_df = closest_gene_df.append(closest_gene_neg_series)
				closest_gene_df = closest_gene_df.append(closest_gene_pos_series)
				convergent_peaks_arr.append(closest_gene_df)
				if ~ignore_conv_peaks:
					downstream_filtered_arr.append(closest_gene_df)

	downstream_filtered_df = pd.DataFrame(columns=downstream_df.columns)
	if len(downstream_filtered_arr) > 0:
		downstream_filtered_df = pd.concat(downstream_filtered_arr)
	downstream_filtered_df = downstream_filtered_df[list(intraButNotUpstream_df.columns[0:6])+list(tts_table.columns)+["distance_from_gene"]]
	downstream_filtered_df.to_csv(outfile, sep="\t", index=False, header=False)

	if len(convergent_peaks_arr) > 0:
		convergent_peaks_df = pd.concat(convergent_peaks_arr)
		convergent_peaks_file = ("%s_convergent_downstream_peaks.txt" % prefix)
		convergent_peaks_df = convergent_peaks_df[list(intraButNotUpstream_df.columns[0:6])+list(tts_table.columns)+["distance_from_gene"]]
		convergent_peaks_df.to_csv(convergent_peaks_file, sep="\t", index=False, header=True)

	return(downstream_filtered_df)


