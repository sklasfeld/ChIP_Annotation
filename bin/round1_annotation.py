#!/usr/bin/env python

import os
import sys
from collections import defaultdict
import pandas as pd
import numpy as np
import math
import upstream_peaks
import subprocess

# print and run bash command
def cmd(cmd_str, verbose):
	if verbose:
		sys.stdout.write("%s\n" % cmd_str)
		sys.stdout.flush()
	os.system(cmd_str)

def wc(file_name, verbose=False, samtools_path=""):
    """count the number of lines in a file or reads in a sam/bam file"""
    if (not os.path.isfile(file_name)):
        sys.stdout.write("Count: 0\n")
        sys.stdout.flush()
        return (0)
    else:
        cmd = ("wc -l %s" % file_name)
    if os.path.splitext(file_name)[1] == ".gz":
        cmd = ("gunzip -c %s | wc -l" % file_name)
    elif os.path.splitext(file_name)[1] == ".bam" or \
                    os.path.splitext(file_name)[1] == ".sam" or \
                    os.path.splitext(file_name)[1] == ".tagAlign":
        cmd = ("%ssamtools view -F 4 -c %s" % (samtools_path, file_name))
    if verbose:
        sys.stdout.write("%s\n" % cmd)
        sys.stdout.flush()
    ps = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT)
    output = ps.communicate()[0]
    break_out = output.split()
    if verbose:
        sys.stdout.write("Count: %i\n" % int(break_out[0]))
        sys.stdout.flush()
    return int(break_out[0])


def r1_annotate(gene_alist, geneBed_file, bed_fname, peaks_df, prefix, \
				dir_name, *positional_parameters, **keyword_parameters):
	"""This function annotates genes based on the following preference order:
	 1) Inside Genes
	 2)Upstream of TSS
	 3) Downstream of TSS. 
	 This is round 1 of annotation. It also output peaks that are not annotated
	 in round 1 to use for round 2 annotation (which involve RNA-sequencing 
	 results.)
	 mandatory parameters:
	 * gene_alist - A file that contains gene IDs and their aliases.
	  The file must be tab-delimited and have atleast 1 column with the 
	  labels "gene_id"
	 the labe
	 * geneBed_file - bed file with gene locations
	 * peak_file - dataframe with all peak info
	 * prefix - prefix for experiment
	 * dir_name - output directory
	 optional parameters:
	 * gene_alist_cols - columns in `gene_alist` used to annotate the genes 
	 (eg. gene_name)
	 * bp_upstream_filter - maximum bp distance upstream of tss that peak can 
	 be annotated to a gene in round 1 (default:1000)
	 * bp_downstream_filter - maximum bp distance downstream of tts that peak can
	 be annotated to a gene in round 1 (default:100)
	 * ignore_conv_peaks - do not output peak information of peaks that annotate
	  to two different peaks (default: False)
	 * bedtools_path - path to bedtools module
	 * verbose - print out counts (default: False)
	"""
	# default keyword parameters
	gene_alist_cols = []
	bp_upstream_filter = 3000
	bp_downstream_filter = 0
	ignore_conv_peaks = False
	bedtools_path=""
	verbose= False
	# user-set keyword parameters
	if ('gene_alist_cols' in keyword_parameters):
		gene_alist_cols = keyword_parameters['gene_alist_cols']
	if ('bp_upstream_filter' in keyword_parameters):
		bp_upstream_filter = keyword_parameters['bp_upstream_filter']
	if ('bp_downstream_filer' in keyword_parameters):
		bp_downstream_filter = keyword_parameters['bp_downstream_filter']
	if ('ignore_conv_peaks' in keyword_parameters):
		ignore_conv_peaks = keyword_parameters['ignore_conv_peaks']
	if ('bedtools_path' in keyword_parameters):
		if len(keyword_parameters['bedtools_path']) > 0:
			bedtools_path = keyword_parameters['bedtools_path']
	if ('verbose' in keyword_parameters):
		verbose = keyword_parameters['verbose']
	narrowPeak_boolean = False

	
	if 'qValue' in peaks_df.columns:
		narrowPeak_boolean = True
	# create a dataframe with gene alias information
	geneAnn_df = pd.read_csv(gene_alist, sep='\t', dtype=str)
	if len(gene_alist_cols)>0:
		for gene_ann_col in gene_alist_cols:
			if not gene_ann_col in list(geneAnn_df.columns):
				sys.exit("%s is missing a column labeled %s" % (gene_alist, gene_ann_col))
	# STEP 1: Collect Peaks that are INTRAGENIC
	if verbose:
		sys.stdout.write("Annotating INTRA genic peaks...\n")
	# first get closest distance of all genes
	closestgenes_file = ("%s/%s_closest.txt" % (dir_name, prefix))
	closestgenes_cmd = ("%sbedtools closest -D b -a %s -b %s > %s" % \
		(bedtools_path, bed_fname, geneBed_file, closestgenes_file))
	cmd(closestgenes_cmd, verbose)
	if wc(closestgenes_file, verbose) == 0:
		sys.exit("\nAwkward!!! There are no peaks close to genes...\n"+ \
			"Failure to pass Round 1 Annotation has closed script...")
	# get bedtools closest output into df
	dis2genes_cols = list(peaks_df.columns[:6]) + \
		["gene_chr", "gene_chromStart", "gene_chromEnd", "gene_id", \
		"gene_score", "gene_strand"] + ["distance_from_gene"]
	closest_df = pd.read_csv(closestgenes_file, sep="\t", header=None, \
		names=dis2genes_cols, index_col=False, dtype={"chr" : object, \
			"start" : np.int64, "stop" : np.int64, "name" : object, \
			"signal" : np.float64, "strand":object, "gene_chr" : object, \
			"gene_chromStart" : np.int64, "gene_chromEnd" : np.int64, \
			"gene_id" : object, "gene_score" : np.float64, \
			"gene_strand":object, "distance_from_gene": np.float64})
	# ignore different types of the same CDS and limit to min_distance=0kb (intragenic)
	intragenic_peaks_df = closest_df.loc[closest_df["distance_from_gene"]==0.0,:].drop_duplicates()
	# add annotations found in peaks_df to intragenic peaks
	intragenic_peaks_df = intragenic_peaks_df.merge(peaks_df, how = 'left', on=list(peaks_df.columns[:6]))
	intragenic_peak_list =  list(intragenic_peaks_df['name'].unique())
	# STEP 2: Collect Peaks that are UPSTREAM of genes
	if verbose:
		sys.stdout.write("Annotating peaks UPSTREAM of genes ...\n")
	# get dataframe of all INTERGENIC peaks (not found inside genes)
	intergenic_peaks_df = peaks_df[~peaks_df["name"].isin(intragenic_peaks_df["name"])]
	# first ignore peaks downstream (-id) of genes since upstream is preferred
	upsteamgenes_file = ("%s/%s_closest_id.txt" % (dir_name, prefix))
	upstream_peaks_df = pd.DataFrame()
	#if not os.path.isfile(upsteamgenes_file):
	if ignore_conv_peaks:
		upstream_peaks_df = upstream_peaks.annotate(prefix, \
			intergenic_peaks_df, geneBed_file, bp_upstream_filter, \
			dir_name, ignore_conv_peaks=ignore_conv_peaks)
	else:
		upstream_peaks_df = upstream_peaks.annotate(prefix, \
			intergenic_peaks_df, geneBed_file, bp_upstream_filter, \
			dir_name)
	#else:
	#	sys.stderr.write("The file %s already exists so I am not going to overwrite it\n" \
	#		% upsteamgenes_file)
	#	# get upstream distance output into df
	#	upstream_peaks_df = pd.read_csv(upsteamgenes_file, sep="\t", header=None, \
	#		names=["chr", "start", "stop", "name", "signal", "strand", \
	#		"tss_chr", "tss_start", "tss_stop", "gene_id", "tss_score", \
	#		"tss_strand", "distance_from_gene"], \
	#		index_col=False, dtype={"chr" : object, \
	#		"start" : np.int64, "stop" : np.int64, "name" : object, \
	#		"signal" : np.float64, "strand":object, "tss_chr" : object, \
	#		"tss_start" : np.int64, "tss_stop" : np.int64, \
	#		"gene_id" : object, "tss_score" : np.float64, \
	#		"tss_strand":object, "distance_from_gene": np.float64})
	
	upstream_peaks_df = upstream_peaks_df.merge(peaks_df, how = 'left', on=list(peaks_df.columns[:6]))
	# STEP 3: Collect Peaks that are DOWNSTREAM of genes
	if verbose:
		sys.stdout.write("Annotating peaks DOWNSTREAM of genes ...\n")
	downstream_peaks_df = pd.DataFrame()
	if bp_downstream_filter != 0:
		# get dataframe of all peaks that have not already been annotated
		# as intragenic nor upstream of a gene
		intraButNotUpstream_df = intergenic_peaks_df[~intergenic_peaks_df["name"].isin(upstream_peaks_df["name"])]	
		# file for downstream peaks
		downstreamgenes_file = ("%s/%s_closest_iu.txt" % (dir_name, prefix))
		if ignore_conv_peaks:
			downstream_peaks_df = downstream_peaks.annotate(prefix, \
				intraButNotUpstream_df, geneBed_file, \
				bp_downstream_filter, dir_name, \
				ignore_conv_peaks=ignore_conv_peaks)
		else:
			downstream_peaks_df = downstream_peaks.annotate(prefix, \
				intraButNotUpstream_df, geneBed_file, \
				bp_downstream_filter, dir_name)
	if len(downstream_peaks_df) > 0:
		downstream_peaks_df = downstream_peaks_df.merge(peaks_df, how = 'left', on=list(peaks_df.columns[:6]))
	# connect upstream and downstream peaks to gene locations rather than
	# tss or tts respectively
	if verbose:
		sys.stdout.write("Collecting info about annotated peaks ...\n")
	peak2gene_cols = list(intragenic_peaks_df.columns)
	gene_table = pd.read_csv(geneBed_file, sep='\t', header=None)
	gene_table.columns=peak2gene_cols[6:12]
	# remove tss info from upstream dataframe
	upstream_peaks_df = upstream_peaks_df.drop(['tss_chr', 'tss_start', \
		'tss_stop', 'tss_score','tss_strand'],axis=1)

	# add gene info to upstream dataframe
	round1_frame = []
	if len(upstream_peaks_df) > 0:
		upstream_peaks_2_genes = upstream_peaks_df.merge(gene_table, \
			how = 'left', on='gene_id')
		if len(downstream_peaks_df) > 0:
			# remove tts info from downstream dataframe
			downstream_peaks_df = downstream_peaks_df.drop(['tts_chr', 'tts_start', \
				'tts_stop', 'tts_score','tts_strand'],axis=1)
			# add gene info to upstream dataframe
			downstream_peaks_2_genes = downstream_peaks_df.merge(gene_table, \
				how = 'left', on='gene_id')
			#downstream_peaks_2_genes = downstream_peaks_2_genes.loc[:,peak2gene_cols]
			# concatenate the annotated peaks into one dataframe
			if len(intragenic_peaks_df) > 0:
				round1_frame = [intragenic_peaks_df, upstream_peaks_2_genes, downstream_peaks_2_genes]
			else:
				round1_frame = [upstream_peaks_2_genes, downstream_peaks_2_genes]
		else:
			if len(intragenic_peaks_df) > 0:
				round1_frame = [intragenic_peaks_df, upstream_peaks_2_genes]
			else:
				round1_frame = [upstream_peaks_2_genes]
	else:
		if len(downstream_peaks_df) > 0:
			# remove tts info from downstream dataframe
			downstream_peaks_df = downstream_peaks_df.drop(['tts_chr', 'tts_start', \
				'tts_stop', 'tts_score','tts_strand'],axis=1)
			# add gene info to upstream dataframe
			downstream_peaks_2_genes = downstream_peaks_df.merge(gene_table, \
				how = 'left', on='gene_id')
			#downstream_peaks_2_genes = downstream_peaks_2_genes.loc[:,peak2gene_cols]
			# concatenate the annotated peaks into one dataframe
			if len(intragenic_peaks_df) > 0:
				round1_frame = [intragenic_peaks_df, downstream_peaks_2_genes]
			else:
				round1_frame = [downstream_peaks_2_genes]
		else:
			if len(intragenic_peaks_df) > 0:
				round1_frame = [intragenic_peaks_df]
			else:
				sys.exit("\nAwkward!!! There are no INTRAGENIC peaks annotated in round 1...\n"+ \
			"Failure to pass Round 1 Annotation has closed script...")
	#round1_df = pd.concat(round1_frame,sort=False)
	round1_df = pd.concat(round1_frame)
	# add gene Alias and gene description information

	round1_alias_df = round1_df.merge(geneAnn_df, how='left', on='gene_id')
	# sort the dataframe by peak location
	round1_alias_df["start"] = round1_alias_df["start"].astype('int64')
	round1_alias_df.sort_values(by=["chr","start"], axis=0, ascending=True, inplace=True)
	# print information for each peak
	round1_ann = ("%s/%s_r1_peak_annotations.txt" % (dir_name, prefix))
	peaks_group_cols = peak2gene_cols[0:6]
	## group by peak info
	peak_groups_df = round1_alias_df.groupby(peaks_group_cols)
	## peak centric columns
	peaks_centric_cols = []
	if narrowPeak_boolean:
		peaks_centric_cols = peaks_group_cols+["qValue"]+list(peaks_df.columns[10:])
	else:
		peaks_centric_cols += list(peaks_df.columns[6:])
	peak_ann_df = round1_alias_df.loc[:,peaks_group_cols]
	## put columns that are not peak centric into a peak context
	peak_nGenes_series = peak_groups_df.apply(lambda x: len(x["gene_id"].unique()))
	peak_nGenes_df = peak_nGenes_series.to_frame().reset_index()
	peak_nGenes_df.columns=peaks_group_cols+['numGenes']
	peak_ann_df = peak_ann_df.merge(peak_nGenes_df,how='outer',on=peaks_group_cols)
	peak_gid_series = peak_groups_df.apply(lambda x: ";".join(str(s) for s in list(x["gene_id"])))
	peak_gid_df = peak_gid_series.to_frame().reset_index()
	peak_gid_df.columns=peaks_group_cols+['gene_id']
	peak_ann_df = peak_ann_df.merge(peak_gid_df,how='outer',on=peaks_group_cols)
	if len(gene_alist_cols)>0:
		for gene_ann_col in gene_alist_cols:
			peak_gname_series = peak_groups_df.apply(lambda x: ";".join(str(s) for s in list(x[gene_ann_col])))
			peak_gname_df = peak_gname_series.to_frame().reset_index()
			peak_gname_df.columns=peaks_group_cols+[gene_ann_col]
			peak_ann_df = peak_ann_df.merge(peak_gname_df,how='outer',on=peaks_group_cols)
	peak2gene_info_cols = ['numGenes','gene_id'] + gene_alist_cols
	column_order = peak_ann_df.columns
	if narrowPeak_boolean:
		column_order= peak2gene_cols[0:6] + ["qValue"] + \
			peak2gene_info_cols+list(peaks_df.columns[10:])	
	else:
		column_order= peak2gene_cols[0:6] + \
			peak2gene_info_cols+list(peaks_df.columns[6:])
	peak_ann_df = peak_ann_df.reindex(columns=column_order)
	pd.set_option('float_format', '{:.2f}'.format)
	peak_ann_df.to_csv(round1_ann, sep="\t", index=False, na_rep="NA")
	# get gene-centric annotation
	gene_groups_df = round1_alias_df.groupby(peak2gene_cols[6:12])
	gene_nPeaks_series = gene_groups_df.apply(lambda x: x.shape[0])
	gene_nPeaks_df = gene_nPeaks_series.to_frame().reset_index()
	gene_nPeaks_df.columns=peak2gene_cols[6:12]+['numPeaks']
	gene_dis_series = gene_groups_df.apply(lambda x: ";".join(str(s) for s in list(x["distance_from_gene"])))
	gene_dis_df = gene_dis_series.to_frame().reset_index()
	gene_dis_df.columns=peak2gene_cols[6:12]+['distance']
	gene_df = gene_nPeaks_df.merge(gene_dis_df,how='outer',on=peak2gene_cols[6:12])
	if narrowPeak_boolean:
		gene_qval_series = gene_groups_df.apply(lambda x: ";".join(str(s) for s in list(x["qValue"])))
		gene_qval_df = gene_qval_series.to_frame().reset_index()
		gene_qval_df.columns=peak2gene_cols[6:12]+['qValue']
		gene_df = gene_df.merge(gene_qval_df,how='outer',on=peak2gene_cols[6:12])
	gene_pname_series = gene_groups_df.apply(lambda x: ";".join(str(s) for s in list(x["name"])))
	gene_pname_df = gene_pname_series.to_frame().reset_index()
	gene_pname_df.columns=peak2gene_cols[6:12]+['peak_name']
	gene_df = gene_df.merge(gene_pname_df,how='outer',on=peak2gene_cols[6:12])
	round1_gene_ann = ("%s/%s_r1_gene_annotations.txt" % (dir_name, prefix))
	gene_df.to_csv(round1_gene_ann, sep="\t", index=False)
	if verbose:
		sys.stdout.write("Number of Unique Genes that are Annotated: %i\n" % \
			gene_df.shape[0])
	return (round1_alias_df)








