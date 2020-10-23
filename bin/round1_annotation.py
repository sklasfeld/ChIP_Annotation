#!/usr/bin/env python

import os
import sys
import pandas as pd
import numpy as np
import bin.upstream_peaks
import pybedtools


def r1_annotate(by, geneBed_file, bed_fname, peaks_df, prefix, \
				dir_name, *positional_parameters, **keyword_parameters):
	"""This function annotates genes based on the following preference order:
	 1) Inside Genes
	 2)Upstream of TSS
	 3) Downstream of TSS. 
	 This is round 1 of annotation. It also output peaks that are not annotated
	 in round 1 to use for round 2 annotation (which involve RNA-sequencing 
	 results.)
	 mandatory parameters:
	 * by - annotate by "peak" or "summit"
	 * geneBed_file - bed file with gene locations
	 * name of bed file if annotating by peaks and narrowPeak if annoting by summits
	 * peak_file - dataframe with all peak info
	 * prefix - prefix for experiment
	 * dir_name - output directory
	 optional parameters:
	 * per_inter_filter - annotates genes that overlap >= this percent of the gene
	 * bp_upstream_filter - maximum bp distance upstream of tss that peak can 
	 be annotated to a gene in round 1 (default:1000)
	 * bp_downstream_filter - maximum bp distance downstream of tts that peak can
	 be annotated to a gene in round 1 (default:100)
	 * ignore_conv_peaks - do not output peak information of peaks that annotate
	  to two different peaks (default: False)
	 * verbose - print out counts (default: False)
	 * round2_bool - if False and nothing is found to annotate in round one
	  then the script must exit.
	"""
	# default keyword parameters
	per_inter_filter = 0
	bp_upstream_filter = 3000
	bp_downstream_filter = 0
	ignore_conv_peaks = False
	verbose= False
	keep_tmps = False
	narrowPeak_boolean = False
	region_cols = list(peaks_df.columns[:6])
	round2_bool = False


	# user-set keyword parameters
	if ('per_inter_filter' in keyword_parameters):
		per_inter_filter = keyword_parameters['per_inter_filter']
	if ('bp_upstream_filter' in keyword_parameters):
		bp_upstream_filter = keyword_parameters['bp_upstream_filter']
	if ('bp_downstream_filer' in keyword_parameters):
		bp_downstream_filter = keyword_parameters['bp_downstream_filter']
	if ('ignore_conv_peaks' in keyword_parameters):
		ignore_conv_peaks = keyword_parameters['ignore_conv_peaks']
	if ('verbose' in keyword_parameters):
		verbose = keyword_parameters['verbose']
	if ('keep_tmps' in keyword_parameters):
		keep_tmps = keyword_parameters['keep_tmps']
	if ('round2_bool' in keyword_parameters):
		round2_bool = keyword_parameters['round2_bool']
	if 'qValue' in peaks_df.columns:
		narrowPeak_boolean = True
	if by == "summit":
		original_np_file = bed_fname
		narrowPeak_df = peaks_df.iloc[:,range(0,10)].copy()
		narrowPeak_df['summit_start'] = narrowPeak_df['start'] + \
			narrowPeak_df['summit']
		narrowPeak_df['summit_stop'] = narrowPeak_df['summit_start'] + 1
		summit_bed = narrowPeak_df.copy()
		summit_bed['start'] = summit_bed['summit_start']
		summit_bed['stop'] = summit_bed['summit_stop']
		summit_bed = summit_bed.loc[:,region_cols]
		temp_summit_bed_file = (("%s/%s_summits.tmp.bed") % (dir_name,prefix))
		summit_bed.to_csv(temp_summit_bed_file, sep="\t", header=False, \
			index=False)
		bed_fname = temp_summit_bed_file
		# get dataframe for table in file  `bed_fname`
		region_df = summit_bed
	else:
		# get dataframe for table in file  `bed_fname`
		region_df = peaks_df


	# gene bed columns
	geneBed_cols = ["gene_chr", "gene_start", "gene_stop", "gene_id", \
		"gene_score", "gene_strand"]

	# create a BedTool for the gene Bed file
	geneBedTool = pybedtools.BedTool(geneBed_file)
	# create a BedTools for the chip Bed/NarrowPeak file
	chipBedTool = pybedtools.BedTool(bed_fname)
	
	
	# STEP 0: first get closest distance of all genes
	# this is simply to check that the peaks
	# are close to any genes
	if bp_upstream_filter > 0 or bp_downstream_filter > 0:
		closestGenes2Peak_BedTool = chipBedTool.closest(geneBedTool,
				s=False, D="b")
		if len(closestGenes2Peak_BedTool) == 0:
			sys.exit("\nAwkward!!! There are no peaks close to genes...\n"+ \
				"Failure to pass Round 1 Annotation has closed script...")

		if keep_tmps:
			closestgenes_file = ("%s/%s_closest_by%s.tsv" % (dir_name, prefix, by))
			closestGenes2Peak_df = closestGenes2Peak_BedTool.to_dataframe(
				names=region_cols+geneBed_cols+["distance"])
			closestGenes2Peak_df = closestGenes2Peak_df["signal"].astype(np.float64)
			closestGenes2Peak_df.to_csv(closestgenes_file, 
				sep="\t", index=False)

	# STEP 1: Collect Peaks that are INTRAGENIC
	if verbose:
		sys.stdout.write("Annotating INTRA genic peaks...\n")

	bedtools_intersect_param={}
	if per_inter_filter > 0:
		bedtools_intersect_param["F"]=per_inter_filter
	bedtoolsIntragenicBedTool = chipBedTool.intersect(
		geneBedTool, wo=True, **bedtools_intersect_param)
	genesOverlap_cols =region_cols + geneBed_cols + ["gene_overlap"]
	if len(bedtoolsIntragenicBedTool) == 0:
		if(bp_upstream_filter == 0 and bp_downstream_filter == 0):
			sys.exit("\nAwkward!!! There are no peaks in genes...\n"+ 
			"Failure to pass Round 1 Annotation has closed script...")
		else:
			sys.stderr.write(("Warning: There are NO INTRAGENIC peaks " + 
				"(within a gene).\n"))
			intragenic_peaks_df=pd.DataFrame(
				columns=genesOverlap_cols+["distance_from_gene"])
	else:

		intragenic_peaks_df = bedtoolsIntragenicBedTool.to_dataframe(
				names=genesOverlap_cols)
		intragenic_peaks_df["signal"] = \
			intragenic_peaks_df["signal"].astype(np.float64)
		if keep_tmps:
			bedtools_intragenic_file = ("%s/%s_intragenic_by%s.txt" % 
				(dir_name, prefix, by))
			intragenic_peaks_df.to_csv(bedtools_intragenic_file, 
				sep="\t", index=False)

		intragenic_peaks_df.loc[:,"gene_overlap"] = (
			intragenic_peaks_df.loc[:,"gene_overlap"] /
			np.float64(intragenic_peaks_df.loc[:,"gene_stop"] - 
				intragenic_peaks_df.loc[:,"gene_start"]))
		intragenic_peaks_df.loc[:,"distance_from_gene"]=0
		if by=="summit":
			orig_col_order = list(intragenic_peaks_df.columns)
			intragenic_peaks_df = intragenic_peaks_df.merge(narrowPeak_df, 
				how = 'left', left_on=region_cols,
				right_on=['chr', 'summit_start', 'summit_stop', 'name', 'signal', 'strand'])
			intragenic_peaks_df = intragenic_peaks_df.drop(
				['start_x','stop_x','fold_change', 'pValue', \
				'qValue', 'summit', 'summit_stop'], axis=1)
			intragenic_peaks_df = intragenic_peaks_df.rename(
				columns = {"start_y":"start", "stop_y":"stop"})
			intragenic_peaks_df = intragenic_peaks_df.loc[:,orig_col_order+['summit_start']]
		print("BEEEEE4444")
		print(peaks_df.loc[:,region_cols+["highMock_MNase2"]].head())
		print(intragenic_peaks_df.loc[:,region_cols+["highMock_MNase2"]].head())
		#print(peaks_df.columns)
		#print(intragenic_peaks_df.columns)
		intragenic_peaks_df = intragenic_peaks_df.merge(peaks_df, how = 'left', 
			on=region_cols,sort=True)
		print("aftaaa")
		print(peaks_df.loc[:,region_cols+["highMock_MNase2"]].head())
		print(intragenic_peaks_df.loc[:,region_cols+["highMock_MNase2"]].head())
		
		#intragenic_peak_list =  list(intragenic_peaks_df['name'].unique())
	

	# STEP 2: Collect Peaks that are UPSTREAM of genes
	if bp_upstream_filter > 0:
		if verbose:
			sys.stdout.write("Annotating peaks UPSTREAM of genes ...\n")
		# get dataframe of all INTERGENIC peaks (not found inside genes)
		if len(intragenic_peaks_df) > 0:
			if by=="peak":
				intergenic_peaks_df = region_df.merge(
					intragenic_peaks_df, how='outer', indicator=True,
					on=list(region_df.columns))

				intergenic_peaks_df = intergenic_peaks_df.loc[ \
					intergenic_peaks_df["_merge"]=="left_only", :]
				intergenic_peaks_df = intergenic_peaks_df.loc[:, list(region_df.columns)]
			else:
				intergenic_peaks_df = region_df.merge(intragenic_peaks_df,
					left_on=['chr', 'start', 'name', 'signal', 'strand'],
					right_on=['chr', 'summit_start', 'name', 'signal', 'strand'],
					how='outer', indicator=True)
				intergenic_peaks_df = intergenic_peaks_df.loc[ \
					intergenic_peaks_df["_merge"]=="left_only", :]
				intergenic_peaks_df = intergenic_peaks_df.rename(
					columns = {"start_x":"start", "stop_x":"stop"})
				intergenic_peaks_df = intergenic_peaks_df.loc[:, list(region_df.columns)]
		else:
			intergenic_peaks_df = region_df.copy()
		# first ignore peaks downstream (-id) of genes since upstream is preferred
		upstream_peaks_df = pd.DataFrame()
		if ignore_conv_peaks:
			upstream_peaks_df = bin.upstream_peaks.annotate(prefix, \
				intergenic_peaks_df, geneBed_file, bp_upstream_filter, \
				dir_name, ignore_conv_peaks=ignore_conv_peaks)
		else:
			upstream_peaks_df = bin.upstream_peaks.annotate(prefix, \
				intergenic_peaks_df, geneBed_file, bp_upstream_filter, \
				dir_name)
		if by=="summit":
			orig_col_order = list(upstream_peaks_df.columns)
			upstream_peaks_df = upstream_peaks_df.merge(narrowPeak_df, 
				how = 'left', left_on=region_cols,
				right_on=['chr', 'summit_start', 'summit_stop', 'name', 'signal', 'strand'])
			upstream_peaks_df = upstream_peaks_df.drop(
				['start_x','stop_x','fold_change', 'pValue', \
				'qValue', 'summit', 'summit_stop'], axis=1)
			upstream_peaks_df = upstream_peaks_df.rename(
				columns = {"start_y":"start", "stop_y":"stop"})
			upstream_peaks_df = upstream_peaks_df.loc[:,orig_col_order+['summit_start']]
		if len(upstream_peaks_df) > 0:
			upstream_peaks_df = upstream_peaks_df.merge(peaks_df, how = 'left', on=list(peaks_df.columns[:6]))
			upstream_peaks_df.loc[:,"gene_overlap"]=0
		# STEP 3: Collect Peaks that are DOWNSTREAM of genes
	if verbose:
		sys.stdout.write("Annotating peaks DOWNSTREAM of genes ...\n")
	downstream_peaks_df = pd.DataFrame()
	if bp_downstream_filter != 0:
		# get dataframe of all peaks that have not already been annotated
		# as intragenic nor upstream of a gene
		if by=="peak":
			intraButNotUpstream_df = intergenic_peaks_df.merge(
				upstream_peaks_df, how='outer', indicator=True,
				on=list(intergenic_peaks_df.columns))
			intraButNotUpstream_df = intraButNotUpstream_df.loc[ \
				intraButNotUpstream_df["_merge"]=="left_only", :]
			intraButNotUpstream_df = intraButNotUpstream_df.loc[:, list(intergenic_peaks_df.columns)]
			#intergenic_peaks_df = intergenic_peaks_df.drop(["_merge"], axis=1)
		else:
			intraButNotUpstream_df = intergenic_peaks_df.merge(intraButNotUpstream_df,
				left_on=['chr', 'start', 'name', 'signal', 'strand'],
				right_on=['chr', 'summit_start', 'name', 'signal', 'strand'],
				how='outer', indicator=True)
			intraButNotUpstream_df = intraButNotUpstream_df.loc[ \
				intraButNotUpstream_df["_merge"]=="left_only", :]
			intraButNotUpstream_df = intraButNotUpstream_df.rename(
				columns = {"start_x":"start", "stop_x":"stop"})
			intraButNotUpstream_df = intraButNotUpstream_df.loc[:, list(intergenic_peaks_df.columns)]	
		# file for downstream peaks
		downstreamgenes_file = ("%s/%s_closest_by%s_iu.txt" % (dir_name, prefix,by))
		if ignore_conv_peaks:
			downstream_peaks_df = downstream_peaks.annotate(prefix, \
				intraButNotUpstream_df, geneBed_file, \
				bp_downstream_filter, dir_name, \
				ignore_conv_peaks=ignore_conv_peaks)
		else:
			downstream_peaks_df = downstream_peaks.annotate(prefix, \
				intraButNotUpstream_df, geneBed_file, \
				bp_downstream_filter, dir_name)
		if by=="summit":
			orig_col_order = list(downstream_peaks_df.columns)
			downstream_peaks_df = downstream_peaks_df.merge(narrowPeak_df, 
				how = 'left', left_on=region_cols,
				right_on=['chr', 'summit_start', 'summit_stop', 'name', 'signal', 'strand'])
			downstream_peaks_df = downstream_peaks_df.drop(
				['start_x','stop_x','fold_change', 'pValue', \
				'qValue', 'summit', 'summit_stop'], axis=1)
			downstream_peaks_df = downstream_peaks_df.rename(
				columns = {"start_y":"start", "stop_y":"stop"})
			downstream_peaks_df = downstream_peaks_df.loc[:,orig_col_order+['summit_start']]
	if len(downstream_peaks_df) > 0:
		downstream_peaks_df = downstream_peaks_df.merge(peaks_df, how = 'left', on=list(peaks_df.columns[:6]))
		downstream_peaks_df.loc[:,"gene_overlap"]=None
	# connect upstream and downstream peaks to gene locations rather than
	# tss or tts respectively
	if verbose:
		sys.stdout.write("Collecting info about annotated peaks ...\n")
	peak2gene_cols = list(intragenic_peaks_df.columns)
	gene_table = pd.read_csv(geneBed_file, sep='\t', header=None)
	gene_table.columns=peak2gene_cols[6:12]
	gene_bedfile_types={"gene_chr" : object, \
		"gene_start" : np.int64, "gene_stop" : np.int64, "gene_id" : object, \
		"gene_score" : np.float64, "gene_strand":object}
	if "." in list(gene_table["gene_score"]):
		gene_table.loc[gene_table["gene_score"]==".","gene_score"]=0
	gene_table = gene_table.astype(gene_bedfile_types)

	if bp_upstream_filter > 0:
		# remove tss info from upstream dataframe
		upstream_peaks_df = upstream_peaks_df.drop(['tss_chr', 'tss_start', \
			'tss_stop', 'tss_score','tss_strand'],axis=1)

	# add gene info to upstream dataframe
	round1_frame = []
	if bp_upstream_filter > 0 and len(upstream_peaks_df) > 0:
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
				if round2_bool:
					sys.stderr.write("\nWarning: There are NO INTERGENIC peaks " + 
						"annotated in round 1...\n")
					sys.stderr.flush()
				else:
					sys.exit("ERROR: no peaks could be annotated to genes.")
	if len(round1_frame) > 0:
		round1_df = pd.concat(round1_frame, sort=True)
	else:
		round1_df = intragenic_peaks_df
	# sort the dataframe by peak location
	round1_df["start"] = round1_df["start"].astype('int64')
	if len(round1_frame) > 0:
		round1_df.sort_values(by=["chr","start"], axis=0, ascending=True, inplace=True)
		# print information for each peak
		round1_ann = ("%s/%s_r1_peak_annotations.txt" % (dir_name, prefix))
		peaks_group_cols = peak2gene_cols[0:6]
		## group by peak info
		peak_groups_df = round1_df.groupby(peaks_group_cols)
		## peak centric columns
		peaks_centric_cols = []
		if narrowPeak_boolean:
			peaks_centric_cols = peaks_group_cols+["qValue"] + \
				list(peaks_df.columns[10:])
		else:
			peaks_centric_cols += list(peaks_df.columns[6:])
		peak_ann_df = round1_df.loc[:,peaks_group_cols]
		## put columns that are not peak centric into a peak context
		peak_nGenes_series = peak_groups_df.apply(
			lambda x: len(x["gene_id"].unique()))
		peak_nGenes_df = peak_nGenes_series.to_frame().reset_index()
		peak_nGenes_df.columns=peaks_group_cols+['numGenes']
		peak_ann_df = peak_ann_df.merge(
			peak_nGenes_df,how='outer',on=peaks_group_cols)
		peak_gid_series = peak_groups_df.apply(lambda x: ";".join(
			str(s) for s in list(x["gene_id"])))
		peak_gid_df = peak_gid_series.to_frame().reset_index()
		peak_gid_df.columns=peaks_group_cols+['gene_id']
		peak_ann_df = peak_ann_df.merge(peak_gid_df,how='outer',on=peaks_group_cols)
		peak2gene_info_cols = ['numGenes','gene_id']
	
		column_order = peak_ann_df.columns
		if narrowPeak_boolean:
			column_order= peak2gene_cols[0:6] + ["qValue"] + \
				peak2gene_info_cols+list(peaks_df.columns[10:])	
		else:
			column_order= peak2gene_cols[0:6] + \
				peak2gene_info_cols+list(peaks_df.columns[6:])
		peak_ann_df = peak_ann_df.reindex(columns=column_order)
		pd.set_option('float_format', '{:.2f}'.format)
		if len(round1_ann) > 0:
			peak_ann_df.to_csv(round1_ann, sep="\t", index=False, na_rep="NA")
		# get gene-centric annotation
		gene_groups_df = round1_df.groupby(peak2gene_cols[6:12])
		gene_nPeaks_series = gene_groups_df.apply(lambda x: x.shape[0])
		gene_nPeaks_df = gene_nPeaks_series.to_frame().reset_index()
		gene_nPeaks_df.columns=peak2gene_cols[6:12]+['numPeaks']
		gene_dis_series = gene_groups_df.apply(lambda x: ";".join(
			str(s) for s in list(x["distance_from_gene"])))
		gene_overlap_series = gene_groups_df.apply(lambda x: ";".join(
			str(s) for s in list(x["gene_overlap"])))
		gene_dis_df = gene_dis_series.to_frame().reset_index()
		gene_dis_df.columns=peak2gene_cols[6:12]+['distance']
		gene_df = gene_nPeaks_df.merge(gene_dis_df,how='outer',on=peak2gene_cols[6:12])
		if narrowPeak_boolean:
			gene_qval_series = gene_groups_df.apply(lambda x: ";".join(
				str(s) for s in list(x["qValue"])))
			gene_qval_df = gene_qval_series.to_frame().reset_index()
			gene_qval_df.columns=peak2gene_cols[6:12]+['qValue']
			gene_df = gene_df.merge(gene_qval_df,how='outer',on=peak2gene_cols[6:12])
		gene_pname_series = gene_groups_df.apply(lambda x: ";".join(
			str(s) for s in list(x["name"])))
		gene_pname_df = gene_pname_series.to_frame().reset_index()
		gene_pname_df.columns=peak2gene_cols[6:12]+['peak_name']
		gene_df = gene_df.merge(gene_pname_df,how='outer',on=peak2gene_cols[6:12])
		round1_gene_ann = ("%s/%s_r1_gene_annotations.txt" % (dir_name, prefix))
		gene_df.to_csv(round1_gene_ann, sep="\t", index=False)
		if verbose:
			sys.stdout.write("Number of Unique Genes that are Annotated: %i\n" % \
				gene_df.shape[0])
	return (round1_df)








