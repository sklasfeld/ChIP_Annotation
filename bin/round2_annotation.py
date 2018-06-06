#!/usr/bin/env python

import pandas as pd

def r2_annotate(gene_alist, gene_df, peak_df, maxdist, out):
	""" calculate distance between features

	mandatory parameters:
	- gene_df: gene dataframe with strand info; distance based on this file
	- peak_df: peaks dataframe to compare with
	- maxdist: maximum distance to report, if no maximum, then input -1
	- out: tab-delimited file with round2 annotations
	"""
	narrowPeak_boolean = False
	if 'qValue' in peak_df.columns:
		narrowPeak_boolean = True
	peaks2genes_df = peak_df.merge(gene_df, left_on="chr", right_on="gene_chr")
	
	# negative strand upstream
	neg_upstream_df = peaks2genes_df[((peaks2genes_df["gene_strand"]=="-") & \
		(peaks2genes_df["start"] >= peaks2genes_df["gene_stop"]))]
	# negative strand downstream 
	neg_downstream_df = peaks2genes_df[((peaks2genes_df["gene_strand"]=="-") & \
		(peaks2genes_df["stop"] <= peaks2genes_df["gene_start"]))]
	# positive strand downstream
	pos_downstream_df = peaks2genes_df[((peaks2genes_df["gene_strand"]=="+") & \
		(peaks2genes_df["start"] >= peaks2genes_df["gene_stop"]))]
	# positve strand upstream
	pos_upstream_df = peaks2genes_df[((peaks2genes_df["gene_strand"]=="+") & \
		(peaks2genes_df["stop"] <= peaks2genes_df["gene_start"]))]
	# intragenic
	intragenic_df = peaks2genes_df[(( \
		(peaks2genes_df["start"] >= peaks2genes_df["gene_start"]) & \
		(peaks2genes_df["start"] < peaks2genes_df["gene_stop"])) | \
		((peaks2genes_df["stop"] > peaks2genes_df["gene_start"]) & \
		(peaks2genes_df["stop"] <= peaks2genes_df["gene_stop"])))]


	# calculate distance for each group
	neg_upstream_df.loc[:,"distance_from_gene"] = ((neg_upstream_df["gene_stop"] - 1) - \
		neg_upstream_df["start"])
	neg_downstream_df.loc[:,"distance_from_gene"] = (neg_downstream_df["gene_start"] - \
		(neg_downstream_df["stop"]-1))
	pos_downstream_df.loc[:,"distance_from_gene"] = (pos_downstream_df["start"] - \
		(pos_downstream_df["gene_stop"]-1))
	pos_upstream_df.loc[:,"distance_from_gene"] = ((pos_upstream_df["stop"] - 1) - \
		pos_upstream_df["gene_start"])
	intragenic_df["distance_from_gene"] = 0

	# filter by maxdist
	neg_upstream_df = neg_upstream_df[-(maxdist) <= neg_upstream_df["distance_from_gene"]]
	neg_downstream_df = neg_downstream_df[maxdist >= neg_downstream_df["distance_from_gene"]]
	pos_upstream_df = pos_upstream_df[-(maxdist) <= pos_upstream_df["distance_from_gene"]]
	pos_downstream_df = pos_downstream_df[maxdist >= pos_downstream_df["distance_from_gene"]]
	
	round2_frames = [neg_upstream_df, neg_downstream_df, pos_upstream_df, \
		pos_downstream_df, intragenic_df]

	round2_df = pd.concat(round2_frames)

	# add alias info about each gene
	geneAnn_df = pd.read_csv(gene_alist, sep='\t', dtype=str)
	round2_alias_df = round2_df.merge(geneAnn_df, how='left', \
		left_on='gene_id', right_on='ID')
	# sort the dataframe by peak location
	round2_alias_df["start"] = round2_alias_df["start"].astype('int64')
	round2_alias_df.sort_values(by=["chr","start"], axis=0, ascending=True, inplace=True)
	

	peaks_group_cols = list(peak_df.columns[0:6])
	if narrowPeak_boolean:
		peaks_group_cols = peaks_group_cols+["qValue"]+list(peaks_df.columns[10:])
	else:
		peaks_group_cols = peaks_group_cols+list(peaks_df.columns[6:])
	peak_groups_df = round2_alias_df.groupby(peaks_group_cols)
	peak_nGenes_series = peak_groups_df.apply(lambda x: len(x["gene_id"].unique()))
	peak_nGenes_df = peak_nGenes_series.to_frame().reset_index()
	peak_nGenes_df.columns=peaks_group_cols+['numGenes']
	peak_gid_series = peak_groups_df.apply(lambda x: ";".join(str(s) for s in list(x["gene_id"])))
	peak_gid_df = peak_gid_series.to_frame().reset_index()
	peak_gid_df.columns=peaks_group_cols+['gene_id']
	peak_gname_series = peak_groups_df.apply(lambda x: ";".join(str(s) for s in list(x["Alias"])))
	peak_gname_df = peak_gname_series.to_frame().reset_index()
	peak_gname_df.columns=peaks_group_cols+['gene_name']
	peak_ann_df = peak_nGenes_df.merge(peak_gid_df,how='outer',on=peaks_group_cols)
	peak_ann_df = peak_ann_df.merge(peak_gname_df,how='outer',on=peaks_group_cols)
	peak2gene_info_cols = ['numGenes','gene_id','gene_name']
	column_order = peak_ann_df.columns
	if narrowPeak_boolean:
		column_order= list(peak_df.columns[0:6]) + ["qValue"] + \
			peak2gene_info_cols+list(peaks_df.columns[10:])
		
	else:
		column_order= list(peak_df.columns[0:6]) + \
			peak2gene_info_cols+list(peaks_df.columns[6:])
	peak_ann_df = peak_ann_df.loc[:,column_order]

	pd.set_option('float_format', '{:.0f}'.format)
	peak_ann_df.to_csv(out, sep="\t", index=False, na_rep="NA")

	return(round2_alias_df)
