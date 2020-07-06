#!/usr/bin/env python

import pandas as pd
import sys

def r2_annotate(gene_df, peaks_df, maxdist, out, \
	*positional_parameters, **keyword_parameters):
	""" calculate distance between features

	mandatory parameters:
	- gene_df: gene dataframe with strand info; distance based on this file
	- peaks_df: peaks dataframe to compare with
	- maxdist: maximum distance to report, if no maximum, then input -1
	- out: tab-delimited file with round2 annotations
	"""

	narrowPeak_boolean = False
	if 'qValue' in peaks_df.columns:
		narrowPeak_boolean = True
	# check the chromosome names match
	peaks_chrs_list = list(peaks_df["chr"].unique())
	genes_chrs_list = list(gene_df["gene_chr"].unique())
	inFile1Not2 = [x for x in peaks_chrs_list if x not in genes_chrs_list]
	inFile1and2 = [x for x in peaks_chrs_list if x in genes_chrs_list]
	if len(inFile1and2) == 0 or len(inFile1Not2) > 0:
		chrInFile1_str = ",".join(peaks_chrs_list)
		chrInFile2_str = ",".join(genes_chrs_list)
		if len(inFile1and2) == 0:
			err_msg = (("WARNING: NONE of the chromosomes columns in " + 
				"the narrowPeak file match the chromosomes " + 
				"in the gene file.\nThe chromosomes in the bed " + 
				"file are: %s.\nThe chromosome in gene file are: %s.") % \
			( chrInFile1_str, chrInFile2_str))
			sys.stderr.write(err_msg)
			return(pd.DataFrame(columns=
				list(peaks_df.columns) + list(gene_df.columns)))
		else:
			err_msg = (("WARNING: The chromosomes columns in the narrowPeak file " + \
				" do not match the chromosomes in the gene file. The " + \
				"chromosomes in the bed file are: %s. " + \
				"The chromosome in gene file are: %s.") % \
			( chrInFile1_str, chrInFile2_str))
			sys.stderr.write(err_msg)


	peaks2genes_df = peaks_df.merge(gene_df, left_on="chr", right_on="gene_chr")
	

	# negative strand upstream
	neg_upstream_df = peaks2genes_df[((peaks2genes_df["gene_strand"]=="-") & \
		(peaks2genes_df["start"] >= peaks2genes_df["gene_stop"]))].copy()
	# negative strand downstream 
	neg_downstream_df = peaks2genes_df[((peaks2genes_df["gene_strand"]=="-") & \
		(peaks2genes_df["stop"] <= peaks2genes_df["gene_start"]))].copy()
	# positive strand downstream
	pos_downstream_df = peaks2genes_df[((peaks2genes_df["gene_strand"]=="+") & \
		(peaks2genes_df["start"] >= peaks2genes_df["gene_stop"]))].copy()
	# positve strand upstream
	pos_upstream_df = peaks2genes_df[((peaks2genes_df["gene_strand"]=="+") & \
		(peaks2genes_df["stop"] <= peaks2genes_df["gene_start"]))].copy()
	# intragenic
	intragenic_df = peaks2genes_df[(( \
		(peaks2genes_df["start"] >= peaks2genes_df["gene_start"]) & \
		(peaks2genes_df["start"] < peaks2genes_df["gene_stop"])) | \
		((peaks2genes_df["stop"] > peaks2genes_df["gene_start"]) & \
		(peaks2genes_df["stop"] <= peaks2genes_df["gene_stop"])))].copy()

	# calculate distance for each group
	#neg_upstream_df.loc[:,"start"].is_copy = False
	neg_upstream_df.loc[:,"distance_from_gene"] = ((neg_upstream_df.loc[:,"gene_stop"] - 1) - 
		(neg_upstream_df.loc[:,"start"]))
	neg_downstream_df.loc[:,"distance_from_gene"] = (neg_downstream_df.loc[:,"gene_start"] - 
		(neg_downstream_df.loc[:,"stop"]-1))
	pos_downstream_df.loc[:,"distance_from_gene"] = (pos_downstream_df.loc[:,"start"] - 
		(pos_downstream_df.loc[:,"gene_stop"]-1))
	pos_upstream_df.loc[:,"distance_from_gene"] = ((pos_upstream_df.loc[:,"stop"] - 1) - 
		pos_upstream_df.loc[:,"gene_start"])
	intragenic_df["distance_from_gene"] = 0

	# filter by maxdist
	neg_upstream_df = neg_upstream_df[-(maxdist) <= neg_upstream_df["distance_from_gene"]]
	neg_downstream_df = neg_downstream_df[maxdist >= neg_downstream_df["distance_from_gene"]]
	pos_upstream_df = pos_upstream_df[-(maxdist) <= pos_upstream_df["distance_from_gene"]]
	pos_downstream_df = pos_downstream_df[maxdist >= pos_downstream_df["distance_from_gene"]]
	
	round2_frames = [neg_upstream_df, neg_downstream_df, pos_upstream_df, \
		pos_downstream_df, intragenic_df]

	#round2_df = pd.concat(round2_frames, sort=False)
	round2_df = pd.concat(round2_frames)

	if len(round2_df) > 0:
		# sort the dataframe by peak location
		round2_df.loc[:,"start"] = round2_df["start"].astype('int64')
		round2_df.sort_values(by=["chr","start"], axis=0, ascending=True, inplace=True)
		# print information for each peak
		peaks_group_cols = list(peaks_df.columns[0:6])
		## group by peak info

		peak_groups_df = round2_df.groupby(peaks_group_cols)
		## peak centric columns
		peaks_centric_cols = peaks_group_cols
		if narrowPeak_boolean:
			peaks_centric_cols = peaks_centric_cols +["qValue"]+list(peaks_df.columns[10:])
		else:
			peaks_centric_cols = peaks_centric_cols + list(peaks_df.columns[6:])
		peak_ann_df = round2_df.loc[:,peaks_centric_cols]
		## put columns that are not peak centric into a peak context	
		peak_nGenes_series = peak_groups_df.apply(lambda x: len(x["gene_id"].unique()))
		peak_nGenes_df = peak_nGenes_series.to_frame().reset_index()
		peak_nGenes_df.columns=peaks_group_cols+['numGenes']
		peak_ann_df = peak_ann_df.merge(peak_nGenes_df,how='outer',on=peaks_group_cols)
		peak_gid_series = peak_groups_df.apply(lambda x: ";".join(str(s) for s in list(x["gene_id"])))
		peak_gid_df = peak_gid_series.to_frame().reset_index()
		peak_gid_df.columns=peaks_group_cols+['gene_id']
		peak_ann_df = peak_ann_df.merge(peak_gid_df,how='outer',on=peaks_group_cols)
		peak2gene_info_cols = ['numGenes','gene_id']
		column_order = peak_ann_df.columns
		if narrowPeak_boolean:
			column_order= list(peaks_df.columns[0:6]) + ["qValue"] + \
				peak2gene_info_cols+list(peaks_df.columns[10:])
			
		else:
			column_order= list(peaks_df.columns[0:6]) + \
				peak2gene_info_cols+list(peaks_df.columns[6:])
		peak_ann_df = peak_ann_df.loc[:,column_order]

		pd.set_option('float_format', '{:.2f}'.format)
		peak_ann_df.to_csv(out, sep="\t", index=False, na_rep="NA")

	return(round2_df)
