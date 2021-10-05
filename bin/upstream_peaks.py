#!/usr/bin/env python

import os
import sys
import pandas as pd
import numpy as np


def annotate(prefix, intergenic_peaks_df, gene_file, maximum_bp, outdir, \
             *positional_parameters, **keyword_parameters):
    """This function annotates upstream peaks to the closest TSS \
	within maximum number of bp. Note that if TSS is within range but on the same \
	strand as a closer TSS than we ignore the further TSS.
	mandatory parameters:
	* prefix - prefix for experiment
	* intergenic_peaks_df - data table that has all peaks that are not within genes
	* gene_file - bed file containing TSS sites.
	* maximum_bp - maximum bp upstream
	* outdir - output directory
	optional parameters:
	* ignore_conv_peaks - remove convergent peak information (peaks that annotate 
	to two different peaks) and put their info in a seperate file (default:False)
    """
    # default keyword parameters
    ignore_conv_peaks = False
    # user-set keyword parameters
    if ('ignore_conv_peaks' in keyword_parameters):
        ignore_conv_peaks = keyword_parameters['ignore_conv_peaks']
    # get tss_table from gene table
    tss_table = pd.read_csv(gene_file, sep='\t', header=None)
    tss_table.columns = ['tss_chr', 'tss_start', 'tss_stop', \
                         'gene_id', 'tss_score', 'tss_strand']
    tss_table.loc[tss_table["tss_strand"] == "+", "tss_stop"] = \
        tss_table.loc[tss_table["tss_strand"] == "+", "tss_start"] + 1
    tss_table.loc[tss_table["tss_strand"] == "-", "tss_start"] = \
        tss_table.loc[tss_table["tss_strand"] == "-", "tss_stop"] - 1
    # get table with joined info on tss and peaks
    peaks2genes_df = intergenic_peaks_df.merge(tss_table, how="inner", right_on="tss_chr", left_on="chr")
    peaks2genes_df["start"] = peaks2genes_df["start"].astype('int64')
    peaks2genes_df["stop"] = peaks2genes_df["stop"].astype('int64')
    peaks2genes_df["tss_start"] = peaks2genes_df["tss_start"].astype('int64')
    peaks2genes_df["tss_stop"] = peaks2genes_df["tss_stop"].astype('int64')
    # Only get rows that are upstream of genes
    # negative strand upstream
    neg_upstream_df = peaks2genes_df[((peaks2genes_df["tss_strand"] == "-") & \
                                      (peaks2genes_df["start"] >= peaks2genes_df["tss_stop"]))].copy()
    # positve strand upstream
    pos_upstream_df = peaks2genes_df[((peaks2genes_df["tss_strand"] == "+") & \
                                      (peaks2genes_df["stop"] <= peaks2genes_df["tss_start"]))].copy()
    # calculate distance for each group
    neg_upstream_df.loc[:, "distance_from_gene"] = ((neg_upstream_df.loc[:, "tss_stop"] - 1) - \
                                                    neg_upstream_df.loc[:, "start"])
    pos_upstream_df.loc[:, "distance_from_gene"] = ((pos_upstream_df.loc[:, "stop"] - 1) - \
                                                    pos_upstream_df.loc[:, "tss_start"])
    # Concatonate the groups together now that their distance is calculated
    upstream_df = neg_upstream_df.append(pos_upstream_df)
    upstream_df = upstream_df[upstream_df["distance_from_gene"] > -maximum_bp]
    # group by peak name
    # * if the peak has only one annotation then keep the annotation
    # * if the peak has 2+ annotations in same direction then
    # annotate the peak to the closest one
    # * if the peak has 2+ annotations in different directions than annotate
    # the peak to the 2 closest peaks if NOT ignoring convergent peaks. Otherwise,
    # do not annotate this peak
    upstream_filtered_arr = []
    convergent_peaks_arr = []
    upstream_pGroups = upstream_df.groupby(list(upstream_df.columns[0:6]))
    for p_info, p_group in upstream_pGroups:
        
        if p_group.shape[0] == 1:
            upstream_filtered_arr.append(p_group)
        else:
            if len(p_group["tss_strand"].unique()) == 1:
                closest_gene_series = p_group.loc[p_group["distance_from_gene"].idxmax(), :]
                closest_gene_df = pd.DataFrame(columns=upstream_df.columns)
                closest_gene_df = closest_gene_df.append(closest_gene_series)
                upstream_filtered_arr.append(closest_gene_df)
            else:
                # seperate the peaks upstream of positive strand genes
                # and upstream of negative strand genes
                neg_upstream_peak_df = p_group[p_group["tss_strand"] == "-"]
                pos_upstream_peak_df = p_group[p_group["tss_strand"] == "+"]
                # get the peaks closest to an upstream gene
                closest_gene_neg_series = \
                    neg_upstream_peak_df.loc[neg_upstream_peak_df["distance_from_gene"].idxmax(), :]
                closest_gene_pos_series = \
                    pos_upstream_peak_df.loc[pos_upstream_peak_df["distance_from_gene"].idxmax(), :]
                closest_gene_df = pd.DataFrame(columns=upstream_df.columns)
                closest_gene_df = closest_gene_df.append(closest_gene_neg_series)
                closest_gene_df = closest_gene_df.append(closest_gene_pos_series)
                convergent_peaks_arr.append(closest_gene_df)
                if not ignore_conv_peaks:
                    upstream_filtered_arr.append(closest_gene_df)
    upstream_filtered_df = pd.DataFrame(columns=upstream_df.columns)
    
    if len(upstream_filtered_arr) > 0:
        upstream_filtered_df = pd.concat(upstream_filtered_arr)
    upstream_filtered_df = upstream_filtered_df[
        list(intergenic_peaks_df.columns[0:6]) + list(tss_table.columns) + ["distance_from_gene"]]
    for upstream_filtered_cols in upstream_filtered_df.columns:
        upstream_filtered_df.loc[:,upstream_filtered_cols] = \
            upstream_filtered_df[upstream_filtered_cols].astype( \
            upstream_df[upstream_filtered_cols].dtype)
    #outfile = ("%s/%s_closest_id.txt" % (outdir, prefix))
    #upstream_filtered_df.to_csv(outfile, sep="\t", index=False, header=False)
    if len(convergent_peaks_arr) > 0:
        convergent_peaks_df = pd.concat(convergent_peaks_arr)
        convergent_peaks_file = ("%s/%s_convergent_upstream_peaks.txt" % (outdir,prefix))
        convergent_peaks_df = convergent_peaks_df[
            list(intergenic_peaks_df.columns[0:6]) + list(tss_table.columns) + ["distance_from_gene"]]
        convergent_peaks_df.to_csv(convergent_peaks_file, sep="\t", index=False, header=True)
    return (upstream_filtered_df)
