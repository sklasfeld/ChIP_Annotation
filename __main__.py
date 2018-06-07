#!/usr/bin/env python

import os
import argparse
import sys
import subprocess
import pandas as pd
import numpy as np
from collections import defaultdict
from collections import Counter
import bin.round1_annotation as round1_annotation
import bin.round2_annotation as round2_annotation
import bin.genome_locations as genome_locations


def cmd(cmd_str, speak):
    """print and run bash command"""
    if speak:
        sys.stdout.write("%s\n" % cmd_str)
        sys.stdout.flush()
    os.system(cmd_str)


def ifnum(in_str):
    """check if string contains numeric value"""
    try:
        num = float(in_str)
    except ValueError:
        return False
    return True


def text2vector(filename):
    """convert list delimited by newlines in a text file into a vector"""
    content = []
    with open(filename) as f:
        for line in f:
            content.append(line.strip())
    return (content)


def uniq_vals_incommon(list1, list2):
    """find unique values in common between two lists"""
    return list(set([x for x in list1 if x in list2]))


def wc(file_name, samtools_path=""):
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
    sys.stdout.write("%s\n" % cmd)
    sys.stdout.flush()
    ps = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT)
    output = ps.communicate()[0]
    break_out = output.split()
    sys.stdout.write("Count: %i\n" % int(break_out[0]))
    sys.stdout.flush()
    return int(break_out[0])


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="This is the main script for gene \
    annotation of ChIP peaks. This script first annotates ChIP peaks by \
    overlapping them with other ChIP peaks, DNase, MNase, or motifs. \
    Then it calls annotate the peaks to genes. The first round annotates peaks \
    by location relative to genes. Priority is given to intragenic peaks and \
    then peaks upstream of genes (limited by `filter_tss_upstream`) and then \
    peaks downstream of genes (limited by `filter_tts_downstream`). \
    The second round is dependent on genes that were found to be significantly \
    differentially expressed (DE). In round 2, peaks are annotated to upstream \
    or downstream DE genes (limited by `outlier_filter`). \
    Last, text files are output that are contain peak-centric and gene-centric \
    annotations")
    parser.add_argument('prefix', help='prefix for output file(eg. LFY)')
    parser.add_argument('dir_name', help='directory where output will go')
    parser.add_argument('bed_file', help='bed file containing ChIP peaks. Warning: \
    peaks should have original peak names.')
    parser.add_argument('gene_alist', help="A file that contains gene IDs \
    and their aliases. The file must be tab-delimited and have \
    atleast 2 columns with the labels `ID` and `Alias` \
    See `Araport11/Araport11_gene_info.txt file as an example \
    of the format. This is the default annotation")
    parser.add_argument('gene_bedfile', help='bed file with gene locations')
    parser.add_argument('-n', '--narrowpeak_file', help='narrowPeak file containing \
    ChIP peaks.')
    parser.add_argument('-bp', '--bedtools_path', help='path to bedtools', default="")
    parser.add_argument('-sp', '--samtools_path', help='path to samtools', default="")
    parser.add_argument('-tss', '--filter_tss_upstream', help='round1 annotation is \
    limited to upstream genes of this many bp (default: 1000)', default=1000, type=int)
    parser.add_argument('-tts', '--filter_tts_downstream', help='round1 annotation is \
    limited to downstream genes of this many bp (default: 100)', default=100, type=int)
    parser.add_argument('-icv', '--ignore_conv_peaks', help='do not output peak information \
    of peaks that annotate to two different peaks ', action='store_true')
    parser.add_argument('-of', '--outlier_filter', required=False, type=int,
                        help='maximum bp distance upstream/downstream of genes of \
                        outlier peaks (peaks not annotated in round 1) \
                        (default:10000)', default=10000)
    parser.add_argument('-ps', '--compareOldPeaks', nargs='*', help='list of peak \
    files that you want to compare with', required=False)
    parser.add_argument('-pn', '--compareOldPeaksNames', nargs='*', help='list of \
    prefixes for the peaks that you want to compare with. This must be equal to \
    "compareOldPeaks" variable (eg.LFY_seedlings)', required=False, default=[])
    parser.add_argument('-ops', '--oldPeakScores', help='if peak overlaps with \
    old peak then report score of old peak. Otherwise put NA', action='store_true')
    parser.add_argument('-rs', '--compareRNAseqs', nargs='*', help='list of text files \
    that contain a newline delimited list of genes that were "differentially expressed" \
     in a sample', required=False)
    parser.add_argument('-rn', '--compareRNAseqsNames', nargs='*', help='list of prefixes \
    for the RNA samples that you want to compare with. This must be equal to \
    "compareRNAseqs" variable (eg.RNA_seedlings)', required=False, default=[])
    parser.add_argument('-mf', '--motifFiles', nargs='*', help='list of bed files with \
    locations of TF motifs', required=False)
    parser.add_argument('-mn', '--motifNames', nargs='*', help='list giving names of \
    the type of motifs given in each bed file in --motifFiles. This list must be of \
    equal length to the "motifFiles" variable (eg.LFY1)', required=False, default=[])
    parser.add_argument('-df', '--dnase_files', nargs='*', help='list of bed files with \
    locations of DNase Hypersensitivity sites', required=False)
    parser.add_argument('-dn', '--dnase_names', nargs='*', help='list giving names for \
    each DNAse experiment corresponding to DNAse bed files in --dnase_files. This \
    list must be of equal length to the "dnase_files" variable (eg.DNAse_flowers)',
                        required=False, default=[])
    parser.add_argument('-nf', '--mnase_files', nargs='*', help='list of bed \
    files with locations of nucleosome binding sites (via MNase)', required=False)
    parser.add_argument('-nn', '--mnase_names', nargs='*', help='list giving names for \
    each MNase experiment corresponding to MNase bed files in --mnase_files. \
    This list must be of equal length to the "mnase_files" variable (eg.HighMock)',
                        required=False, default=[])
    parser.add_argument('-idr', '--globalIDRpval', help='global IDR pvalue used on \
    to get this output', default="")
    parser.add_argument('--keep_tmps', help='keep temp files', action='store_true')
    parser.add_argument('--verbose', help='echo processing', action='store_true')

    args = parser.parse_args()

    if args.compareOldPeaks or args.compareOldPeaksNames:
        if len(args.compareOldPeaks) != len(args.compareOldPeaksNames):
            sys.exit("compareOldPeaks and compareOldPeaksNames must be of equal length!")
    if args.compareRNAseqs or args.compareRNAseqsNames:
        if len(args.compareRNAseqs) != len(args.compareRNAseqsNames):
            sys.exit("compareRNAseqs and compareRNAseqsNames must be of equal length!")
    if args.motifFiles or args.motifNames:
        if len(args.motifFiles) != len(args.motifNames):
            sys.exit("motifFiles and motifNames must be of equal length!")
    if args.dnase_files or args.dnase_names:
        if len(args.dnase_files) != len(args.dnase_names):
            sys.exit("dnase_files and dnase_names must be of equal length!")
    if args.mnase_files or args.mnase_names:
        if len(args.mnase_files) != len(args.mnase_names):
            sys.exit("mnase_files and mnase_names must be of equal length!")

    dir_name = os.path.abspath(args.dir_name)
    gene_alist = args.gene_alist

    count_file = ("%s/%s_counts.txt" % (dir_name, args.prefix))

    counts_list = list()
    # roundOfAnnotation = {} # roundOfAnnotation[peak_name] = round_annotated
    # 1 = round 1 ()
    # 2 = round 2 (adopted peaks)
    # 3 = round 3 (orphan peaks)



    # ROUND 0: Characterize Peaks (independent of Gene Annotation)
    #

    ## print number of total peaks
    counts_list.append("number of total peaks\t%i" % wc(args.bed_file))
    peaks_df = pd.read_csv(args.bed_file, sep='\t', header=None, dtype=str)
    all_peak_names = list(set(peaks_df[3]))

    # print average length of peaks
    length_series = peaks_df.iloc[:, 2].astype('int64') - \
                    peaks_df.iloc[:, 1].astype('float64')
    counts_list.append("average length of peaks\t%f" % length_series.mean())

    # get peak dataframe
    peak_info_columns = ["chr", "start", "stop", "name", "signal", "strand", \
                         "fold_change", "pValue", "qValue", "summit"]
    bed_cols = peak_info_columns[:6]
    cur_peakfile = args.bed_file
    if args.narrowpeak_file:
        cur_peakfile = args.narrowpeak_file
        peaks_df = pd.read_csv(cur_peakfile, sep="\t", header=None, index_col=False, dtype=str)
        peaks_df.columns = peak_info_columns
    else:
        peaks_df.columns = bed_cols

    # compare experiments ChIP peaks with other ChIP experiments
    if args.compareOldPeaks:
        for oldpeakidx in range(0, len(args.compareOldPeaks)):
            old_peak_file = args.compareOldPeaks[oldpeakidx]  # bed file to compare with
            old_peak_prefix = args.compareOldPeaksNames[oldpeakidx]  # prefix for bed file to compare with
            compare_peak_file = ("%s/%s_peaks_found_in_%s.txt" % \
                                 (dir_name, args.prefix, old_peak_prefix))
            overlap_df_x = genome_locations.compare_bedfiles(args.bed_file, \
                                                             old_peak_file, compare_peak_file, verbal=args.verbose, \
                                                             bedtools_path=args.bedtools_path)
            overlap_df_x.columns = ["new_chr", "new_start", "new_stop", "new_name", \
                                    "new_score", "new_strand", "old_chr", "old_start", "old_stop", \
                                    "old_name", "old_score", "old_strand", "overlap"]
            # count number of peaks that overlap this other ChIP experiment
            overlap_count_x = len(overlap_df_x["new_name"].unique())
            counts_list.append("Number of peaks overlap %s\t%i" % \
                               (old_peak_prefix, overlap_count_x))
            newPeaksInOldPeaks = list(overlap_df_x["new_name"])
            peaks_df[old_peak_prefix] = overlap_df_x.apply(lambda row: row["new_name"] in newPeaksInOldPeaks, axis=1)
            peaks_df[old_peak_prefix].fillna(False, inplace=True)
            if not args.keep_tmps:
                os.remove(compare_peak_file)

        overlap_all_counts = len(peaks_df.loc[peaks_df.apply( \
            lambda row: all(row[args.compareOldPeaksNames]), axis=1), \
                                              "name"].unique())
        counts_list.append("Number of peaks overlap all old peaks\t%i" % \
                           (overlap_all_counts))

    # add motif info to peakwise files
    if args.motifFiles:
        for m_idx in range(0, len(args.motifFiles)):
            mFile = args.motifFiles[m_idx]
            mPrefix = args.motifNames[m_idx]
            motif_info_file = ("%s/peaks_in_motif_%s.txt" % (dir_name, \
                                                             mPrefix))
            motif_table = genome_locations.compare_bedfiles(args.bed_file, \
                                                            mFile, motif_info_file, verbal=args.verbose, \
                                                            bedtools_path=args.bedtools_path)
            ### START HERE!!!
            motif_table.columns = bed_cols + ['motif_chr', 'motif_start', \
                                              'motif_end', 'motif_name', 'motif_score', 'motif_strand', 'overlap']
            motif_table["location"] = motif_table["motif_chr"].map(str) + "_" \
                                      + motif_table["motif_start"].map(str) + "_" \
                                      + motif_table["motif_end"].map(str)
            motifsInPeak_df = motif_table.groupby(bed_cols)
            motifsLoc_series = motifsInPeak_df.apply( \
                lambda x: ";".join(str(s) for s in list(x["location"])))
            motifsLoc_df = motifsLoc_series.to_frame().reset_index()
            motifsLoc_df = motifsLoc_df.rename(columns={0: (mPrefix)})
            peaks_df = peaks_df.merge(motifsLoc_df, how="left", on=bed_cols)
            peakswithmotif_count = len(motifsLoc_df.loc[:, bed_cols].drop_duplicates())
            counts_list.append("number of peaks with motif %s \t%i" \
                               % (mPrefix, peakswithmotif_count))
            if not (args.keep_tmps):
                os.remove(motif_info_file)

    # add dnase info to peakwise files
    if args.dnase_files:
        for d_idx in range(0, len(args.dnase_files)):
            dFile = args.dnase_files[d_idx]
            dPrefix = args.dnase_names[d_idx]
            dnase_info_file = ("%s/peaks_in_dnase_%s.txt" % (dir_name, \
                                                             dPrefix))
            dnase_table = genome_locations.compare_bedfiles(args.bed_file, \
                                                            dFile, dnase_info_file, verbal=args.verbose, \
                                                            bedtools_path=args.bedtools_path)
            dnase_table_cols_needed = range(0, 9) + [dnase_table.shape[1] - 1]
            dnase_table = dnase_table.iloc[:, dnase_table_cols_needed]
            dnase_table.columns = bed_cols + ["dnase_chr", "dnase_start", \
                                              "dnase_stop", "overlap"]
            peaksInDHS = list(dnase_table["name"].unique())
            peaks_df.loc[:, dPrefix] = peaks_df.apply( \
                lambda row: row["name"] in peaksInDHS, axis=1)
            peaks_df[dPrefix].fillna(False, inplace=True)
            peakswithDHS_count = len(dnase_table.loc[:, bed_cols].drop_duplicates())
            counts_list.append("number of peaks within DNAse sites (%s) \t%i" \
                               % (dPrefix, peakswithDHS_count))
            if not (args.keep_tmps):
                os.remove(dnase_info_file)

    # add mnase info to peakwise files
    if args.mnase_files:
        for mn_idx in range(0, len(args.mnase_files)):
            mnFile = args.mnase_files[mn_idx]
            mnPrefix = args.mnase_names[mn_idx]
            mnase_info_file = ("%s/peaks_in_mnase_%s.txt" % (dir_name, \
                                                             mnPrefix))
            mnase_table = genome_locations.compare_bedfiles(args.bed_file, \
                                                            mnFile, mnase_info_file, verbal=args.verbose, \
                                                            bedtools_path=args.bedtools_path)
            mnase_table_cols_needed = range(0, 9) + [mnase_table.shape[1] - 1]
            mnase_table = mnase_table.iloc[:, mnase_table_cols_needed]
            mnase_table.columns = bed_cols + ["mnase_chr", "mnase_start", \
                                              "mnase_stop", "overlap"]
            peaksInNucleosomes = list(mnase_table["name"].unique())
            peaks_df.loc[:, mnPrefix] = peaks_df.apply( \
                lambda row: row["name"] in peaksInNucleosomes, axis=1)
            peaks_df[mnPrefix].fillna(False, inplace=True)
            peakswithNuc_count = len(mnase_table.loc[:, bed_cols].drop_duplicates())
            counts_list.append("number of peak with summits within MNAse sites (%s) \t%i" \
                               % (mnPrefix, peakswithNuc_count))
            if not (args.keep_tmps):
                os.remove(mnase_info_file)

    # ROUND 1

    ## annotate peaks that are :
    ##		intragenic
    ##		${noFilter_tss_upstream} bp upstream
    ##		${noFilter_TTS_downstream} bp downstream

    round1_peaks = round1_annotation.r1_annotate(gene_alist, args.gene_bedfile, \
                                                 args.bed_file, peaks_df, args.prefix, dir_name, \
                                                 bp_upstream_filter=args.filter_tss_upstream, \
                                                 bp_downstream_filter=args.filter_tts_downstream, \
                                                 ignore_conv_peak=args.ignore_conv_peaks, \
                                                 bedtools_path=args.bedtools_path, \
                                                 verbose=args.verbose)

    count_intragenic_genes = len( \
        round1_peaks.loc[round1_peaks["distance_from_gene"] == 0, \
                         bed_cols].drop_duplicates())
    count_upstream_genes_r1 = len( \
        round1_peaks.loc[round1_peaks["distance_from_gene"] < 0, \
                         bed_cols].drop_duplicates())
    count_downstream_genes_r1 = len( \
        round1_peaks.loc[round1_peaks["distance_from_gene"] > 0, \
                         bed_cols].drop_duplicates())
    if args.verbose:
        # number intragenic
        sys.stdout.write("Number of Peaks INSIDE Genes: %i\n" % \
                         count_intragenic_genes)
        # number upstream
        sys.stdout.write("Number of Peaks 1-%ibp UPSTREAM of Genes: %i\n" % \
                         (args.filter_tss_upstream, count_upstream_genes_r1))
        if args.filter_tts_downstream != 0:
            # number downstream
            sys.stdout.write("Number of Peaks 1-%ibp DOWNSTREAM of Genes: %i\n" % \
                             (args.filter_tts_downstream, count_downstream_genes_r1))

    ### print Number of Peaks INSIDE Genes
    round1_peaks["distance_from_gene"] = round1_peaks["distance_from_gene"].astype('float64')
    counts_list.append("peaks INSIDE genes in round 1:\t%i" % count_intragenic_genes)
    ### print Number of Peaks Upstream of Genes
    counts_list.append("peaks 1-%ibp UPSTREAM of genes in round 1:\t%i" % \
                       (args.filter_tss_upstream, count_upstream_genes_r1))
    ### print Number of Peaks Downstream of Genes
    if args.filter_tts_downstream > 0:
        counts_list.append("peaks 1-%ibp DOWNSTREAM of genes in round 1:\t%i" % \
                           (args.filter_tts_downstream, count_downstream_genes_r1))
    round1_peaks["roundOfAnnotation"] = 1

    ### Print number of peaks annotated in round 1
    round1_count = len(round1_peaks.loc[:, bed_cols].drop_duplicates())
    counts_list.append("peaks annotated in round1\t%i" % \
                       (round1_count))

    upstream_convergent_peaks_file = ("%s_convergent_upstream_peaks.txt" % args.prefix)
    if (not os.path.isfile(upstream_convergent_peaks_file)):
        counts_list.append("number of UPSTREAM convergent peaks in round 1:\t0")
    else:
        up_convergent_peaks_df = pd.read_csv(upstream_convergent_peaks_file, sep='\t', dtype=str)
        counts_list.append("number of UPSTREAM convergent peaks in round 1:\t%i" % \
                           len(list((up_convergent_peaks_df.loc[:, bed_cols].drop_duplicates()))))

    downstream_convergent_peaks_file = ("%s_convergent_downstream_peaks.txt" % args.prefix)
    if (not os.path.isfile(downstream_convergent_peaks_file)):
        counts_list.append("number of DOWNSTREAM convergent peaks in round 1:\t0")
    else:
        down_convergent_peaks_df = pd.read_csv(downstream_convergent_peaks_file, sep='\t', dtype=str)
        counts_list.append("number of DOWNSTREAM convergent peaks ignored in round 1:\t%i" % \
                           len(list((down_convergent_peaks_df.loc[:, bed_cols].drop_duplicates()))))

    # echo "ROUND 2"
    all_peaks_df = pd.DataFrame()

    # get peaks not found in round1_peaks
    outlier_df = peaks_df[~peaks_df["name"].isin(round1_peaks["name"])]

    # compare gene locations to outliers
    gene_bedfile_df = pd.read_csv(args.gene_bedfile, sep='\t', header=None)
    gene_bedfile_df.columns = ["gene_chr", "gene_start", "gene_stop", "gene_id", "gene_score", "gene_strand"]

    # annotate peaks that have not been annotated yet to DE genes within a bp limit
    round2_peaks = pd.DataFrame(columns=list(peaks_df.columns))
    orphan_peaks_df = pd.DataFrame(columns=list(peaks_df.columns))
    all_de_genes = []
    if args.compareRNAseqs:
        # RNAseqs_dict = defaultdict(list) # RNAseqs_dict[sampleName]=[genes]
        for rna_samp in range(0, len(args.compareRNAseqs)):
            rna_sample_file = args.compareRNAseqs[rna_samp]
            # rna_sample_name = args.compareRNAseqsNames[rna_samp]
            de_genes = text2vector(rna_sample_file)
            # RNAseqs_dict[rna_sample_name] = de_genes
            all_de_genes += de_genes
        all_de_genes = list(set(all_de_genes))
        de_genes_df = gene_bedfile_df[gene_bedfile_df['gene_id'].isin(all_de_genes)]

        round2_ann_file = ("%s/%s_r2_peak_annotations.txt" % (dir_name, args.prefix))
        round2_peaks = round2_annotation.r2_annotate(gene_alist, de_genes_df, outlier_df,
                                                     args.outlier_filter, round2_ann_file)
        round2_peaks["roundOfAnnotation"] = 2
        orphan_peaks_df = outlier_df[~outlier_df["name"].isin(round2_peaks["name"])]

        ### Print number of peaks annotated in round 2
        round2_count = len(round2_peaks.loc[:, bed_cols].drop_duplicates())
        counts_list.append("peaks annotated in round2\t%i" % \
                           (round2_count))
        ### Print number of peaks annotated in round 1 AND 2
        counts_list.append("total peaks annotated\t%i" % (round1_count + round2_count))
        # get dataframe with all peak annotations including unassigned peaks
        all_peaks_frame = [round1_peaks, round2_peaks, orphan_peaks_df]
        all_peaks_df = pd.concat(all_peaks_frame)

    else:
        orphan_peaks_df = outlier_df
        # get dataframe with all peak annotations including unassigned peaks
        all_peaks_frame = [round1_peaks, orphan_peaks_df]
        all_peaks_df = pd.concat(all_peaks_frame)

    # reorganize columns
    all_peaks_col_order = bed_cols + ['roundOfAnnotation']
    if args.narrowpeak_file:
        all_peaks_col_order += ['qValue', 'summit']
    all_peaks_col_order = all_peaks_col_order + \
                          ['gene_id', 'Alias'] + \
                          ['distance_from_gene'] + \
                          args.compareOldPeaksNames + \
                          args.motifNames + \
                          args.dnase_names + \
                          args.mnase_names

    all_peaks_df = all_peaks_df.loc[:, all_peaks_col_order]
    all_peaks_df = all_peaks_df.rename(columns={'Alias': 'gene_name'})

    # ROUND 3: print out peaks that remain unassigned (orphan peaks)
    orphan_bed_f = ("%s/%s_unassigned.bed" % (dir_name, args.prefix))
    orphan_peaks_df.iloc[:, 0:6].to_csv(orphan_bed_f, sep="\t", header=False, index=False)
    if args.narrowpeak_file:
        orphan_np_f = ("%s/%s_unassigned.narrowPeak" % (dir_name, args.prefix))
        orphan_peaks_df.iloc[:, 0:10].to_csv(orphan_np_f, sep="\t", header=False, index=False)

    # Assign DE genes to peaks

    ## How many genes were annotated in each round (indepedent of DE)
    gene_count_r1 = len(all_peaks_df.loc[all_peaks_df["roundOfAnnotation"] == 1, "gene_id"].unique())
    counts_list.append("genes annotated in round1\t%i" % \
                       (gene_count_r1))
    if args.compareRNAseqs:
        gene_count_r2 = len(all_peaks_df.loc[all_peaks_df["roundOfAnnotation"] == 2, "gene_id"].unique())
        counts_list.append("genes annotated in round2\t%i" % \
                           (gene_count_r2))
        counts_list.append("total genes annotated\t%i" % \
                           (gene_count_r1 + gene_count_r2))

        ## Which genes are DE according to the RNA analysis?
        for rna_samp in range(0, len(args.compareRNAseqs)):
            rna_sample_file = args.compareRNAseqs[rna_samp]
            rna_sample_name = args.compareRNAseqsNames[rna_samp]
            de_genes = text2vector(rna_sample_file)
            all_peaks_df[rna_sample_name] = all_peaks_df.apply(lambda row: row["gene_id"] in de_genes, axis=1)
            all_peaks_df[rna_sample_name].fillna(False, inplace=True)
            peaks_ann_to_de_counts = len(all_peaks_df.loc[all_peaks_df[rna_sample_name] \
                                                          == True, "gene_id"].unique())
            round1_splice = all_peaks_df[all_peaks_df['roundOfAnnotation'] == 1]
            r1_peaks_ann_to_de_counts = len(round1_splice.loc[round1_splice[rna_sample_name] \
                                                              == True, "gene_id"].unique())
            counts_list.append("total differentially expressed genes in %s\t%i" % \
                               (rna_sample_name, len(de_genes)))
            counts_list.append("# genes annotated to %s in round1\t%i" % \
                               (rna_sample_name, r1_peaks_ann_to_de_counts))
            counts_list.append("# genes annotated to %s in both rounds\t%i" % \
                               (rna_sample_name, peaks_ann_to_de_counts))
        ## see how many gene are found in all rna samples
        counts_list.append("differentially expressed genes given all RNA samples\t%i" % \
                           len(all_de_genes))
        de_round1 = len([x for x in \
                         all_peaks_df.loc[( \
                                              (all_peaks_df[args.compareRNAseqsNames].any(axis=1)) & \
                                              (all_peaks_df[u'roundOfAnnotation'] == 1)), \
                                          "gene_id"].unique() if str(x) != 'nan'])
        counts_list.append("# DE genes annotated given all RNA samples in round1\t%i" % \
                           de_round1)
        de_genes_total = len([x for x in \
                              all_peaks_df.loc[all_peaks_df[args.compareRNAseqsNames].all(axis=1), \
                                               "gene_id"].unique() if str(x) != 'nan'])
        counts_list.append("# DE genes annotated given all RNA samples in both rounds\t%i" % \
                           de_genes_total)

    # Print out Peak-centric datatable...
    peak_group_cols = bed_cols + ['roundOfAnnotation']
    if args.narrowpeak_file:
        peak_group_cols += ['qValue', 'summit']
    peak_group_cols = peak_group_cols + \
                      args.compareOldPeaksNames + \
                      args.motifNames + \
                      args.dnase_names + \
                      args.mnase_names
    ## note: non-peak_centric columns = gene_id, Alias, distance_from_gene
    peak_grouped_df = all_peaks_df.groupby(peak_group_cols)
    peak_centric_cols_df = all_peaks_df.loc[:, peak_group_cols]
    ## list number of genes that annotate to peak
    peak_nGenes_series = peak_grouped_df.apply(lambda x: len(x["gene_id"].unique()))
    peak_nGenes_df = peak_nGenes_series.to_frame().reset_index()
    peak_nGenes_df.columns = peak_group_cols + ['numGenes']
    peak_ann_df = peak_centric_cols_df.merge(peak_nGenes_df, how="left", \
                                             on=peak_group_cols)
    ## list gene ids that annotate to peaks
    peak_gid_series = peak_grouped_df.apply(lambda x: ";".join(str(s) for s in list(x["gene_id"])))
    peak_gid_df = peak_gid_series.to_frame().reset_index()
    peak_gid_df.columns = peak_group_cols + ['gene_id']
    peak_ann_df = peak_ann_df.merge(peak_gid_df, how='left', on=peak_group_cols)
    ## list gene names that annotate to peaks
    peak_gname_series = peak_grouped_df.apply(lambda x: ";".join(str(s) for s in list(x["gene_name"])))
    peak_gname_df = peak_gname_series.to_frame().reset_index()
    peak_gname_df.columns = peak_group_cols + ['gene_name']
    peak_ann_df = peak_ann_df.merge(peak_gname_df, how='outer', on=peak_group_cols)
    ## list distances of genes that annotate to peaks
    peak_dist_series = peak_grouped_df.apply(lambda x: ";".join(str(s) for s in list(x["distance_from_gene"])))
    peak_dist_df = peak_dist_series.to_frame().reset_index()
    peak_dist_df.columns = peak_group_cols + ['distance_from_gene']
    peak_ann_df = peak_ann_df.merge(peak_dist_df, how='outer', on=peak_group_cols)
    ## list True if any of the genes that the peak is annotated to is DE
    ## in each sample
    if args.compareRNAseqs:
        for rna_name in args.compareRNAseqsNames:
            peak_degene_series = peak_grouped_df.apply(lambda x: x[rna_name].any())
            peak_degene_df = peak_degene_series.to_frame().reset_index()
            peak_degene_df.columns = peak_group_cols + [rna_name]
            peak_ann_df = peak_ann_df.merge(peak_degene_df, how='left', \
                                            on=peak_group_cols)
            peak_ann_df[rna_name].fillna(False, inplace=True)
    ## reorganize data-table to show gene columns closer to peak columns and
    ## extra info towards the later columns
    peak2gene_info_cols = ['numGenes', 'gene_id', 'gene_name','distance_from_gene']
    peak_col_order = peak_group_cols[0:7] + peak2gene_info_cols + \
                     args.compareRNAseqsNames + peak_group_cols[7:]
    peak_ann_df = peak_ann_df.loc[:, peak_col_order]
    pd.set_option('float_format', '{:.0f}'.format)
    peak_out = ("%s/%s_peakwise_ann.txt" % (dir_name, args.prefix))
    peak_ann_df = peak_ann_df.drop_duplicates()
    peak_ann_df.to_csv(peak_out, sep="\t", index=False, na_rep="NA")

    # Print out Gene-centric datatable...
    gene_group_cols = ["gene_id", "gene_name"] + \
                      args.compareRNAseqsNames
    gene_groups_df = all_peaks_df.groupby(gene_group_cols)
    ## list number of genes that annotate to peak
    gene_nPeaks_series = gene_groups_df.apply(lambda x: x.shape[0])
    gene_nPeaks_df = gene_nPeaks_series.to_frame().reset_index()
    gene_nPeaks_df.columns = gene_group_cols + ['numPeaks']
    ## list peak that annotate to gene
    all_peaks_df["peaks"] = all_peaks_df["name"].map(str) + "=" + \
                            all_peaks_df["chr"].map(str) + ":" \
                            + all_peaks_df["start"].map(str) + "-" \
                            + all_peaks_df["stop"].map(str)
    gene_pname_series = gene_groups_df.apply(lambda x: ";".join(str(s) for s in list(x["peaks"])))
    gene_pname_df = gene_pname_series.to_frame().reset_index()
    gene_pname_df.columns = gene_group_cols + ['peaks']
    gene_ann_df = gene_nPeaks_df.merge(gene_pname_df, how='outer', on=gene_group_cols)
    # round of annotation
    gene_roundAnn_series = gene_groups_df.apply(lambda x: ";".join(str(s) for s in list(x['roundOfAnnotation'])))
    gene_roundAnn_df = gene_roundAnn_series.to_frame().reset_index()
    gene_roundAnn_df.columns = gene_group_cols + ['roundOfAnnotation']
    gene_ann_df = gene_ann_df.merge(gene_roundAnn_df, how='outer', on=gene_group_cols)
    # distance
    gene_dist_series = gene_groups_df.apply(lambda x: ";".join(str(s) for s in list(x['distance_from_gene'])))
    gene_dist_df = gene_dist_series.to_frame().reset_index()
    gene_dist_df.columns = gene_group_cols + ['distance_to_peak']
    gene_ann_df = gene_ann_df.merge(gene_dist_df, how='outer', on=gene_group_cols)
    # list qvalues of peaks that annotate to the gene
    if args.narrowpeak_file:
        gene_qval_series = gene_groups_df.apply(lambda x: ";".join(str(s) for s in list(x['qValue'])))
        gene_qval_df = gene_qval_series.to_frame().reset_index()
        gene_qval_df.columns = gene_group_cols + ['qValue']
        gene_ann_df = gene_ann_df.merge(gene_qval_df, how='outer', on=gene_group_cols)
    ## list True if any of the peaks that the gene is annotated to is found in a
    ## a different ChIP sample (compareOldPeaks)
    if args.compareOldPeaks:
        for peak_sample in args.compareOldPeaksNames:
            gene_peakcomp_series = gene_groups_df.apply(lambda x: x[peak_sample].any())
            gene_peakcomp_df = gene_peakcomp_series.to_frame().reset_index()
            gene_peakcomp_df.columns = gene_group_cols + [peak_sample]
            gene_ann_df = gene_ann_df.merge(gene_peakcomp_df, how='left', \
                                            on=gene_group_cols)
            gene_ann_df[peak_sample].fillna(False, inplace=True)
    ## list True if any of the peaks that the gene is annotated to is found in a
    ## a motif (motifFiles)
    if args.motifFiles:
        for motif_col in args.motifNames:
            gene_motif_series = gene_groups_df.apply(lambda x: x[motif_col].notnull().any())
            gene_motif_df = gene_motif_series.to_frame().reset_index()
            gene_motif_df.columns = gene_group_cols + [motif_col]
            gene_ann_df = gene_ann_df.merge(gene_motif_df, how='left', \
                                            on=gene_group_cols)
            gene_ann_df[motif_col].fillna("NA", inplace=True)
    ## list True if any of the peaks that the gene is annotated to is found in a
    ## a DHS (dnase_files)
    if args.dnase_files:
        for dhs_col in args.dnase_names:
            gene_dhs_series = gene_groups_df.apply(lambda x: x[dhs_col].any())
            gene_dhs_df = gene_dhs_series.to_frame().reset_index()
            gene_dhs_df.columns = gene_group_cols + [dhs_col]
            gene_ann_df = gene_ann_df.merge(gene_dhs_df, how='left', \
                                            on=gene_group_cols)
            gene_ann_df[dhs_col].fillna(False, inplace=True)
    ## list True if any of the peaks that the gene is annotated to is found in a
    ## a nucleosome (mnase_files)
    if args.mnase_files:
        for mnase_samp in args.mnase_names:
            gene_mnase_series = gene_groups_df.apply(lambda x: x[mnase_samp].any())
            gene_mnase_df = gene_mnase_series.to_frame().reset_index()
            gene_mnase_df.columns = gene_group_cols + [mnase_samp]
            gene_ann_df = gene_ann_df.merge(gene_mnase_df, how='left', \
                                            on=gene_group_cols)
            gene_ann_df[mnase_samp].fillna(False, inplace=True)

    gene_out = ("%s/%s_genewise_ann.txt" % (dir_name, args.prefix))
    gene_ann_df = gene_ann_df.drop_duplicates()
    gene_ann_df.to_csv(gene_out, sep="\t", index=False, na_rep="NA")

    olog_file = open(count_file, 'w')
    for i in range(0, len(counts_list)):
        olog_file.write("%s\n" % (counts_list[i]))

    olog_file.close()
