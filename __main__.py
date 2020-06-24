#!/usr/bin/env python

import os
import argparse
import sys
import pandas as pd
import numpy as np
from collections import defaultdict
from collections import Counter
import bin.round1_annotation as round1_annotation
import bin.round2_annotation as round2_annotation
import bin.genome_locations as genome_locations
pd.options.mode.chained_assignment = 'raise'


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

def listOrFalse(res_series):
    false_list = pd.Series([True if x == False or x == "False" else False for x in res_series])
    if false_list.all():
        return False
    else:
        string_list = ";".join([str(y) for y in res_series if (y != False and y != "False")])
        return string_list

def restricted_float(x):
    try:
        x = float(x)
    except ValueError:
        raise argparse.ArgumentTypeError("%r not a floating-point literal" % (x,))

    if x < 0.0 or x > 1.0:
        raise argparse.ArgumentTypeError("%r not in range [0.0, 1.0]"%(x,))
    return x

def gene_overlap_dnase(row):
    if row.start_b <= row.gene_start:
        if row.stop_b <= row.gene_start:
            return(0)
        elif row.stop_b >= row.gene_stop:
            return (row.gene_stop-row.gene_start)
        else:
            return (row.stop_b-row.gene_start)
    elif row.start_b < row.gene_stop:
        if row.stop_b >= row.gene_stop:
            return (row.gene_stop-row.start_b)
        else:
            return (row.stop_b-row.start_b)
    else:
        return (0)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="At its core the purpose of \
        this code is to annotate ChIP peak regions to genes. \
        Other parameters allow for users to annotate both the peaks and the \
        genes further with external datasets.")
    parser.add_argument('prefix', help='prefix for output file(eg. LFY)')
    parser.add_argument('dir_name', help='directory where output will go')
    parser.add_argument('bed_file', help='bed file containing ChIP peaks. \
        Warning: peaks should have original peak names.')
    parser.add_argument('gene_bedfile', help='bed file with gene locations')
    parser.add_argument('-n', '--narrowpeak_file', help='narrowPeak file \
        containing ChIP peaks.')
    parser.add_argument('-pi', '--percent_inter', help='round1 \
        annotation only annotates genes that overlap at least this \
        percent of the gene. For example if this is set to 1 then \
        the peak must overlap the gene 100%% to annotate to it. \
        If it is set to .5 then the peak must overlap \
        at least half. (default:0)', default=0, type=restricted_float)
    parser.add_argument('-tss', '--filter_tss_upstream', help='round1 \
        annotation is limited to upstream genes of this many bp \
        (default: 3000)', default=3000, type=int)
    parser.add_argument('-tts', '--filter_tts_downstream', help='round1 \
        annotation is limited to downstream genes of this many bp \
        (default: 100)', default=100, type=int)
    parser.add_argument('-icv', '--ignore_conv_peaks', help='do not output \
        peak information of peaks that annotate to two different peaks ',
        action='store_true')
    parser.add_argument('-r2', '--round2ann', help='Annotate outlier peaks \
        (peaks not annotated in round 1) to differentially expressed genes. \
        For this parameter, set the RNA differential expression output files',
        nargs='*')
    parser.add_argument('-rl', '--rnaDESignificance', nargs=3, \
        help='using the tab-delimited text files set by the `round2ann` parameter, \
        we can limit which genes are defined as differentially expressed. \
        This parameter requires 3 values:[column] [min/max] [value]. For example, \
        if you want to limit DE genes to genes output from DESeq2 with adjp<=0.01 \
        then set this to `--rnaDESignificance adjp max 0.01`. Another example is \
        if you want to limit DE genes to genes output from DESeq2 with \
        log2FoldChange>=1 then set this to `--rnaDESignificance log2FoldChange min 1`.', \
        required=False)
    parser.add_argument('-of', '--outlier_filter', required=False, type=int,
        help='If `--round2ann` is set, this is the \
        maximum bp distance upstream/downstream of genes of outlier peaks \
        (peaks not annotated \
        in round 1) (default:10000).', default=10000)
    parser.add_argument('-pb', '--comparePeaksToBeds', nargs='*', help='list of \
        bed files that you want to compare with the peaks in this sample.', \
        required=False)
    parser.add_argument('-pbn', '--comparePeaksToBedsNames', nargs='*', \
        help='list of prefixes for the bed files that you want to compare with. \
        This must be equal to "comparePeaksToBeds" variable \
        (eg.LFY_seedlings)', required=False, default=[])
    parser.add_argument('-pbf', '--compBedScores', help='if peak overlaps with \
        a row in one of the bed files in `comparePeaksToBeds` then report a \
        feature of the old peak. Otherwise put NA. Any \
        column in a narrowPeak file can be reported (see \
        https://genome.ucsc.edu/FAQ/FAQformat.html#format12 for description \
        of each column).', choices=["chrom", "chromStart", "chromEnd", "name", \
        "score", "strand", \
        "singalValue", "pValue", "qValue", "peak"])
    parser.add_argument('-sb', '--compareSumRegToBed', nargs='*', 
        help='list of bed files that you want to compare with the region \
        around the ChIP summit site. Note that for this function to work, \
        the following other parameters must also be set correctly: \
        `narrowpeak_file`, `compareSumRegToBedColNames`, and \
        `summitRegion`', required=False)
    parser.add_argument('-sbn', '--compareSumRegToBedColNames', nargs='*', 
        help='list of column names for the output of the comparison \
        between the features in the bed file(s) in `compareSumRegToBed` \
        and the ChIP summit regions (eg.LFY1_motif)', required=False, default=[])
    parser.add_argument('-sbr', '--summitRegion', nargs='*', type=np.int64, 
        help='To run `compareSumRegToBed` the user must set the \
        boundaries around the peak summit to compare to. The specific \
        format for this parameter is: \
        `[compareSumRegToBedColNames],[bp_upstream],[bp_downstream]`.  \
        For example, if you want to compare a motif bed file \
        (given the column name "LFY1_motif") to -50/+100 bp from the summit \
        and a bed file with MNase nucleosome (given the column name \
        "High_MNase") to the actual summit site you would set \
        `--summitRegion LFY1_motif,50,100 High_MNase,0,0`')
    parser.add_argument('-pt', '--comparePeaksToText', nargs='*', help='list of \
        files that contain tables with headers. To compare peakwise the table must \
        include at least the following columns: chrom, start, stop', \
        required=False)
    parser.add_argument('-ptn', '--comparePeaksToTextNames', nargs='*', \
        help='list of prefixes for the text files that you want to compare with. \
        The number of values set to this parameter must be equal to \
        the number of values set to parameter `comparePeaksToText`. \
        (eg.dexVmock_db)', required=False, default=[])
    parser.add_argument('-ptf', '--addPeaksToTextFeatures', help='In each \
        peak-centric file given in `comparePeaksToText` there are other \
        columns that you may want to include in this analysis. If so you \
        need to set this in a specific format: \
        `[comparePeaksToTextName]:[column_of_interest_x],[column_of_interest_y],[...]`. \
        This data will be output in columns formatted: \
        `[comparePeaksToTextName]:[column_of_interest]`. For example, if we set \
        `--comparePeaksToTextNames timepoint1_dexVmock_db timepoint2_dexVmock_db \
        --addPeaksToTextFeatures timepoint1_dexVmock_db:logFC \
        timepoint2_dexVmock_db:logFC,adjP` then the following columns will be \
        written: `timepoint1_dexVmock_db:logFC`, `timepoint2_dexVmock_db:logFC`, \
        and `timepoint2_dexVmock_db:adjP`', required=False, default=[], \
        nargs='*')
    parser.add_argument('-gt', '--compareGenesToText', help='tab-delimited \
        text file(s) that contain gene information. Each row must contain a \
        column containing gene IDs that are labeled with the header "gene_id". \
        If this is set, `compareGenesToTextNames` and `addGeneToTextFeatures` \
        must also be set.', \
        required=False, default=[], nargs='*')
    parser.add_argument('-gtn', '--compareGenesToTextNames', help='prefix names for \
        each genewise text file given in `compareGenesToText`', \
        required=False, default=[], nargs='*')
    parser.add_argument('-gtf', '--addGeneToTextFeatures', help='Given the \
        gene-centric file(s) in `compareGenesToText` the script will report \
        True or False whether the genes annotated appear in the text file \
        ("gene-match") and/or report other columns from the table in \
        `compareGenesToText`. To specify the columns ([col1],[col2][..][coln]) \
        that are output set this parameter to this specific format: \
        `[compareGenesToTextNames]:[col1],[col2][..][coln]`. \
        To report the True/False "gene-match" set one of the columns to the \
        dash symbol "-". This data will be output in columns labeled by the \
        `compareGenesToTextName". For example, if this parameter is set: \
        `--addGeneToTextFeatures RPM:sample1,sample2 \
        Sample1vsSample2_DE:-,logFC,adjp` then the output will contain columns: \
        "RPM:sample1", "RPM:sample2", "Sample1vsSample2_DE" (which gives the \
        gene-match), "Sample1vsSample2_DE:logFC", and \
        "Sample1vsSample2_DE:adjp"', required=False, default=[], nargs='*')
    parser.add_argument('-df', '--dnase_files', nargs='*', help='list of bed \
        files with locations of DNase Hypersensitivity sites', required=False)
    parser.add_argument('-dn', '--dnase_names', nargs='*', help='list giving \
        names for each DNAse experiment corresponding to DNAse bed files in \
        --dnase_files. This list must be of equal length to the "dnase_files" \
        variable (eg.DNAse_flowers)', required=False, default=[])
    parser.add_argument('-idr', '--globalIDRpval', help='global IDR pvalue \
        used on to get this output', default="")
    parser.add_argument('--keep_tmps', help='keep temp files', \
        action='store_true')
    parser.add_argument('--verbose', help='echo processing', \
        action='store_true')

    args = parser.parse_args()
    if args.rnaDESignificance:
        if args.rnaDESignificance[1].lower() != "max" and \
            args.rnaDESignificance[1].lower() != "min":
            err_msg = ("ERROR in `rnaDESignificance` parameter format." + \
                "The second value should be max or min. Three values are required.")
            sys.exit(err_msg)
        else:
            args.rnaDESignificance[1] = args.rnaDESignificance[1].lower()
    if args.comparePeaksToBeds or args.comparePeaksToBedsNames:
        if len(args.comparePeaksToBeds) != len(args.comparePeaksToBedsNames):
            err_msg = "comparePeaksToBeds and comparePeaksToBedsNames must " + \
                "be of equal length!"
            err_msg = ("%s\ncomparePeaksToBeds length = %i" % (err_msg, \
                len(args.comparePeaksToBeds)))
            err_msg = ("%s\ncomparePeaksToBedsNames length = %i\n" % (err_msg, \
                len(args.comparePeaksToBedsNames)))
            sys.exit(err_msg)
    if args.comparePeaksToText or args.comparePeaksToTextNames:
        if len(args.comparePeaksToText) != len(args.comparePeaksToTextNames):
            err_msg = "comparePeaksToText and comparePeaksToTextNames must " + \
                "be of equal length!"
            err_msg = ("%s\ncomparePeaksToText length = %i" % (err_msg, \
                len(args.comparePeaksToText)))
            err_msg = ("%s\ncomparePeaksToTextNames length = %i\n" % (err_msg, \
                len(args.comparePeaksToTextNames)))
            sys.exit(err_msg)

    addPeaksToTextFeatures_dic = defaultdict(list) 
    if len(args.addPeaksToTextFeatures) > 0:
        for chip_feat in args.addPeaksToTextFeatures:
            prefixNCols =chip_feat.split(":")
            colList = prefixNCols[1].split(",")
            addPeaksToTextFeatures_dic[prefixNCols[0]]=colList
    if (len(args.compareGenesToText) != len(args.compareGenesToTextNames) or
        len(args.compareGenesToText) != len(args.addGeneToTextFeatures)):
            length1= len(args.compareGenesToText)
            length2= len(args.compareGenesToTextNames)
            length3= len(args.addGeneToTextFeatures)
            err_msg = (("ERROR: To perform genewise comparisons with " +
                "tab-delimited text files, the user must set all three " + 
                "parameters: `compareGenesToText`(user-set:%i), " +
                "`compareGenesToTextNames` (user-set:%i), and " + 
                "`addGeneToTextFeatures` (user-set:%i)") \
                 % (length1, length2,length3))
            sys.exit(err_msg)
    # addGeneToTextFeatures_dic[compareGenesToTextNames]=[col1,col2...,col3]
    addGeneToTextFeatures_dic = defaultdict(list) 
    if len(args.addGeneToTextFeatures) > 0:
        for chip_feat in args.addGeneToTextFeatures:
            prefixNCols =chip_feat.split(":")
            colList = prefixNCols[1].split(",")
            addGeneToTextFeatures_dic[prefixNCols[0]]=colList

    if args.compareSumRegToBed and not args.narrowpeak_file:
        err_msg = ("ERROR: A narrowPeak file must be set in the " +
            "`narrowpeak_file` parameter to perform the " + 
            "`compareSumRegToBed` function. Currently it is not set.")
    if args.compareSumRegToBed and not args.summitRegion:
        err_msg = ("ERROR: To perform the `compareSumRegToBed` function, " +
            "the user must set the `summitRegion` parameter. Currently " +
            "it is not set")
    if args.compareSumRegToBed or args.compareSumRegToBedColNames:
        if len(args.compareSumRegToBed) != len(args.compareSumRegToBedColNames):
            err_msg = (("ERROR: Parameters `compareSumRegToBed` and " + 
                "`compareSumRegToBedColNames` must be of equal length!\n" +
                "compareSumRegToBed length = %i\n" +
                "compareSumRegToBedColNames length = %i\n") % 
                (len(args.compareSumRegToBed),
                len(args.compareSumRegToBedColNames)))
            sys.exit(err_msg)
    if args.compareSumRegToBed or args.summitRegion:
        if len(args.compareSumRegToBed) != len(args.summitRegion):
            err_msg = (("ERROR: Parameters `compareSumRegToBed` and " + 
                "`summitRegion` must be of equal length!\n" +
                "compareSumRegToBed length = %i\n" +
                "summitRegion length = %i\n") % 
                (len(args.compareSumRegToBed),
                len(args.summitRegion)))
            sys.exit(err_msg)
    if args.dnase_files or args.dnase_names:
        if len(args.dnase_files) != len(args.dnase_names):
            err_msg = "dnase_files and dnase_names must be of equal length!"
            err_msg = ("%s\ndnase_files length = %i" % (err_msg, \
                len(args.dnase_files)))
            err_msg = ("%s\ndnase_names length = %i\n" % (err_msg, \
                len(args.dnase_names)))
            sys.exit(err_msg)
    summit_region_downstream={}
    summit_region_upstream={}
    if args.summitRegion:
        for sr in args.summitRegion:
            sr_err = (('ERROR: In correct format for parameter ' +
                '`summitRegion`. %s is written in correctly. ' +
                'The specific format for this parameter is: ' + 
                '`[compareSumRegToBedColNames],[bp_upstream],[bp_downstream]`. ' +
                'For example, if you want to compare a motif bed file ' +
                '(given the column name "LFY1_motif") to -50/+100 bp from ' + 
                'the summit and a bed file with MNase nucleosome ' + 
                '(given the column name "High_MNase") to the actual ' + 
                'summit site you would set ' + 
                '`--summitRegion LFY1_motif,50,100 High_MNase,0,0`')
                % (sr))
            if not sr.find(","):
                sys.exit(sr_err)
            summitRegArgList = sr.split(",")
            if len(summitRegArgList) != 3:
                sys.exit(sr_err)
            summit_region_downstream[summitRegArgList[0]]=summitRegArgList[1]
            summit_region_upstream[summitRegArgList[0]]=summitRegArgList[1]


    dir_name = os.path.abspath(args.dir_name)

    count_file = ("%s/%s_counts.txt" % (dir_name, args.prefix))

    counts_list = list()

    # import gene bed file
    gene_bedfile_df_cols = ["gene_chr", "gene_start", "gene_stop", \
        "gene_id", "gene_score", "gene_strand"]
    gene_bedfile_df = pd.read_csv(args.gene_bedfile, sep='\t', header=None, \
        names=gene_bedfile_df_cols)
    gene_bedfile_types={"gene_chr" : object, \
        "gene_start" : np.int64, "gene_stop" : np.int64, "gene_id" : object, \
        "gene_score" : np.float64, "gene_strand":object}
    gene_bedfile_df.loc[gene_bedfile_df["gene_score"]==".","gene_score"]=0
    gene_bedfile_df = gene_bedfile_df.astype(gene_bedfile_types)

    # roundOfAnnotation = {} # roundOfAnnotation[peak_name] = round_annotated
    # 1 = round 1 ()
    # 2 = round 2 (adopted peaks)
    # 3 = round 3 (orphan peaks)



    # ROUND 0: Characterize Peaks (independent of Gene Annotation)
    
    # get peak dataframe
    peak_info_columns = ["chr", "start", "stop", "name", "signal", "strand", \
                         "fold_change", "pValue", "qValue", "summit"]
    bed_cols = peak_info_columns[:6]
    peaks_df = pd.read_csv(args.bed_file, sep='\t', header=None, \
        names=bed_cols, dtype={"chr" : object, "start" : np.int64, \
        "stop" : np.int64, "name" : object, "signal" : np.float64, \
        "strand":object})

    ## print number of total peaks
    counts_list.append("number of total peaks\t%i" % len(peaks_df))
    
    all_peak_names = list(peaks_df["name"])

    # print average length of peaks
    length_series = peaks_df.loc[:, "stop"]- peaks_df.loc[:, "start"]
    counts_list.append("average length of peaks\t%f" % length_series.mean())

    cur_peakfile = args.bed_file
    if args.narrowpeak_file:
        cur_peakfile = args.narrowpeak_file
        narrowPeak_df = pd.read_csv(args.narrowpeak_file, 
            sep="\t", header=None, \
            index_col=False, names = peak_info_columns, \
            dtype={"chr" : object, "start" : np.int64, \
            "stop" : np.int64, "name" : object, "signal" : np.float64, \
            "strand" : object, "fold_change": np.float64, \
            "pValue" : np.float64, "qValue" : np.float64, \
            "summit" : np.int64})
        peaks_df = narrowPeak_df.copy()


    # compare experiments ChIP peaks with other bed files
    comparePeaksToBeds_colnames=args.comparePeaksToBedsNames
    if args.comparePeaksToBeds:
        for oldpeakidx in range(0, len(args.comparePeaksToBeds)):
            old_peak_file = args.comparePeaksToBeds[oldpeakidx]  # bed file to compare with
            old_peak_prefix = args.comparePeaksToBedsNames[oldpeakidx]  # prefix for bed file to compare with
            overlap_df_x = genome_locations.compare_bedfiles(args.bed_file, \
                old_peak_file, verbal=args.verbose)
            # count number of peaks that overlap this other ChIP experiment
            overlap_count_x = len(overlap_df_x["name"].unique())
            counts_list.append("Number of peaks overlap %s\t%i" % \
                               (old_peak_prefix, overlap_count_x))
            newPeaksInOldPeaks = list(overlap_df_x["name"])
            if args.compBedScores:
                peakb_feat_dic = {"chrom" : "chr_b", \
                    "chromStart" : "start_b", "chromEnd" : "stop_b", \
                    "name": "name_b", "score": "signal_b", \
                    "strand": "strand_b", "signalValue" : "fold_change_b", \
                    "pValue" : "pValue_b", "qValue" : "qValue_b", 
                    "peak": "summit_b"}
                if not peakb_feat_dic[args.compBedScores] in list(overlap_df_x.columns):
                    score_col_unavailable = (("\nWARNING: %s column is not available in %s." + \
                        " You may want to remove/edit your --peakScore parameter" +
                        " or use a different file.\n") % (args.compBedScores, old_peak_file))
                    sys.stderr.write(score_col_unavailable)
                else:
                    peakBScoresInPeakA_df =overlap_df_x.groupby(bed_cols)
                    peakBScores_series = peakBScoresInPeakA_df.apply( \
                        lambda x: ";".join(str(s) for s in list(x[peakb_feat_dic[args.compBedScores]])))
                    peakBScores_df = peakBScores_series.to_frame().reset_index()
                    newcolname=(("%s:%s") % (old_peak_prefix, args.compBedScores))
                    comparePeaksToBeds_colnames.append(newcolname)
                    peakBScores_df = \
                        peakBScores_df.rename(columns={0: (newcolname)})
                    peaks_df = peaks_df.merge(peakBScores_df, how="left", on=bed_cols)
                    peaks_df.loc[:,newcolname] = \
                        peaks_df[newcolname].fillna("False")
            peaks_df.loc[:,old_peak_prefix] = \
                np.where(peaks_df["name"].isin(newPeaksInOldPeaks),True, False)
            peaks_df.loc[:,old_peak_prefix] = \
                peaks_df[old_peak_prefix].fillna(False)
        overlap_all_counts = len(peaks_df.loc[peaks_df.apply( \
            lambda row: all(row[args.comparePeaksToBedsNames]), axis=1), \
                                              "name"].unique())
        counts_list.append("Number of peaks overlap all old peaks\t%i" % \
                           (overlap_all_counts))
    comparePeaksToText_colnames=args.comparePeaksToTextNames
    # compare experiments ChIP peaks with other peakwise text-files
    if args.comparePeaksToText:
        for textf_index in range(0, len(args.comparePeaksToText)):
            textf = args.comparePeaksToText[textf_index]  # text file to compare with
            textf_prefix = args.comparePeaksToTextNames[textf_index]  # prefix for text file to compare with
            # convert text file to a pandas dataframe
            textf_df = pd.read_csv(textf,sep='\t', header=0, dtype=str)
            if "chrom" not in list(textf_df.columns) or \
                "start" not in list(textf_df.columns) or \
                "stop" not in list(textf_df.columns):
                sys.exit(("ERROR: files listed in `comparePeaksToText` parameters " +
                    "must be tab-delimited tables that contain the following columns: " +
                    "`chrom`, `start`, and `stop`. This is not found in the following " +
                    "file: %s\n") % textf)
            
            bedColsInTextDf=["chrom","start","stop"]
            add=1
            removeBeforeMerge=[]
            for i in range(3,len(bed_cols)):
                if bed_cols[i] in list(textf_df.columns):
                    if add==1:
                        bedColsInTextDf.append(bed_cols[i])
                    else:
                        removeBeforeMerge.append(bed_cols[i])
                else:
                    add=0
            # change column names in text file to match names in comparison
            # text file
            colsList_for_file2 = ["chr_b", "start_b", "stop_b", "name_b", "signal_b", \
                "strand_b", "fold_change_b", "pValue_b", "qValue_b", "summit_b"]
            renameBedColsInTextDf_dic={}
            renameBedColsInTextDf_list=[]
            for i in range(0,len(bedColsInTextDf)):
                renameBedColsInTextDf_dic[bedColsInTextDf[i]]=colsList_for_file2[i]
                renameBedColsInTextDf_list.append(colsList_for_file2[i])
            textf_df = textf_df.rename(columns=renameBedColsInTextDf_dic)

            # set column types of textf_df
            col_dtype_dict = {"chr_b" : object, "start_b" : np.int64, \
                "stop_b" : np.int64, "name_b" : object, "signal_b" : np.float64, \
                "strand_b" : object, "fold_change_b": np.float64, \
                "pValue_b" : np.float64, "qValue_b" : np.float64, \
                "summit_b" : np.int64, "overlap":np.int64}
            for textdf_col in renameBedColsInTextDf_list:
                textf_df[textdf_col] = \
                    textf_df[textdf_col].astype(col_dtype_dict[textdf_col])


            comp_bed_df = textf_df.loc[:,renameBedColsInTextDf_list].copy()
            temp_bed = (("%s/%s.bed.tmp") % (dir_name,textf_prefix))
            comp_bed_df.to_csv(temp_bed, sep="\t", header=False, \
                    index=False)
            overlap_df_x = genome_locations.compare_bedfiles(args.bed_file, \
                temp_bed, verbal=args.verbose)
            if not args.keep_tmps:
                os.remove(temp_bed)
            # count number of peaks that overlap this other ChIP experiment
            overlap_count_x = len(overlap_df_x["name"].unique())
            counts_list.append("Number of peaks overlap %s\t%i" % \
                               (textf_prefix, overlap_count_x))
            newPeaksInOldPeaks = list(overlap_df_x["name"])
            if args.addPeaksToTextFeatures:
                
                if len(removeBeforeMerge) >0 :
                    textf_df = textf_df.drop(removeBeforeMerge, axis=1)
                
                colsList_for_file2 = colsList_for_file2[0:len(renameBedColsInTextDf_list)]
                overlap_df_x = overlap_df_x.merge(textf_df, how="left", 
                    left_on=colsList_for_file2,
                    right_on=renameBedColsInTextDf_list)
                for textColOfInterest in addPeaksToTextFeatures_dic[textf_prefix]:

                    if not textColOfInterest in list(overlap_df_x.columns):
                        score_col_unavailable = (("\nERROR: %s column is not available in %s." + \
                            " You may want to remove/edit your --addPeaksToTextFeatures parameter" +
                            " or use a different file.\n") % (textColOfInterest, textf))
                        sys.exit(score_col_unavailable)
                    peakBScoresInPeakA_df =overlap_df_x.groupby(bed_cols)
                    peakBScores_series = peakBScoresInPeakA_df.apply( \
                        lambda x: ";".join(str(s) for s in list(x[textColOfInterest])))
                    peakBScores_df = peakBScores_series.to_frame().reset_index()
                    newcolname=(("%s:%s") % (textf_prefix, textColOfInterest))
                    comparePeaksToText_colnames.append(newcolname)
                    peakBScores_df = \
                        peakBScores_df.rename(columns={0: (newcolname)})
                    peaks_df = peaks_df.merge(peakBScores_df, how="left", on=bed_cols)
                    peaks_df.loc[:,newcolname] = \
                    peaks_df[newcolname].fillna("False")
                
            peaks_df.loc[:,textf_prefix] = \
                np.where(peaks_df["name"].isin(newPeaksInOldPeaks),True, False)
            peaks_df.loc[:,textf_prefix] = \
                peaks_df[textf_prefix].fillna(False)
        overlap_all_counts = len(peaks_df.loc[peaks_df.apply( \
            lambda row: all(row[args.comparePeaksToBedsNames]), axis=1), \
                                              "name"].unique())
        counts_list.append("Number of peaks overlap all old peaks\t%i" % \
                           (overlap_all_counts))


    # compare bed files to summit region
    # summitRegion, compareSumRegToBed, compareSumRegToBedColNames
    if args.compareSumRegToBed:
        

        for m_idx in range(0, len(args.compareSumRegToBed)):
            mFile = args.compareSumRegToBed[m_idx]
            mPrefix = args.compareSumRegToBedColNames[m_idx]

            # set the summit region
            bp_downstream_of_summit=summit_region_downstream[mPrefix]
            bp_upstream_of_summit=summit_region_upstream[mPrefix]

            # create a bed file with each feature being a summit region
            summitRegionBedFile = ("%s/%s_summitPlus%iMinus%i.txt" % (dir_name, \
                args.prefix, bp_downstream_of_summit, bp_upstream_of_summit))
            summit_df = peaks_df.copy()
            summit_df["start"]=(peaks_df["start"] + 
                peaks_df["summit"] - bp_downstream_of_summit)
            summit_df["stop"]=(peaks_df["start"] + 
                peaks_df["summit"] + bp_upstream_of_summit + 1)
            summit_df = summit_df.loc[:,bed_cols]
            summit_df.to_csv(summitRegionBedFile, sep="\t", header=False, \
                index=False)


            bedtools_table = genome_locations.compare_bedfiles(summitRegionBedFile, \
                mFile, verbal=args.verbose)

            bedtools_table["location"] = bedtools_table["chr_b"].map(str) + "_" \
                                      + bedtools_table["start_b"].map(str) + "_" \
                                      + bedtools_table["stop_b"].map(str)
            featuresInSummitRegion_df = bedtools_table.groupby(bed_cols)
            bedLoc_series = featuresInSummitRegion_df.apply( \
                lambda x: ";".join(str(s) for s in list(x["location"])))
            bedLoc_df = bedLoc_series.to_frame().reset_index()
            bedLoc_df = bedLoc_df.rename(columns={0: (mPrefix)})
            bedLoc_df = bedLoc_df.loc[:,["name", mPrefix]]
            peaks_df = peaks_df.merge(bedLoc_df, how="left", on="name")
            if not args.keep_tmps:
                os.remove(summitRegionBedFile)


    # add dnase info to peakwise files
    dnase_peakwise_colnames=args.dnase_names
    if args.dnase_files:
        for d_idx in range(0, len(args.dnase_files)):
            dFile = args.dnase_files[d_idx]
            dPrefix = args.dnase_names[d_idx]
            dnase_table = genome_locations.compare_bedfiles(args.bed_file, \
                dFile, verbal=args.verbose)
            #dnase_table_cols_needed = list(range(0, 6)) + [dnase_table.shape[1] - 1]
            #dnase_table = dnase_table.iloc[:, dnase_table_cols_needed]
            peaksInDHS = list(dnase_table["name"].unique())
            peaks_df.loc[:, dPrefix] = peaks_df.apply( \
                lambda row: row["name"] in peaksInDHS, axis=1)
            peaks_df.loc[:,dPrefix].fillna(False, inplace=True)
            peakswithDHS_count = \
                len(dnase_table.loc[:, bed_cols].drop_duplicates())
            counts_list.append("number of peaks within DNAse sites (%s) \t%i" \
                               % (dPrefix, peakswithDHS_count))
            
            # count what percent of the peak contains DNAse regions
            dnase_table.astype(
                {'start':np.float64,
                'stop':np.float64,
                'overlap':np.float64})
            dnase_table["percOverlap"] = \
                dnase_table["overlap"] / (
                    dnase_table["stop"] - dnase_table["start"])
            dnase_table_grouped =dnase_table.groupby(bed_cols)
            peakPercOverlap_series = dnase_table_grouped.apply( \
                lambda x: x["percOverlap"].sum())
            peakPercOverlap_df = peakPercOverlap_series.to_frame().reset_index()
            newcolname=(("%s:percentDNAseOverlapsPeak") % (dPrefix))
            dnase_peakwise_colnames.append(newcolname)
            peakPercOverlap_df = \
                peakPercOverlap_df.rename(columns={0: (newcolname)})
            peaks_df = peaks_df.merge(peakPercOverlap_df, how="left", on=bed_cols)
            peaks_df.loc[:,newcolname] = \
            peaks_df[newcolname].fillna("False")

    # ROUND 1

    ## annotate peaks that are :
    ##		intragenic
    ##		${noFilter_tss_upstream} bp upstream
    ##		${noFilter_TTS_downstream} bp downstream
    round1_peaks = round1_annotation.r1_annotate('peak',
        args.gene_bedfile, args.bed_file, 
        peaks_df, args.prefix, dir_name,
        per_inter_filter=args.percent_inter,
        bp_upstream_filter=args.filter_tss_upstream,
        bp_downstream_filter=args.filter_tts_downstream,
        ignore_conv_peaks=args.ignore_conv_peaks,
        verbose=args.verbose, keep_tmps=args.keep_tmps)
    round1_peaks["distance_from_gene"] = \
        round1_peaks["distance_from_gene"].astype('float64')
    count_intragenic_genes = len( \
        round1_peaks.loc[round1_peaks["distance_from_gene"] == 0.0, \
                         bed_cols].drop_duplicates())
    count_upstream_genes_r1 = len( \
        round1_peaks.loc[round1_peaks["distance_from_gene"] < 0, \
                         bed_cols].drop_duplicates())
    count_downstream_genes_r1 = len( \
        round1_peaks.loc[round1_peaks["distance_from_gene"] > 0, \
                         bed_cols].drop_duplicates())

    ## simultaneously annotate summits
    ### Note that `round1_peaks` has a row for each
    ### peak-gene pairing. Therefore peaks with multiple
    ### gene annotations will be annotated in multiple rows.
    if args.narrowpeak_file:
        round1_summits = round1_annotation.r1_annotate('summit',
            args.gene_bedfile, args.narrowpeak_file, peaks_df, 
            args.prefix, dir_name, \
            per_inter_filter=args.percent_inter, \
            bp_upstream_filter=args.filter_tss_upstream, \
            bp_downstream_filter=args.filter_tts_downstream, \
            ignore_conv_peaks=args.ignore_conv_peaks, 
            verbose=args.verbose, keep_tmps=args.keep_tmps)
        round1_summits = round1_summits.loc[:, (
            bed_cols + gene_bedfile_df_cols + ["summit_start"])]

        round1_peaks = round1_peaks.merge(
            round1_summits, on=bed_cols+ gene_bedfile_df_cols,
            how='left')
        
        round1_peaks["summit_ann"] = round1_peaks.apply(
            lambda row: (row['start'] + 
                row["summit"]) == 
                row["summit_start"], axis=1)
        round1_peaks = round1_peaks.drop(["summit_start"],axis=1)

    if args.verbose:
        # number intragenic
        sys.stdout.write("Number of Peaks INSIDE Genes: %i\n" % \
                         count_intragenic_genes)
        # number upstream
        sys.stdout.write("Number of Peaks 1-%ibp UPSTREAM of Genes: %i\n" % \
                         (args.filter_tss_upstream, count_upstream_genes_r1))
        if args.filter_tts_downstream != 0:
            # number downstream
            sys.stdout.write("Number of Peaks 1-%ibp DOWNSTREAM of Genes: %i\n" \
                % (args.filter_tts_downstream, count_downstream_genes_r1))

    ### print Number of Peaks INSIDE Genes
    counts_list.append("peaks INSIDE genes in round 1:\t%i" % \
        count_intragenic_genes)
    ### print Number of Peaks Upstream of Genes
    counts_list.append("peaks 1-%ibp UPSTREAM of genes in round 1:\t%i" % \
                       (args.filter_tss_upstream, count_upstream_genes_r1))
    ### print Number of Peaks Downstream of Genes
    if args.filter_tts_downstream > 0:
        counts_list.append("peaks 1-%ibp DOWNSTREAM of genes in round 1:\t%i" \
                           % (args.filter_tts_downstream, count_downstream_genes_r1))
    round1_peaks["roundOfAnnotation"] = 1

    ### Print number of peaks annotated in round 1
    round1_count = len(round1_peaks.loc[:, bed_cols].drop_duplicates())
    counts_list.append("peaks annotated in round1\t%i" % \
                       (round1_count))

    upstream_convergent_peaks_file = ("%s/%s_convergent_upstream_peaks.txt" \
        % (dir_name,args.prefix))
    if (not os.path.isfile(upstream_convergent_peaks_file)):
        counts_list.append("number of UPSTREAM convergent peaks in round 1:\t0")
    else:
        up_convergent_peaks_df = pd.read_csv(upstream_convergent_peaks_file, \
            sep='\t', dtype=str)
        up_conv_peak_count = len(up_convergent_peaks_df)/2
        counts_list.append(("number of UPSTREAM convergent" + \
            " peaks in round 1:\t%i") % up_conv_peak_count)

    if args.filter_tts_downstream > 0:
        downstream_convergent_peaks_file = \
            ("%s/%s_convergent_downstream_peaks.txt" % (dir_name,args.prefix))
        if (not os.path.isfile(downstream_convergent_peaks_file)):
            counts_list.append("number of DOWNSTREAM convergent peaks in round 1:\t0")
        else:
            down_convergent_peaks_df = pd.read_csv( \
                downstream_convergent_peaks_file, sep='\t', dtype=str)
            down_conv_peak_count = len(down_convergent_peaks_df)/2
            counts_list.append(("number of DOWNSTREAM convergent peaks " + \
                "ignored in round 1:\t%i") % down_conv_peak_count)

    # echo "ROUND 2"
    all_peaks_df = pd.DataFrame()

    # get peaks not found in round1_peaks
    outlier_df = peaks_df[~peaks_df["name"].isin(round1_peaks["name"])]

    # compare gene locations to outliers

    # annotate peaks that have not been annotated yet to DE genes within a bp limit
    #### HERE I AM
    round2_peaks = pd.DataFrame(columns=list(peaks_df.columns))
    orphan_peaks_df = pd.DataFrame(columns=list(peaks_df.columns))
    all_de_genes = []
    if args.round2ann:
        # RNAseqs_dict = defaultdict(list) # RNAseqs_dict[sampleName]=[genes]
        for rna_samp in range(0, len(args.round2ann)):
            rna_diffExp_forRound2Ann_file = args.round2ann[rna_samp]
            rna_diffExp_forRound2Ann_df = pd.read_csv(rna_diffExp_forRound2Ann_file, sep='\t')
            if "gene_id" not in list(rna_diffExp_forRound2Ann_df.columns):
                rna_cols_str = ",".join(list(rna_diffExp_forRound2Ann_df.columns))
                geneIDerr= (("\nThe `gene_id` column cannot be found in %s." + \
                    "\nPlease check the column names.\n" + \
                    "COLUMNS FOUND: %s\n") % (rna_diffExp_forRound2Ann_file, rna_cols_str))
                sys.exit(geneIDerr)
            else:
                if len(rna_diffExp_forRound2Ann_df["gene_id"]) != \
                    len(rna_diffExp_forRound2Ann_df["gene_id"].unique()):
                    geneIDerr= (("\nSomething is wrong!\n" + \
                        "The `gene_id` column must be unique in %s.") % \
                    (rna_diffExp_forRound2Ann_file))
                    sys.exit(geneIDerr)
            
            de_genes = list(rna_diffExp_forRound2Ann_df.loc[:,"gene_id"])
            if args.rnaDESignificance:
                rna_de_sig_col=args.rnaDESignificance[0]
                rna_de_sig_val=np.float64(args.rnaDESignificance[2])
                if rna_de_sig_col in list(rna_diffExp_forRound2Ann_df.columns):
                    if args.rnaDESignificance[1] == "max":

                        de_genes=list(rna_diffExp_forRound2Ann_df.loc[ \
                            rna_diffExp_forRound2Ann_df[rna_de_sig_col]<=rna_de_sig_val,"gene_id"])
                    else:
                        de_genes=list(rna_diffExp_forRound2Ann_df.loc[ \
                            rna_diffExp_forRound2Ann_df[rna_de_sig_col]>=rna_de_sig_val,"gene_id"])
                else:
                    rnaScoreerr= (("\nThe %s column cannot be found " + \
                        "in %s.\nPlease check the column names.\n" + \
                        "COLUMNS FOUND: %s\n") % \
                        (rna_de_sig_col, rna_diffExp_forRound2Ann_file, list(rna_diffExp_forRound2Ann_df.columns)))
                    sys.exit(rnaScoreerr)

            # RNAseqs_dict[rna_sample_name] = de_genes
            all_de_genes += de_genes
        all_de_genes = list(set(all_de_genes))
        de_genes_df = gene_bedfile_df[gene_bedfile_df['gene_id'].isin(all_de_genes)]
        if args.round2ann:
            round2_ann_file = ("%s/%s_r2_peak_annotations.txt" 
                % (dir_name, args.prefix))
            round2_peaks = round2_annotation.r2_annotate(de_genes_df, 
                outlier_df, args.outlier_filter, round2_ann_file)
            round2_peaks["roundOfAnnotation"] = 2
            orphan_peaks_df = \
                outlier_df[~outlier_df["name"].isin(round2_peaks["name"])]

            ### Print number of peaks annotated in round 2
            round2_count = len(round2_peaks.loc[:, bed_cols].drop_duplicates())
            counts_list.append("peaks annotated in round2\t%i" % \
                               (round2_count))
            ### Print number of peaks annotated in round 1 AND 2
            counts_list.append("total peaks annotated\t%i" \
                % (round1_count + round2_count))
            # get dataframe with all peak annotations including unassigned peaks
            all_peaks_frame = [round1_peaks, round2_peaks, orphan_peaks_df]
            #all_peaks_df = pd.concat(all_peaks_frame, sort=False)
        else:
            orphan_peaks_df = outlier_df
            all_peaks_frame = [round1_peaks, orphan_peaks_df]
        all_peaks_df = pd.concat(all_peaks_frame, sort=True)

    else:
        orphan_peaks_df = outlier_df
        # get dataframe with all peak annotations including unassigned peaks
        all_peaks_frame = [round1_peaks, orphan_peaks_df]
        all_peaks_df = pd.concat(all_peaks_frame, sort=True)
        #all_peaks_df = pd.concat(all_peaks_frame)

    # reorganize columns
    all_peaks_col_order = bed_cols + ['roundOfAnnotation']
    if args.narrowpeak_file:
        all_peaks_col_order += ['qValue', 'summit']
    all_peaks_col_order = all_peaks_col_order + ['gene_id']
    if args.narrowpeak_file:
        all_peaks_col_order += ['summit_ann'] 
    all_peaks_col_order = all_peaks_col_order + \
            ['gene_overlap', 'distance_from_gene'] + \
            comparePeaksToBeds_colnames + comparePeaksToText_colnames + \
            args.compareSumRegToBedColNames + dnase_peakwise_colnames
    

    all_peaks_df = all_peaks_df.loc[:, all_peaks_col_order]
    #all_peaks_df = all_peaks_df.rename(columns={'Alias': 'gene_name'})

    # ROUND 3: print out peaks that remain unassigned (orphan peaks)
    orphan_bed_f = ("%s/%s_unassigned.bed" % (dir_name, args.prefix))
    orphan_peaks_df.iloc[:, 0:6].to_csv(orphan_bed_f, sep="\t", header=False, \
        index=False)
    if args.narrowpeak_file:
        orphan_np_f = ("%s/%s_unassigned.narrowPeak" % (dir_name, args.prefix))
        orphan_peaks_df.iloc[:, 0:10].to_csv(orphan_np_f, sep="\t", \
            header=False, index=False)

   

    ## How many genes were annotated in each round (indepedent of DE)
    gene_count_r1 = \
        len(all_peaks_df.loc[all_peaks_df["roundOfAnnotation"] == 1, \
        "gene_id"].unique())
    counts_list.append("genes annotated in round1\t%i" % \
                       (gene_count_r1))

    # For DNAse peaks that overlap with intragenic peaks,
    # list what percentage of the gene does the DNAse region
    # overlap
    dnase_genewise_cols=[]
    if args.dnase_files:
        intragenic_peaks_df = all_peaks_df.loc[(
            all_peaks_df['distance_from_gene']==0),bed_cols+['gene_id']].copy()
        intragenic_peaks_df = \
            intragenic_peaks_df.merge(
                gene_bedfile_df, on='gene_id',
                how='left')
        for d_idx in range(0, len(args.dnase_files)):
            dFile = args.dnase_files[d_idx]
            dPrefix = args.dnase_names[d_idx]
            gene_dnase_df = genome_locations.compare_bedfiles(dFile, \
                args.gene_bedfile, verbal=args.verbose)
            gene_dnase_df['gene_overlap'] = gene_dnase_df.apply(lambda row:
                row['overlap'] / (row['stop_b']- row['start_b']),
                axis=1)
            gene_dnase_df = gene_dnase_df.loc[:,['name_b','gene_overlap']]
            gene_dnase_grouped =gene_dnase_df.groupby("name_b")
            genePercOverlap_series = gene_dnase_grouped.apply( \
                lambda x: x["gene_overlap"].sum())
            genePercOverlap_df = genePercOverlap_series.to_frame().reset_index()

            newcolname=(("%s:percentDNAseInGene") % (dPrefix))
            dnase_genewise_cols.append(newcolname)
            genePercOverlap_df.columns = ["gene_id",newcolname]
            all_peaks_df = all_peaks_df.merge(genePercOverlap_df, 
                on=['gene_id'], how='left')
            

    # Assign DE genes to peaks
    if args.round2ann:
        gene_count_r2 = \
            len(all_peaks_df.loc[all_peaks_df["roundOfAnnotation"] == 2, \
                "gene_id"].unique())
        counts_list.append("genes annotated in round2\t%i" % \
                           (gene_count_r2))
        counts_list.append("total genes annotated\t%i" % \
                           (gene_count_r1 + gene_count_r2))


    # Add info from other ChIP genewise annotations
    addGeneToTextFeatures_cols=[]
    if args.compareGenesToText:
        for other_genewise_idx in range(0,len(args.compareGenesToText)):
            other_genewise_file = args.compareGenesToText[other_genewise_idx]
            other_genewise_prefix = args.compareGenesToTextNames[other_genewise_idx]
            other_genewise_df = pd.read_csv(other_genewise_file, sep="\t",
                index_col=False)
            # Report TRUE/FALSE whether both samples annotate to the same
            # genes
            all_peaks_df.loc[:,other_genewise_prefix] = \
                all_peaks_df["gene_id"].isin(other_genewise_df["gene_id"])
            # Add other features of other ChIP if other ChIP annotates to 
            # the respective genes 
            if args.addGeneToTextFeatures:
                gene_feature_list=addGeneToTextFeatures_dic[other_genewise_prefix]
                if "-" in gene_feature_list:
                    # Report TRUE/FALSE whether both samples 
                    # annotate to the same genes
                    all_peaks_df.loc[:,other_genewise_prefix] = \
                        all_peaks_df["gene_id"].isin(
                        other_genewise_df["gene_id"])
                    addGeneToTextFeatures_cols.append(other_genewise_prefix)
                    gene_feature_list.remove("-")
                column_check = \
                    [x for x in gene_feature_list if x not in list(other_genewise_df.columns)]
                if len(column_check) > 1:
                    sys.exit(("\nERROR: At least one feature in " + 
                        "`addGeneToTextFeatures` is NOT found.\n" + 
                        "Missing features in %s: %s\n") %
                        (other_genewise_file, ", ".join(column_check)))
                subset_genewise_df = other_genewise_df.loc[:,["gene_id"] + \
                    gene_feature_list]
                ochip_feature_cols = [other_genewise_prefix + ":" + x \
                    for x in gene_feature_list]
                subset_genewise_df.columns =  ["gene_id"] + ochip_feature_cols
                all_peaks_df = all_peaks_df.merge(subset_genewise_df, \
                    how="left", on="gene_id")
                addGeneToTextFeatures_cols = addGeneToTextFeatures_cols + \
                ochip_feature_cols



    # Print out Peak-centric datatable...
    peak_group_cols = bed_cols + ['roundOfAnnotation']
    if args.narrowpeak_file:
        peak_group_cols += ['qValue', 'summit']
    peak_group_cols = peak_group_cols + \
                      comparePeaksToBeds_colnames + \
                      comparePeaksToText_colnames + \
                      args.compareSumRegToBedColNames + \
                      dnase_peakwise_colnames
    
    ## note: non-peak_centric columns = gene_id, summit_ann, 
    # Alias, distance_from_gene
    peak_grouped_df = all_peaks_df.groupby(bed_cols)
    peak_centric_cols_df = all_peaks_df.loc[:, peak_group_cols]
    
    
    ## list number of genes that annotate to peak
    peak_nGenes_series = peak_grouped_df.apply(lambda x: len(x["gene_id"].unique()))
    peak_nGenes_df = peak_nGenes_series.to_frame().reset_index()
    peak_nGenes_df.columns = bed_cols + ['numGenes']
    peak_ann_df = peak_centric_cols_df.merge(peak_nGenes_df, how="left", \
                                             on=bed_cols)


    ## list gene ids that annotate to peaks
    peak_gid_series = peak_grouped_df.apply(lambda x: ";".join(str(s) for s in list(x["gene_id"].unique())))
    peak_gid_df = peak_gid_series.to_frame().reset_index()
    peak_gid_df.columns = bed_cols + ['gene_id']
    peak_ann_df = peak_ann_df.merge(peak_gid_df, how='left', on=bed_cols)

    ## list percent of gene that overlaps with the peak
    peak_gOverlap_series = peak_grouped_df.apply(lambda x: ";".join(str(s) for s in list(x["gene_overlap"].unique())))
    peak_gOverlap_df = peak_gOverlap_series.to_frame().reset_index()
    peak_gOverlap_df.columns = bed_cols + ['gene_overlap']
    peak_ann_df = peak_ann_df.merge(peak_gOverlap_df, 
        how='outer', on=bed_cols)

    ## list distances of genes that annotate to peaks
    peak_dist_series = peak_grouped_df.apply(lambda x: ";".join(str(s) for s in list(x["distance_from_gene"].unique())))
    peak_dist_df = peak_dist_series.to_frame().reset_index()
    peak_dist_df.columns = bed_cols + ['distance_from_gene']
    peak_ann_df = peak_ann_df.merge(peak_dist_df, 
        how='outer', on=bed_cols)
    
    ## if applicable, list whether the summit also annotates the gene
    if args.narrowpeak_file:
        peak_summitAnn_series = peak_grouped_df.apply(
            lambda x: ";".join(
            str(s) for s in list(x.loc[:,bed_cols+["summit_ann"]].
            drop_duplicates()["summit_ann"])))
        peak_summitAnn_df = peak_summitAnn_series.to_frame().reset_index()
        peak_summitAnn_df.columns = bed_cols + ['summit_ann']
        peak_ann_df = peak_ann_df.merge(peak_summitAnn_df, 
            how='left', on=bed_cols)
    
    ## if applicable, list percentage that DNAse (which overlaps with
    ## intragenic peak) overlaps gene
    if args.dnase_files:
        for dnase_pg_col in dnase_genewise_cols:
            peak_dnase2Gene_series = peak_grouped_df.apply(lambda x: ";".join(
                str(s) for s in list(x[dnase_pg_col].unique())))
            peak_dnase2Gene_df = peak_dnase2Gene_series.to_frame().reset_index()
            peak_dnase2Gene_df.columns = bed_cols + [dnase_pg_col]
            peak_ann_df = peak_ann_df.merge(peak_dnase2Gene_df, how='left', \
                                            on=bed_cols)
            
    ## if applicable, list features of any other ChIP experiments that were 
    ## annotated to the same genes
    if len(addGeneToTextFeatures_cols)>0:
        for other_chip_feats in addGeneToTextFeatures_cols:
            other_peak_features_series = peak_grouped_df.apply(
                lambda x: ";".join(
                    str(s) for s in list(x.loc[:,bed_cols+[other_chip_feats]].
                drop_duplicates()[other_chip_feats])))
            other_peak_features_df = other_peak_features_series.to_frame().reset_index()
            other_peak_features_df.columns = bed_cols + [other_chip_feats]
            peak_ann_df = peak_ann_df.merge(other_peak_features_df, how='left', \
                                            on=bed_cols)



    ## reorganize data-table to show gene columns closer to peak columns and
    ## extra info towards the later columns
    peak2gene_info_cols = ['numGenes', 'gene_id'] + \
        ['gene_overlap', 'distance_from_gene', 'summit_ann']
    peak_col_order = peak_group_cols[0:7] + peak2gene_info_cols + \
        peak_group_cols[7:] + dnase_genewise_cols + \
        addGeneToTextFeatures_cols 

    peak_ann_df = peak_ann_df.loc[:, peak_col_order]
    pd.set_option('float_format', '{:.2f}'.format)
    peak_out_tsv = ("%s/%s_peakwise_ann.tsv" % (dir_name, args.prefix))
    peak_out_csv = ("%s/%s_peakwise_ann.csv" % (dir_name, args.prefix))
    
    

    peak_ann_df = peak_ann_df.drop_duplicates()

    peak_ann_df.to_csv(peak_out_tsv, sep="\t", index=False, na_rep="NA")
    peak_ann_df.to_csv(peak_out_csv, index=False, na_rep="NA")

    # Print out Gene-centric datatable...
    gene_group_cols = ["gene_id"] + \
        dnase_genewise_cols + \
        addGeneToTextFeatures_cols
    all_ann_peaks = all_peaks_df.loc[all_peaks_df["gene_id"].notnull(),:].copy()
    geneCentricColsOnly_df = all_ann_peaks.loc[:,gene_group_cols].drop_duplicates()

    gene_groups_df = all_ann_peaks.groupby("gene_id")
    ## list number of genes that annotate to peak
    gene_nPeaks_series = gene_groups_df.apply(lambda x: x.shape[0])

    gene_nPeaks_df = gene_nPeaks_series.to_frame().reset_index()
    gene_nPeaks_df.columns = ["gene_id", "numPeaks"]
    ## list peak that annotate to gene
    all_ann_peaks.loc[:,"peaks"] = all_ann_peaks["name"].map(str) + "=" + \
                            all_ann_peaks["chr"].map(str) + ":" \
                            + all_ann_peaks["start"].map(str) + "-" \
                            + all_ann_peaks["stop"].map(str)

    gene_pname_series = gene_groups_df.apply(lambda x: ";".join(str(s) for s in list(x["peaks"].unique())))
    gene_pname_df = gene_pname_series.to_frame().reset_index()
    gene_pname_df.columns = ['gene_id', 'peaks']
    gene_ann_df = gene_nPeaks_df.merge(gene_pname_df, how='outer', on='gene_id')
    # round of annotation
    gene_roundAnn_series = gene_groups_df.apply(lambda x: ";".join(str(s) for s in list(x['roundOfAnnotation'].unique())))
    gene_roundAnn_df = gene_roundAnn_series.to_frame().reset_index()
    gene_roundAnn_df.columns = ['gene_id', 'roundOfAnnotation']
    gene_ann_df = gene_ann_df.merge(gene_roundAnn_df, how='outer', on='gene_id')
    # gene_overlap
    gene_overlap_series = gene_groups_df.apply(lambda x: ";".join(str(s) for s in list(x['gene_overlap'].unique())))
    gene_overlap_df = gene_overlap_series.to_frame().reset_index()
    gene_overlap_df.columns = ['gene_id', 'gene_overlap']
    gene_ann_df = gene_ann_df.merge(gene_overlap_df, how='outer', on='gene_id')
    # distance
    gene_dist_series = gene_groups_df.apply(lambda x: ";".join(str(s) for s in list(x['distance_from_gene'].unique())))
    gene_dist_df = gene_dist_series.to_frame().reset_index()
    gene_dist_df.columns = ['gene_id', 'distance_to_peak']
    gene_ann_df = gene_ann_df.merge(gene_dist_df, how='outer', on='gene_id')
    # list qvalues of peaks that annotate to the gene
    if args.narrowpeak_file:
        gene_qval_series = gene_groups_df.apply(lambda x: ";".join(str(s) for s in list(x['qValue'].unique())))
        gene_qval_df = gene_qval_series.to_frame().reset_index()
        gene_qval_df.columns = ['gene_id', 'qValue']
        gene_ann_df = gene_ann_df.merge(gene_qval_df, how='outer', on='gene_id')
        gene_summitAnn_series = gene_groups_df.apply(lambda x: ";".join(str(s) for s in list(x['summit_ann'].unique())))
        gene_summitAnn_df = gene_summitAnn_series.to_frame().reset_index()
        gene_summitAnn_df.columns = ['gene_id', 'summit_ann']
        gene_ann_df = gene_ann_df.merge(gene_summitAnn_df, how='outer', on='gene_id')
    ## list True if any of the peaks that the gene is annotated to is found in a
    ## a different ChIP sample (comparePeaksToBeds)
    if args.comparePeaksToBeds:
        for peak_sample in comparePeaksToBeds_colnames:
            #gene_peakcomp_series = gene_groups_df.apply(lambda x: x[peak_sample].any() != False)
            #gene_peakcomp_series = gene_groups_df.apply(lambda x: ";".join([str(y) for y in x[peak_sample] if y != False and y != "False"]))
            gene_peakcomp_series = gene_groups_df.apply(lambda x:listOrFalse(x[peak_sample]))
            gene_peakcomp_df = gene_peakcomp_series.to_frame().reset_index()
            gene_peakcomp_df.columns = ['gene_id', peak_sample]
            gene_ann_df = gene_ann_df.merge(gene_peakcomp_df, how='left', \
                                            on='gene_id')
            gene_ann_df.loc[:,peak_sample].fillna(False, inplace=True)

    ## list True if any of the peaks that the gene is annotated to is found in a
    ## a different ChIP sample (comparePeaksToBeds)
    if args.comparePeaksToText:
        for peak_sample in comparePeaksToText_colnames:
            #gene_peakcomp_series = gene_groups_df.apply(lambda x: x[peak_sample].any() != False)
            #gene_peakcomp_series = gene_groups_df.apply(lambda x: ";".join([str(y) for y in x[peak_sample] if y != False and y != "False"]))
            gene_peakcomp_series = gene_groups_df.apply(lambda x:listOrFalse(x[peak_sample]))
            gene_peakcomp_df = gene_peakcomp_series.to_frame().reset_index()
            gene_peakcomp_df.columns = ['gene_id', peak_sample]
            gene_ann_df = gene_ann_df.merge(gene_peakcomp_df, how='left', \
                                            on='gene_id')
            gene_ann_df.loc[:,peak_sample].fillna(False, inplace=True)
            
    ## list True if any of the peaks that the gene is annotated to 
    ## contains a feature in their summit region (compareSumRegToBed)
    if args.compareSumRegToBed:
        for sumRegion_col in args.compareSumRegToBedColNames:
            gene_sumRegion_series = gene_groups_df.apply(lambda x: x[sumRegion_col].notnull().any())
            gene_sumRegion_df = gene_sumRegion_series.to_frame().reset_index()
            gene_sumRegion_df.columns = ['gene_id', sumRegion_col]
            gene_ann_df = gene_ann_df.merge(gene_sumRegion_df, how='left', \
                                            on='gene_id')
            gene_ann_df.loc[:,sumRegion_col].fillna("NA", inplace=True)
    ## list True if any of the peaks that the gene is annotated to is found in a
    ## a DHS (dnase_files)
    if args.dnase_files:
        for dhs_col in dnase_peakwise_colnames:
            if dhs_col in args.dnase_names:
                gene_dhs_series = gene_groups_df.apply(lambda x: x[dhs_col].any())
            else:
                gene_dhs_series = gene_groups_df.apply(
                    lambda x: ";".join(str(s) for s in list(x[dhs_col].unique())))
            gene_dhs_df = gene_dhs_series.to_frame().reset_index()
            gene_dhs_df.columns = ['gene_id', dhs_col]
            gene_ann_df = gene_ann_df.merge(gene_dhs_df, how='left', \
                                            on='gene_id')
            gene_ann_df.loc[:,dhs_col].fillna(False, inplace=True)

    gene_ann_df = geneCentricColsOnly_df.merge(gene_ann_df, how='left', on='gene_id')
    gene_out_tsv = ("%s/%s_genewise_ann.tsv" % (dir_name, args.prefix))
    gene_out_csv = ("%s/%s_genewise_ann.csv" % (dir_name, args.prefix))
    gene_ann_df = gene_ann_df.drop_duplicates()
    gene_ann_df.to_csv(gene_out_tsv, sep="\t", index=False, na_rep="NA")
    gene_ann_df.to_csv(gene_out_csv,  index=False, na_rep="NA")

    olog_file = open(count_file, 'w')
    for i in range(0, len(counts_list)):
        olog_file.write("%s\n" % (counts_list[i]))

    olog_file.close()

