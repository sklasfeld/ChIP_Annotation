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
pd.options.mode.chained_assignment = 'raise'


def cmd(cmd_str, speak):
    """print and run bash command"""
    if speak:
        sys.stdout.write("%s\n" % cmd_str)
        sys.stdout.flush()
    os.system(cmd_str)


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

    parser = argparse.ArgumentParser(description="This is the main script for \
        gene annotation of ChIP peaks. This script first annotates ChIP peaks \
        by overlapping them with other ChIP peaks, DNase, MNase, or motifs. \
        Then it annotate the peaks to genes. The first round annotates \
        peaks by location relative to genes. Priority is given to intragenic \
        peaks (limited by `percent_inter`) and then peaks upstream of genes \
        (limited by `filter_tss_upstream`) and then peaks downstream of genes \
        (limited by `filter_tts_downstream`). The second round is \
        dependent on genes that were found to be significantly differentially \
        expressed (DE). In round 2, peaks are annotated to upstream or \
        downstream DE genes (limited by `outlier_filter`). Last, text files \
        are output that are contain peak-centric and gene-centric annotations")
    parser.add_argument('prefix', help='prefix for output file(eg. LFY)')
    parser.add_argument('dir_name', help='directory where output will go')
    parser.add_argument('bed_file', help='bed file containing ChIP peaks. \
        Warning: peaks should have original peak names.')
    parser.add_argument('gene_alist', help="A file that contains gene IDs \
        and their aliases. The file must be tab-delimited and have \
        atleast 1 column with the label `gene_id` \
        See `Araport11/Araport11_gene_info.txt file as an example \
        of the format. This is the default annotation")
    parser.add_argument('gene_bedfile', help='bed file with gene locations')
    parser.add_argument('-n', '--narrowpeak_file', help='narrowPeak file \
        containing ChIP peaks.')
    parser.add_argument('-ga', '--gene_alist_cols', nargs='*', help='columns \
        in `gene_alist` file that will be used to annotate each of the \
        gene ids. (eg. gene_name)', required=False, default=[])
    parser.add_argument('-bp', '--bedtools_path', help='path to bedtools', \
        default="")
    parser.add_argument('-sp', '--samtools_path', help='path to samtools', \
        default="")
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
    parser.add_argument('-of', '--outlier_filter', required=False, type=int,
        help='If `--round2ann` is set, this is the \
        maximum bp distance upstream/downstream of genes of outlier peaks \
        (peaks not annotated \
        in round 1) (default:10000).', default=10000)
    parser.add_argument('-ps', '--comparePeaksToBeds', nargs='*', help='list of \
        bed files that you want to compare with the peaks in this sample.', \
        required=False)
    parser.add_argument('-pn', '--comparePeaksToBedsNames', nargs='*', \
        help='list of prefixes for the bed files that you want to compare with. \
        This must be equal to "comparePeaksToBeds" variable \
        (eg.LFY_seedlings)', required=False, default=[])
    parser.add_argument('-pc', '--compBedScores', help='if peak overlaps with \
        a row in one of the bed files in `comparePeaksToBeds` then report a \
        feature of the old peak. Otherwise put NA. Any \
        column in a narrowPeak file can be reported (see \
        https://genome.ucsc.edu/FAQ/FAQformat.html#format12 for description \
        of each column).', choices=["chrom", "chromStart", "chromEnd", "name", \
        "score", "strand", \
        "singalValue", "pValue", "qValue", "peak"])
    parser.add_argument('-tf', '--comparePeaksToText', nargs='*', help='list of \
        files that contain tables with headers. To compare peakwise the table must \
        include at least the following columns: chrom, start, stop', \
        required=False)
    parser.add_argument('-tn', '--comparePeaksToTextNames', nargs='*', \
        help='list of prefixes for the text files that you want to compare with. \
        The number of values set to this parameter must be equal to \
        the number of values set to parameter `comparePeaksToText`. \
        (eg.dexVmock_db)', required=False, default=[])
    parser.add_argument('-tc', '--addTextFeatures', help='In each \
        peak-centric file given in `comparePeaksToText` there are other \
        columns that you may want to include in this analysis. If so you \
        need to set this in a specific format: \
        `[comparePeaksToTextName]:[column_of_interest_x],[column_of_interest_y],[...]`. \
        This data will be output in columns formatted: \
        `[comparePeaksToTextName]:[column_of_interest]`. For example, if we set \
        `--comparePeaksToTextNames timepoint1_dexVmock_db timepoint2_dexVmock_db \
        --addTextFeatures timepoint1_dexVmock_db:logFC \
        timepoint2_dexVmock_db:logFC,adjP` then the following columns will be \
        written: `timepoint1_dexVmock_db:logFC`, `timepoint2_dexVmock_db:logFC`, \
        and `timepoint2_dexVmock_db:adjP`', required=False, default=[], \
        nargs='*')
    parser.add_argument('-oc', '--otherChipGeneAnn', help='tab-delimited \
        text file(s) that contain gene-wise annotation information of other \
        ChIP experiments to compare with genes annotated in this sample. \
        Must contain columns labeled `gene_id` and `peaks`. \
        If this is set, `otherChipPrefix` must also be set and then script \
        will output columns (TRUE/FALSE) if both experiments had peaks that \
        were assigned to the same gene. Instead of TRUE/FALSE, you can set \
        `--otherChipGeneName` to get the names/location of the peak in the \
        `otherChipGeneAnn` file that annotates to the experiment of interest', \
        required=False, default=[], nargs='*')
    parser.add_argument('-op', '--otherChipPrefix', help='prefix names for \
        each genewise file given in `otherChipGeneAnn`. (eg. set \
        `-op K27_ABA_4hr` for H3K27me3 ChIP treated for 4hr with ABA)', \
        required=False, default=[], nargs='*')
    parser.add_argument('-on', '--otherChipGeneName', help='Instead of \
        TRUE/FALSE, this get the names/location of the peak in the \
        `otherChipGeneAnn` file that annotates to the experiment of interest \
        or FALSE')
    parser.add_argument('-ot', '--otherChipFeatures', help='In each \
        gene-centric file given in `otherChipGeneAnn` there are other \
        columns that you may want to include in this analysis. If so you \
        need to set this in a specific format: \
        `otherChipPrefix`:`col1,col2..coln`. This data will be output in \
        columns labeled by the `otherChipPrefix` and previous file name. For \
        example, if set such that \
        `--otherChipFeatures K27_ABA_4hr:4hr_ABA_v_none_logFC,4hr_ABA_v_none_adjP \
        K27_mock_4hr:4hr_mock_v_none_logFC,4hr_mock_v_none_adjp` then the \
        output will include columns: K27_ABA_4hr:4hr_ABA_v_none_logFC, \
        K27_ABA_4hr:4hr_ABA_v_none_adjP, K27_mock_4hr:4hr_mock_v_none_logFC, \
        and K27_mock_4hr:4hr_mock_v_none_adjP', required=False, default=[], \
        nargs='*')
    parser.add_argument('-rf', '--RNAcounts', help='tab-delimited \
    text files that contains counts (estimated counts, TPM, ect) for each gene. \
    Each column represents a different RNA sample. Note that the code expects \
    the column for each row name has the column header `gene_id` rather than \
    being empty.', nargs='*')
    parser.add_argument('-ra', '--RNAcounts_suffix', help='By default, \
        this script just names the RNA count results in the final annotation \
        file the same names as the column names given in the `RNAcounts` file. \
        You may want to add a suffix to these names to distinguish them as \
        counts. (eg. `_rnaTPM`)', required=False)
    parser.add_argument('-rs', '--compareRNAdiffExp', nargs='*', \
        help='tab-delimited text files that contains "differentially \
        expressed" gene-ids in a columns with the header `gene_ids`. ', \
        required=False, default=[])
    parser.add_argument('-rn', '--compareRNAdiffExpNames', nargs='*', \
        help='list of prefixes for the RNA samples that you want to compare \
        with. This must be equal to "compareRNAdiffExp" variable \
        (eg.RNA_wtVmutant)', required=False, default=[])
    parser.add_argument('-rl', '--rnaDESignificance', nargs=3, \
        help='using the tab-delimited text files set by the `round2ann` parameter, \
        we can limit which genes are defined as differentially expressed. \
        Note that this will only affect which genes are used for round 2 annotation \
        since round 2 only annotates based on DE genes. It will not \
        affect the columns that show differential expression values. This \
        parameter requires 3 values:[column] [min/max] [value]. For example, \
        if you want to limit DE genes to genes output from DESeq2 with adjp<=0.01 \
        then set this to `--rnaDESignificance adjp max 0.01`. Another example is \
        if you want to limit DE genes to genes output from DESeq2 with \
        log2FoldChange>=1 then set this to `--rnaDESignificance log2FoldChange min 1`.', \
        required=False)
    parser.add_argument('-rc', '--rnaScores', help='Set this parameter to the \
        column name(s) you would like to report from each \
        `compareRNAdiffExp` files. If you want to print the same columns \
        from each file in `compareRNAdiffExp` then use a comma-delimited \
        list. For example, to print `log2FoldChange` and `adjp` from each file \
        use `--rnaScores log2FoldChange,adjp`. If you want to print different \
        columns from different file then make comma-delimited lists seperated \
        by spaces for each file. For example, `--compareRNAdiffExp file1.txt \
        file2.txt --rnaScores log2FoldChange,adjp log2fc` will report \
        `log2FoldChange` and `adjp` from file1.txt and `log2fc` from \
        file2.txt', nargs='*', default=[])
    parser.add_argument('-mf', '--motifFiles', nargs='*', help='list of bed \
        files with locations of TF motifs. WARNING: each row in the bed files \
        must have a unique name in the 4th column.', required=False)
    parser.add_argument('-mn', '--motifNames', nargs='*', help='list giving \
        names of the type of motifs given in each bed file in --motifFiles. \
        This list must be of equal length to the "motifFiles" variable \
        (eg.LFY1)', required=False, default=[])
    parser.add_argument('-ms', '--callMotifBySummit', nargs=2, type=np.int64, 
        help='If `narrowpeak_file` is defined, then call the \
        presence of motifs by using $D bp downstream and $U bp upstream of the \
        summit rather than the peak region. For example to use -100/+250 bp \
        from the summit you would set `--callMotifBySummit 100 250`')
    parser.add_argument('-df', '--dnase_files', nargs='*', help='list of bed \
        files with locations of DNase Hypersensitivity sites', required=False)
    parser.add_argument('-dn', '--dnase_names', nargs='*', help='list giving \
        names for each DNAse experiment corresponding to DNAse bed files in \
        --dnase_files. This list must be of equal length to the "dnase_files" \
        variable (eg.DNAse_flowers)', required=False, default=[])
    parser.add_argument('-nf', '--mnase_files', nargs='*', help='list of bed \
    files with locations of nucleosome binding sites (via MNase)', \
    required=False)
    parser.add_argument('-nn', '--mnase_names', nargs='*', help='list giving \
        names for each MNase experiment corresponding to MNase bed files in \
        --mnase_files. This list must be of equal length to the "mnase_files" \
        variable (eg.HighMock)',required=False, default=[])
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
    if args.compareRNAdiffExp or args.compareRNAdiffExpNames:
        if len(args.compareRNAdiffExp) != len(args.compareRNAdiffExpNames):
            err_msg = "`compareRNAdiffExp` and `compareRNAdiffExpNames` " + \
                "must be of equal length!"
            err_msg = ("%s\n`compareRNAdiffExp` length = %i" % (err_msg, \
                len(args.compareRNAdiffExp)))
            err_msg = ("%s\n`compareRNAdiffExpNames` length = %i\n" % \
                (err_msg, len(args.compareRNAdiffExpNames)))
            sys.exit(err_msg)
    rna_score_cols = {}
    compareRNAdiffExpCols = []
    compareRNAdiffExpCols_dict = {}
    if len(args.rnaScores) > 1 and args.compareRNAdiffExp:
        if len(args.compareRNAdiffExp) != len(args.rnaScores):
            err_msg = "If `rnaScores` is set, it must contain 1 item or be " + \
                "of equal length to `compareRNAdiffExp`!"
            err_msg = ("%s\n`compareRNAdiffExp` length = %i" % (err_msg, \
                len(args.compareRNAdiffExp)))
            err_msg = ("%s`rnaScores` length = %i\n" % (err_msg, \
                len(args.rnaScores)))
            sys.exit(err_msg)
        else:
            for rna_col in range(0,len(args.compareRNAdiffExp)):
                rna_filename = args.compareRNAdiffExp[rna_col]
                rna_prefix = args.compareRNAdiffExpNames[rna_col]
                rna_score_cols[rna_filename] = args.rnaScores[rna_col].split(",") 
                rna_sample_col_names = [rna_prefix + ":" + x \
                    for x in rna_score_cols[rna_filename]]
                compareRNAdiffExpCols = compareRNAdiffExpCols + rna_sample_col_names
                compareRNAdiffExpCols_dict[rna_filename]=rna_sample_col_names
    elif len(args.rnaScores) == 1 and args.compareRNAdiffExp:
        for rna_col in range(0,len(args.compareRNAdiffExp)):
                rna_filename = args.compareRNAdiffExp[rna_col]
                rna_prefix = args.compareRNAdiffExpNames[rna_col]
                rna_score_cols[rna_filename] = args.rnaScores[0].split(",") 
                rna_sample_col_names = [rna_prefix + ":" + x \
                    for x in rna_score_cols[rna_filename]]
                compareRNAdiffExpCols = compareRNAdiffExpCols + rna_sample_col_names
                compareRNAdiffExpCols_dict[rna_filename]=rna_sample_col_names
    else:
        if args.compareRNAdiffExp:
            compareRNAdiffExpCols = args.compareRNAdiffExpNames

    addTextFeatures_dic = defaultdict(list) 
    if len(args.addTextFeatures) > 0:
        for chip_feat in args.addTextFeatures:
            prefixNCols =chip_feat.split(":")
            colList = prefixNCols[1].split(",")
            addTextFeatures_dic[prefixNCols[0]]=colList

    # set paths 
    if len(args.bedtools_path) > 0:
        args.bedtools_path = os.path.expanduser(args.bedtools_path)+"/"
    if len(args.samtools_path) > 0:
        args.samtools_path = os.path.expanduser(args.bedtools_path)+"/"


    if len(args.otherChipGeneAnn)>0 and len(args.otherChipPrefix)>0:
        if len(args.otherChipGeneAnn) != len(args.otherChipPrefix):
            length1= len(args.otherChipGeneAnn)
            length2= len(args.otherChipPrefix)
            err_msg = (("ERROR: User must set a prefix in `otherChipPrefix` " + \
                "corresponding exactly to each file given in " + \
                "`otherChipGeneAnn`. Instead, user set %i files " + \
                "to `otherChipGeneAnn` and %i prefixes to `otherChipPrefix`.") \
                 % (length1, length2))
            sys.exit(err_msg)
    elif len(args.otherChipGeneAnn)>0 and len(args.otherChipPrefix)==0:
        err_msg = (("ERROR: User must set a prefix in `otherChipPrefix` " + \
                "corresponding exactly to each of the files given in " + \
                "`otherChipGeneAnn`: %s.") % (",".join(args.otherChipGeneAnn)))
        sys.exit(err_msg)
    # otherChipFeatures_dic[otherChipPrefix]=[col1,col2...,col3]
    otherChipFeatures_dic = defaultdict(list) 
    if len(args.otherChipFeatures) > 0:
        for chip_feat in args.otherChipFeatures:
            prefixNCols =chip_feat.split(":")
            colList = prefixNCols[1].split(",")
            otherChipFeatures_dic[prefixNCols[0]]=colList


    if args.motifFiles or args.motifNames:
        if len(args.motifFiles) != len(args.motifNames):
            err_msg = "motifFiles and motifNames must be of equal length!"
            err_msg = ("%s\nmotifFiles length = %i" % (err_msg, \
                len(args.motifFiles)))
            err_msg = ("%s\nmotifNames length = %i\n" % (err_msg, \
                len(args.motifNames)))
            sys.exit(err_msg)
    if args.dnase_files or args.dnase_names:
        if len(args.dnase_files) != len(args.dnase_names):
            err_msg = "dnase_files and dnase_names must be of equal length!"
            err_msg = ("%s\ndnase_files length = %i" % (err_msg, \
                len(args.dnase_files)))
            err_msg = ("%s\ndnase_names length = %i\n" % (err_msg, \
                len(args.dnase_names)))
            sys.exit(err_msg)
    if args.mnase_files or args.mnase_names:
        if len(args.mnase_files) != len(args.mnase_names):
            err_msg = "mnase_files and mnase_names must be of equal length!"
            err_msg = ("%smnase_files length = %i" % (err_msg, \
                len(args.mnase_files)))
            err_msg = ("%smnase_names length = %i\n" % (err_msg, \
                len(args.mnase_names)))
            sys.exit(err_msg)
    if len(args.gene_alist_cols)<1:
        sys.stderr.write("WARNING: no columns are being used to annotate " + \
            "each gene!\n")

    dir_name = os.path.abspath(args.dir_name)
    gene_alist = args.gene_alist

    count_file = ("%s/%s_counts.txt" % (dir_name, args.prefix))

    counts_list = list()

    # import gene bed file
    gene_bedfile_df_cols = ["gene_chr", "gene_start", "gene_stop", \
        "gene_id", "gene_score", "gene_strand"]
    gene_bedfile_df = pd.read_csv(args.gene_bedfile, sep='\t', header=None, \
        names=gene_bedfile_df_cols, dtype={"gene_chr" : object, \
        "gene_start" : np.int64, "gene_stop" : np.int64, "gene_id" : object, \
        "gene_score" : np.float64, "gene_strand":object})

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
            compare_peak_file = ("%s/%s_peaks_found_in_%s.txt" % \
                                 (dir_name, args.prefix, old_peak_prefix))
            overlap_df_x = genome_locations.compare_bedfiles(args.bed_file, \
                old_peak_file, compare_peak_file, verbal=args.verbose, \
                bedtools_path=args.bedtools_path)
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
            if not args.keep_tmps:
                os.remove(compare_peak_file)
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

            compare_peak_file = ("%s/%s_peaks_found_in_%s.txt" % \
                                 (dir_name, args.prefix, textf_prefix))
            overlap_df_x = genome_locations.compare_bedfiles(args.bed_file, \
                temp_bed, compare_peak_file, verbal=args.verbose, \
                bedtools_path=args.bedtools_path)
            if not args.keep_tmps:
                os.remove(temp_bed)
            # count number of peaks that overlap this other ChIP experiment
            overlap_count_x = len(overlap_df_x["name"].unique())
            counts_list.append("Number of peaks overlap %s\t%i" % \
                               (textf_prefix, overlap_count_x))
            newPeaksInOldPeaks = list(overlap_df_x["name"])
            if args.addTextFeatures:
                
                if len(removeBeforeMerge) >0 :
                    textf_df = textf_df.drop(removeBeforeMerge, axis=1)
                
                colsList_for_file2 = colsList_for_file2[0:len(renameBedColsInTextDf_list)]
                overlap_df_x = overlap_df_x.merge(textf_df, how="left", 
                    left_on=colsList_for_file2,
                    right_on=renameBedColsInTextDf_list)
                for textColOfInterest in addTextFeatures_dic[textf_prefix]:

                    if not textColOfInterest in list(overlap_df_x.columns):
                        score_col_unavailable = (("\nERROR: %s column is not available in %s." + \
                            " You may want to remove/edit your --addTextFeatures parameter" +
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
            if not args.keep_tmps:
                os.remove(compare_peak_file)
        overlap_all_counts = len(peaks_df.loc[peaks_df.apply( \
            lambda row: all(row[args.comparePeaksToBedsNames]), axis=1), \
                                              "name"].unique())
        counts_list.append("Number of peaks overlap all old peaks\t%i" % \
                           (overlap_all_counts))
    # compare experiments ChIP peaks with RNA DE experiments
    for rna_samp in range(0, len(args.compareRNAdiffExp)):
        rna_diffExp_file = args.compareRNAdiffExp[rna_samp]
        rna_diffExp_df = pd.read_csv(rna_diffExp_file, sep='\t')
        if "gene_id" not in list(rna_diffExp_df.columns):
            rna_cols_str = ",".join(list(rna_diffExp_df.columns))
            geneIDerr= (("\nThe `gene_id` column cannot be found in %s." + \
                "\nPlease check the column names.\n" + \
                "COLUMNS FOUND: %s\n") % (rna_diffExp_file, rna_cols_str))
            sys.exit(geneIDerr)
        else:
            if len(rna_diffExp_df["gene_id"]) != \
                len(rna_diffExp_df["gene_id"].unique()):
                geneIDerr= (("\nSomething is wrong!\n" + \
                    "The `gene_id` column must be unique in %s.") % \
                (rna_diffExp_file))
                sys.exit(geneIDerr)
        if len(args.rnaScores) > 0:
            for rna_sampl_col in rna_score_cols[rna_diffExp_file]:
                if rna_sampl_col not in list(rna_diffExp_df.columns):
                    rna_cols_str = ",".join(list(rna_diffExp_df.columns))
                    rnaScoreerr= (("\nThe %s column cannot be found " + \
                        "in %s.\nPlease check the column names.\n" + \
                        "COLUMNS FOUND: %s\n") % \
                        (rna_sampl_col, rna_diffExp_file, rna_cols_str))
                    sys.exit(rnaScoreerr) 

    # add motif info to peakwise files
    if args.motifFiles:
        motifSearchPeaks=""
        if args.callMotifBySummit and args.narrowpeak_file:
                bp_downstream_of_summit=args.callMotifBySummit[0]
                bp_upstream_of_summit=args.callMotifBySummit[1]
                motifSearchPeaks = ("%s/%s_summitPlus%iMinus%i.txt" % (dir_name, \
                    args.prefix, bp_downstream_of_summit, bp_upstream_of_summit))
                summit_df = peaks_df.copy()
                summit_df["start"]=peaks_df["start"] + peaks_df["summit"] - bp_downstream_of_summit
                summit_df["stop"]=peaks_df["start"] + peaks_df["summit"] + bp_upstream_of_summit
                summit_df = summit_df.loc[:,bed_cols]
                summit_df.to_csv(motifSearchPeaks, sep="\t", header=False, \
                    index=False)
        else:
            motifSearchPeaks=args.bed_file

        for m_idx in range(0, len(args.motifFiles)):
            mFile = args.motifFiles[m_idx]
            mPrefix = args.motifNames[m_idx]

            motif_info_file = ("%s/%s_peaks_in_motif_%s.txt" % (dir_name, \
                                                             args.prefix, \
                                                             mPrefix))
            motif_table = genome_locations.compare_bedfiles(motifSearchPeaks, \
                mFile, motif_info_file, verbal=args.verbose, \
                bedtools_path=args.bedtools_path)
            ### START HERE!!!
            motif_table["location"] = motif_table["chr_b"].map(str) + "_" \
                                      + motif_table["start_b"].map(str) + "_" \
                                      + motif_table["stop_b"].map(str)
            motifsInPeak_df = motif_table.groupby(bed_cols)
            motifsLoc_series = motifsInPeak_df.apply( \
                lambda x: ";".join(str(s) for s in list(x["location"])))
            motifsLoc_df = motifsLoc_series.to_frame().reset_index()
            motifsLoc_df = motifsLoc_df.rename(columns={0: (mPrefix)})
            motifsLoc_df = motifsLoc_df.loc[:,["name", mPrefix]]
            peaks_df = peaks_df.merge(motifsLoc_df, how="left", on="name")
            peakswithmotif_count = \
                len(motifsLoc_df.loc[:, bed_cols].drop_duplicates())
            counts_list.append("number of peaks with motif %s \t%i" \
                               % (mPrefix, peakswithmotif_count))
            if not (args.keep_tmps):
                os.remove(motif_info_file)
        if not (args.keep_tmps) and args.callMotifBySummit and args.narrowpeak_file:
            os.remove(motifSearchPeaks)


    # add dnase info to peakwise files
    dnase_peakwise_colnames=args.dnase_names
    if args.dnase_files:
        for d_idx in range(0, len(args.dnase_files)):
            dFile = args.dnase_files[d_idx]
            dPrefix = args.dnase_names[d_idx]
            dnase_info_file = ("%s/peaks_in_dnase_%s.txt" % (dir_name, \
                                                             dPrefix))
            dnase_table = genome_locations.compare_bedfiles(args.bed_file, \
                dFile, dnase_info_file, verbal=args.verbose, \
                bedtools_path=args.bedtools_path)
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
            
            # count what percent of the DNAse region overlaps the peak
            dnase_table.astype(
                {'start':np.float64,
                'stop':np.float64,
                'overlap':np.float64})
            dnase_table["percOverlap"] = \
                dnase_table["overlap"] / (
                    dnase_table["stop"] - dnase_table["start"])
            dnase_table_grouped =dnase_table.groupby(bed_cols)
            peakPercOverlap_series = dnase_table_grouped.apply( \
                lambda x: ";".join(str(s) for s in list(x["percOverlap"])))
            peakPercOverlap_df = peakPercOverlap_series.to_frame().reset_index()
            newcolname=(("%s:percentDNAseOverlapsPeak") % (dPrefix))
            dnase_peakwise_colnames.append(newcolname)
            peakPercOverlap_df = \
                peakPercOverlap_df.rename(columns={0: (newcolname)})
            peaks_df = peaks_df.merge(peakPercOverlap_df, how="left", on=bed_cols)
            peaks_df.loc[:,newcolname] = \
            peaks_df[newcolname].fillna("False")
            


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
            peaksInNucleosomes = list(mnase_table["name"].unique())
            peaks_df.loc[:, mnPrefix] = peaks_df.apply( \
                lambda row: row["name"] in peaksInNucleosomes, axis=1)
            peaks_df.loc[:,mnPrefix].fillna(False, inplace=True)
            peakswithNuc_count = \
                len(mnase_table.loc[:, bed_cols].drop_duplicates())
            counts_list.append("number of peak with summits within MNAse sites (%s) \t%i" \
                               % (mnPrefix, peakswithNuc_count))
            if not (args.keep_tmps):
                os.remove(mnase_info_file)

    # ROUND 1

    ## annotate peaks that are :
    ##		intragenic
    ##		${noFilter_tss_upstream} bp upstream
    ##		${noFilter_TTS_downstream} bp downstream
    round1_peaks = round1_annotation.r1_annotate('peak',
        gene_alist, args.gene_bedfile,
        args.bed_file, peaks_df, args.prefix, dir_name,
        gene_alist_cols=args.gene_alist_cols,
        per_inter_filter=args.percent_inter,
        bp_upstream_filter=args.filter_tss_upstream,
        bp_downstream_filter=args.filter_tts_downstream,
        ignore_conv_peaks=args.ignore_conv_peaks,
        bedtools_path=args.bedtools_path, verbose=args.verbose)


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
            gene_alist, args.gene_bedfile, \
            args.narrowpeak_file, peaks_df, args.prefix, dir_name, \
            gene_alist_cols=args.gene_alist_cols, \
            per_inter_filter=args.percent_inter, \
            bp_upstream_filter=args.filter_tss_upstream, \
            bp_downstream_filter=args.filter_tts_downstream, \
            ignore_conv_peaks=args.ignore_conv_peaks, \
            bedtools_path=args.bedtools_path, verbose=args.verbose)
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
            # rna_sample_name = args.compareRNAdiffExpNames[rna_samp]
            
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
            round2_ann_file = ("%s/%s_r2_peak_annotations.txt" % (dir_name, args.prefix))
            round2_peaks = round2_annotation.r2_annotate(gene_alist, de_genes_df, 
                outlier_df, args.outlier_filter, round2_ann_file, \
                gene_alist_cols=args.gene_alist_cols)
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
    all_peaks_col_order = all_peaks_col_order + ['gene_id'] + \
                          args.gene_alist_cols 
    if args.narrowpeak_file:
        all_peaks_col_order += ['summit_ann'] 
    all_peaks_col_order = all_peaks_col_order + \
            ['gene_overlap', 'distance_from_gene'] + \
            comparePeaksToBeds_colnames + comparePeaksToText_colnames + \
            args.motifNames + dnase_peakwise_colnames + args.mnase_names
    

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
    
    rna_counts_names=[]
    if args.RNAcounts:
        for rna_counts_file in args.RNAcounts:
            rna_counts_df = pd.read_csv(rna_counts_file, sep='\t')
            new_rna_counts_names = list(rna_counts_df.columns)[1:]
            if args.RNAcounts_suffix:
                new_rna_counts_names = [x  + args.RNAcounts_suffix \
                    for x in new_rna_counts_names]
                rna_counts_df.columns = ['gene_id'] + new_rna_counts_names
            rna_counts_names = rna_counts_names + new_rna_counts_names
            all_peaks_df = all_peaks_df.merge(rna_counts_df, on="gene_id", \
                how = "left")

    # For DNAse peaks that overlap with intragenic peaks,
    # list what percentage of the gene does the DNAse region
    # overlap
    dnase_genewise_cols=[]
    dnase_peakAndGenewise_cols=[]
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
            dnase_info_file = ("%s/peaks_in_dnase_%s.txt" % (dir_name, \
                                                             dPrefix))
            dnase_table = genome_locations.compare_bedfiles(args.bed_file, \
                dFile, dnase_info_file, verbal=args.verbose, \
                bedtools_path=args.bedtools_path)
            dnase_table = intragenic_peaks_df.merge(dnase_table, 
                on=bed_cols, how="left")
            
            dnase_table['gene_overlap'] = dnase_table.apply(
                lambda x: gene_overlap_dnase(x), 
                axis=1).astype(np.float64)
            dnase_table['perc_gene_overlap'] = dnase_table.apply( lambda row: 
                row['gene_overlap'] / (row['gene_stop']- row['gene_start']),
                axis=1)
            dnase_table = dnase_table.loc[:,bed_cols+['gene_id','perc_gene_overlap']].copy()
            newcolname=(("%s:IntraGeneOverlap") % (dPrefix))
            dnase_peakAndGenewise_cols.append(newcolname)
            dnase_table = dnase_table.rename(index=str, 
                columns={"perc_gene_overlap": newcolname})
            all_peaks_df = all_peaks_df.merge(dnase_table, 
                on=bed_cols+['gene_id'], how='left')


            gene_dnase_file = ("%s/peaks_in_dnase_%s.txt" % (dir_name, \
                                                             dPrefix))
            gene_dnase_df = genome_locations.compare_bedfiles(dFile, \
                args.gene_bedfile, gene_dnase_file, verbal=args.verbose, \
                bedtools_path=args.bedtools_path)
            gene_dnase_df['gene_overlap'] = gene_dnase_df.apply(lambda row:
                row['overlap'] / (row['stop_b']- row['start_b']),
                axis=1)
            gene_dnase_df = gene_dnase_df.loc[:,['name_b','gene_overlap']]
            newcolname=(("%s:percentDNAseInGene") % (dPrefix))
            dnase_genewise_cols.append(newcolname)
            gene_dnase_df.columns = ["gene_id",newcolname]
            all_peaks_df = all_peaks_df.merge(gene_dnase_df, 
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

    ## Which genes are DE according to the RNA analysis?
    for rna_samp in range(0, len(args.compareRNAdiffExp)):
        rna_diffExp_file = args.compareRNAdiffExp[rna_samp] # RNA DE file
        rna_diffExp_df = pd.read_csv(rna_diffExp_file, sep='\t')
        if len(args.rnaScores) > 0:
            # new column name for RNA DE file
            new_rna_diffExp_columns = compareRNAdiffExpCols_dict[rna_diffExp_file]
            old_rna_diffExp_columns = rna_score_cols[rna_diffExp_file]
            sub_rna_diffExp_df= rna_diffExp_df.loc[:,["gene_id"] + \
                old_rna_diffExp_columns]
            sub_rna_diffExp_df.columns = ["gene_id"] + \
                    new_rna_diffExp_columns
            all_peaks_df = all_peaks_df.merge(sub_rna_diffExp_df, \
                    on = "gene_id", how = "left")
        else:
            new_rna_diffExp_column = compareRNAdiffExpCols[rna_samp]
            de_genes = list(rna_diffExp_df.loc[:,"gene_id"].unique())
            all_peaks_df.loc[:, new_rna_diffExp_column] = all_peaks_df.apply( \
                lambda row: row["gene_id"] in de_genes, axis=1)
            #all_peaks_df.loc[:, compareRNAdiffExpCols].fillna(False, inplace=True)
        counts_list.append(("total differentially expressed genes " + \
        "in %s\t%i") % (args.compareRNAdiffExpNames, \
            len(rna_diffExp_df.loc[:,"gene_id"].unique())))
    if len(args.compareRNAdiffExp) > 0:
        peaks_ann_to_de_counts = len(all_peaks_df.loc[ \
            all_peaks_df[compareRNAdiffExpCols].notnull().sum(axis=1) > 0, "gene_id"].unique())
        round1_splice = all_peaks_df[all_peaks_df['roundOfAnnotation'] == 1]
        r1_peaks_ann_to_de_counts = len(round1_splice.loc[ \
            round1_splice[compareRNAdiffExpCols].notnull().sum(axis=1) > 0, "gene_id"].unique())
        
        counts_list.append("# genes annotated to sample in round1\t%i" % \
                           (r1_peaks_ann_to_de_counts))
        counts_list.append("# genes annotated to sample in both rounds\t%i" % \
                           (peaks_ann_to_de_counts))
        if len(args.rnaScores) == 0:
            all_peaks_df = all_peaks_df.loc[:, compareRNAdiffExpCols].fillna(False)
        de_round1 = len([x for x in \
            all_peaks_df.loc[( \
            (all_peaks_df[compareRNAdiffExpCols].any(axis=1)) & \
            (all_peaks_df[u'roundOfAnnotation'] == 1)), \
            "gene_id"].unique() if str(x) != 'nan'])
        counts_list.append(("# DE genes annotated given all RNA " + \
            "samples in round1\t%i") % de_round1)
        de_genes_total = len([x for x in \
            all_peaks_df.loc[all_peaks_df[compareRNAdiffExpCols].all(axis=1), \
            "gene_id"].unique() if str(x) != 'nan'])
        counts_list.append(("# DE genes annotated given all RNA samples " + \
            "in both rounds\t%i") % de_genes_total)

    # Add info from other ChIP genewise annotations
    ### HERE!
    otherChIPFeatures_cols=[]
    if args.otherChipGeneAnn:
        for other_genewise_idx in range(0,len(args.otherChipGeneAnn)):
            other_genewise_file = args.otherChipGeneAnn[other_genewise_idx]
            other_genewise_prefix = args.otherChipPrefix[other_genewise_idx]
            other_genewise_df = pd.read_csv(other_genewise_file, sep="\t")
            # Add location of other ChIP if other ChIP annotates to this gene, 
            # Otherwise False
            if args.otherChipGeneName:
                subset_genewise_df = other_genewise_df.loc[:,["gene_id","peaks"]]
                subset_genewise_df.rename(index=str, \
                    columns={"peaks": other_genewise_prefix})
                all_peaks_df = all_peaks_df.merge(subset_genewise_df, how="left", \
                                             on="gene_id")
            # Add True if other ChIP annotates to this gene, 
            # Otherwise False
            else:
                all_peaks_df.loc[:,other_genewise_prefix] = \
                    all_peaks_df["gene_id"].isin(other_genewise_df["gene_id"])
            # Add Other features of other ChIP if other ChIP annotates to this gene, 
            # Otherwise False
            if args.otherChipFeatures:
                subset_genewise_df = other_genewise_df.loc[:,["gene_id"] + \
                    otherChipFeatures_dic[other_genewise_prefix]]
                ochip_feature_cols = [other_genewise_prefix + ":" + x \
                    for x in otherChipFeatures_dic[other_genewise_prefix]]
                subset_genewise_df.columns =  ["gene_id"] + ochip_feature_cols
                all_peaks_df = all_peaks_df.merge(subset_genewise_df, \
                    how="left", on="gene_id")
                otherChIPFeatures_cols = otherChIPFeatures_cols + \
                ochip_feature_cols



    # Print out Peak-centric datatable...
    peak_group_cols = bed_cols + ['roundOfAnnotation']
    if args.narrowpeak_file:
        peak_group_cols += ['qValue', 'summit']
    peak_group_cols = peak_group_cols + \
                      comparePeaksToBeds_colnames + \
                      comparePeaksToText_colnames + \
                      args.motifNames + \
                      dnase_peakwise_colnames + \
                      args.mnase_names
    
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
    ##
    ## list gene names that annotate to peaks
    if len(args.gene_alist_cols)>0:
        for gene_ann_col in args.gene_alist_cols:
            peak_gname_series = peak_grouped_df.apply(lambda x: ";".join(str(s) for s in list(x[gene_ann_col].unique())))
            peak_gname_df = peak_gname_series.to_frame().reset_index()
            peak_gname_df.columns=bed_cols+[gene_ann_col]
            peak_ann_df = peak_ann_df.merge(peak_gname_df,how='outer',on=bed_cols)
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
    

    ## if applicable, list counts (eg.TPMs) of the gene expression for 
    ## each RNA sample
    if args.RNAcounts:
        for rna_name in rna_counts_names:
            peak_tpmGene_series = peak_grouped_df.apply(lambda x: ";".join(str(s) for s in list(x[rna_name].unique())))
            peak_tpmGene_df = peak_tpmGene_series.to_frame().reset_index()
            peak_tpmGene_df.columns = bed_cols + [rna_name]
            peak_ann_df = peak_ann_df.merge(peak_tpmGene_df, how='left', \
                                            on=bed_cols)
    
    ## if applicable, list percentage that DNAse (which overlaps with
    ## intragenic peak) overlaps gene
    if args.dnase_files:
        for dnase_pg_col in dnase_peakAndGenewise_cols + dnase_genewise_cols:
            peak_dnase2Gene_series = peak_grouped_df.apply(lambda x: ";".join(
                str(s) for s in list(x[dnase_pg_col].unique())))
            peak_dnase2Gene_df = peak_dnase2Gene_series.to_frame().reset_index()
            peak_dnase2Gene_df.columns = bed_cols + [dnase_pg_col]
            peak_ann_df = peak_ann_df.merge(peak_dnase2Gene_df, how='left', \
                                            on=bed_cols)
    ## if applicable, list True (or RNA DE value) if any of the genes 
    ## that the peak is annotated to is DE in each sample
    if args.compareRNAdiffExp:
        for rna_name in compareRNAdiffExpCols:
            peak_degene_series = peak_grouped_df.apply(lambda x: x[rna_name].any())
            if len(args.rnaScores) > 0:
                peak_degene_series = peak_grouped_df.apply( \
                    lambda x: ";".join(str(s) for s in list(x[rna_name].unique())))
            peak_degene_df = peak_degene_series.to_frame().reset_index()
            peak_degene_df.columns = bed_cols + [rna_name]
            peak_ann_df = peak_ann_df.merge(peak_degene_df, how='left', \
                                            on=bed_cols)
            peak_ann_df.loc[:,rna_name].fillna(False, inplace=True)
    ## if applicable, list True (or location) of any other ChIP experiments that were 
    ## annotated to the same genes
    if args.otherChipGeneAnn:
        for other_chip_name in args.otherChipPrefix:
            both_peaks_ann2Gene_series = peak_grouped_df.apply(lambda x: ";".join(str(s) for s in list(x[other_chip_name].unique())))
            both_peaks_ann2Gene_df = both_peaks_ann2Gene_series.to_frame().reset_index()
            both_peaks_ann2Gene_df.columns = bed_cols + [other_chip_name]
            peak_ann_df = peak_ann_df.merge(both_peaks_ann2Gene_df, how='left', \
                                            on=bed_cols)
            
    ## if applicable, list features of any other ChIP experiments that were 
    ## annotated to the same genes
    if len(otherChIPFeatures_cols)>0:
        for other_chip_feats in otherChIPFeatures_cols:
            other_peak_features_series = peak_grouped_df.apply(lambda x: ";".join(str(s) for s in list(x[other_chip_name].unique())))
            other_peak_features_df = other_peak_features_series.to_frame().reset_index()
            other_peak_features_df.columns = bed_cols + [other_chip_feats]
            peak_ann_df = peak_ann_df.merge(other_peak_features_df, how='left', \
                                            on=bed_cols)



    ## reorganize data-table to show gene columns closer to peak columns and
    ## extra info towards the later columns
    peak2gene_info_cols = ['numGenes', 'gene_id'] + args.gene_alist_cols + \
        ['gene_overlap', 'distance_from_gene', 'summit_ann']
    peak_col_order = peak_group_cols[0:7] + peak2gene_info_cols + \
        rna_counts_names + compareRNAdiffExpCols + \
        dnase_genewise_cols + peak_group_cols[7:] + \
        args.otherChipPrefix + otherChIPFeatures_cols \
        + dnase_peakAndGenewise_cols

    peak_ann_df = peak_ann_df.loc[:, peak_col_order]
    pd.set_option('float_format', '{:.2f}'.format)
    peak_out_tsv = ("%s/%s_peakwise_ann.tsv" % (dir_name, args.prefix))
    peak_out_csv = ("%s/%s_peakwise_ann.csv" % (dir_name, args.prefix))
    
    

    peak_ann_df = peak_ann_df.drop_duplicates()

    peak_ann_df.to_csv(peak_out_tsv, sep="\t", index=False, na_rep="NA")
    peak_ann_df.to_csv(peak_out_csv, index=False, na_rep="NA")

    # Print out Gene-centric datatable...
    gene_group_cols = ["gene_id"] + args.gene_alist_cols + \
        rna_counts_names + compareRNAdiffExpCols + \
        dnase_genewise_cols + args.otherChipPrefix + \
        otherChIPFeatures_cols
    all_ann_peaks = all_peaks_df.loc[all_peaks_df["gene_id"].notnull(),:].copy()
    geneCentricColsOnly_df = all_ann_peaks.loc[:,gene_group_cols].copy()
    

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
            
    ## list True if any of the peaks that the gene is annotated to is found in a
    ## a motif (motifFiles)
    if args.motifFiles:
        for motif_col in args.motifNames:
            gene_motif_series = gene_groups_df.apply(lambda x: x[motif_col].notnull().any())
            gene_motif_df = gene_motif_series.to_frame().reset_index()
            gene_motif_df.columns = ['gene_id', motif_col]
            gene_ann_df = gene_ann_df.merge(gene_motif_df, how='left', \
                                            on='gene_id')
            gene_ann_df.loc[:,motif_col].fillna("NA", inplace=True)
    ## list True if any of the peaks that the gene is annotated to is found in a
    ## a DHS (dnase_files)
    if args.dnase_files:
        for dhs_col in dnase_peakwise_colnames + dnase_peakAndGenewise_cols:
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
    ## list True if any of the peaks that the gene is annotated to is found in a
    ## a nucleosome (mnase_files)
    if args.mnase_files:
        for mnase_samp in args.mnase_names:
            gene_mnase_series = gene_groups_df.apply(lambda x: x[mnase_samp].any())
            gene_mnase_df = gene_mnase_series.to_frame().reset_index()
            gene_mnase_df.columns = ['gene_id', mnase_samp]
            gene_ann_df = gene_ann_df.merge(gene_mnase_df, how='left', \
                                            on='gene_id')
            gene_ann_df.loc[:,mnase_samp].fillna(False, inplace=True)

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

