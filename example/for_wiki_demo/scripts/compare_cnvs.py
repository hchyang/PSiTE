#!/usr/bin/env python
""" Script to compare two sets of CNV calls.


    Written 2015-2017 by Regina Bohnert, Nora Rieber

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.

"""

import sys
import os
import argparse
import subprocess
import pybedtools
import pandas as pd

def parse_args(argv):
    """Parse options from the command line.

    Parameters
    ----------
    argv : list
        List of command line arguments.

    Returns
    -------
    args : argparse object
        Parsed command line arguments.
    """

    parser = argparse.ArgumentParser(description='Compare two sets of CNV calls.')
    parser.add_argument('-e1', '--exp1', default='.', type=str, required=True, help='name of file with CNV calls from experiment 1', metavar='FILE', dest='fname1')
    parser.add_argument('-e2', '--exp2', default='.', type=str, required=True, help='name of file with CNV calls from experiment 2 / ground truth', metavar='FILE', dest='fname2')
    parser.add_argument('-g', '--genome', default='.', type=str, required=True, help='name of genome file (contig and corresponding contig size)', metavar='FILE', dest='genome')
    parser.add_argument('-e', '--evalgenebased', default='.', type=str, required=True, help='do a gene-based evaluation? y or n. (This requires the previous generation of a gene-based text file both for ground truth and caller results; see gene_centered_analysis.py)', metavar='FILE', dest='evalgenebased')
    args = parser.parse_args()
    return args

def Jaccard_stats(bed_fname1, bed_fname2, genome):
    """Compute Jaccard index.

    Parameters
    ----------
    bed_fname1 : string
        Name of file with CNV calls from experiment 1 in BED format.
    bed_fname2 : string
        Name of file with CNV calls from experiment 2 in BED format.
    genome : string
        Name of genome file.

    Returns
    -------
    jacc_idx : float
        Jaccard index.
    """

    exp1 = pybedtools.BedTool(bed_fname1).merge()
    exp2 = pybedtools.BedTool(bed_fname2).merge()
    res = exp1.jaccard(exp2, sorted=True, g=genome)
    pybedtools.cleanup(remove_all=True)
    jacc_idx = res['jaccard']
    return jacc_idx


def assess_intersect(intersect_bed_fname):
    """Given a BED file generated with eval_intervals(), compare the amp/del status of the intersected intervals.

    Parameters
    ----------
    intersect_bed_fname : string
        BED file with 8 columns as specified above

    Returns
    -------
    cnt : int
        number of unique -wa intervals with coherent amp/del status
    cnt_inverted : int
        number of unique -wa intervals with inverted amp/del status
    cnt_del : int
        number of unique -wa intervals with coherent amp/del status - deletions only
    cnt_amp : int
        number of unique -wa intervals with coherent amp/del status - amplifications only
    """

    cnt = 0
    cnt_inverted = 0
    cnt_del = 0
    cnt_amp = 0
    previouslineconcat = ''

    with open(intersect_bed_fname, 'r') as fid_in:
        for line in fid_in:
            line = line.strip().split('\t')

            assert(line[3] == 'amplification' or line[3] == 'deletion')
            assert(line[7] == 'amplification' or line[7] == 'deletion')

            currlineconcat = line[0] + '\t' + line[1] + '\t' + line[2] + '\t' + line[3]
            if currlineconcat != previouslineconcat:
                if line[3] == line[7]:
                    cnt += 1

                    if line[3] == 'amplification':
                        cnt_amp += 1
                    elif line[3] == 'deletion':
                        cnt_del += 1
                    else:
                        print 'Error: unexpected input. Fourth column is neither amplification nor deletion!'
                        assert(False)
                else:
                    cnt_inverted += 1
            previouslineconcat = currlineconcat

    return cnt, cnt_inverted, cnt_del, cnt_amp



def assess_intersect_bp(intersect_bed_fname):
    """Given a BED file generated with eval_nt(), compare the amp/del status of the intersected intervals.

    Parameters
    ----------
    intersect_bed_fname : string
        BED file with 8 columns as specified above.

    Returns
    -------
    cnt : int
        number of bp in the overlap of two intervals with coherent amp/del status
    cnt_inverted : int
        number of bp in the overlap of two intervals with inverted amp/del status
    cnt_del : int
        number of bp in the overlap of two intervals with coherent amp/del status - deletions only
    cnt_amp : int
        number of bp in the overlap of two intervals with coherent amp/del status - amplifications only
    """

    cnt = 0
    cnt_inverted = 0
    cnt_del = 0
    cnt_amp = 0

    with open(intersect_bed_fname, 'r') as fid_in:
        for line in fid_in:
            line = line.strip().split('\t')

            assert(line[3] == 'amplification' or line[3] == 'deletion')
            assert(line[7] == 'amplification' or line[7] == 'deletion')

            if line[3] == line[7]:
                cnt += int(line[8])

                if line[3] == 'amplification':
                    cnt_amp += int(line[8])
                elif line[3] == 'deletion':
                    cnt_del += int(line[8])
                else:
                    print 'Error: unexpected input. Fourth column is neither amplification nor deletion!'
                    assert(False)
            else:
                cnt_inverted += int(line[8])

    return cnt, cnt_inverted, cnt_del, cnt_amp


def eval_intervals(pred_fname, gt_fname, genome, frac):
    """Evaluate with respect to intervals.

    Parameters
    ----------
    pred_fname : string
        Name of file with predicted CNV calls in BED format (chr, start, stop, amplification/deletion).
    gt_fname : string
        Name of file with ground truth CNV calls in BED format (chr, start, stop, amplification/deletion).
    genome : string
        Name of genome file.
    frac : float
        Fraction of interval overlap.

    Returns
    -------
    tp : int
        Number of true intervals.
    pp : int
        Number of predicted intervals.
    td : int
        Number of true discovered intervals.
    gtp : int
        Number of ground truth intervals.
    pr : float
        Interval precision defined as TP/(TP+FP).
    rc : float
        Interval recall defined as TP/(TP+FN).
    pr_inverted : float
        Interval precision without requiring matching ampdel_status
    rc_inverted : float
        Interval recall without requiring matching ampdel_status
    pr_amp : float
        Interval precision for amplifications only
    rc_amp : float
        Interval recall for amplifications only
    pr_del : float
        Interval precision for deletions only
    rc_del : float
        Interval recall for deletions only
    pp_amp : int
        Number of predicted amplifications
    pp_del : int
        Number of predicted deletions
    """

    pred = pd.read_table(pred_fname,header=None,index_col=(3),names=['contig','start','stop','ampdel_status'])
    pp = len(pred.index)
    pp_amp = len(pred.query('ampdel_status == \'amplification\''))
    pp_del = len(pred.query('ampdel_status == \'deletion\''))

    gt = pd.read_table(gt_fname,header=None,index_col=(3),names=['contig','start','stop','ampdel_status'])
    gtp = len(gt.index)
    gtp_amp = len(gt.query('ampdel_status == \'amplification\''))
    gtp_del = len(gt.query('ampdel_status == \'deletion\''))

    bedtools_intersect_fname1 = pred_fname.replace('.bed', '_'+str(frac)+'_pred_intersected_w_groundtruth.bed')
    bedtools_intersect_cmd = 'bedtools intersect -a ' + pred_fname + ' -b ' + gt_fname + ' -wa -wb -f ' + str(frac) + ' -g ' + genome + ' -sorted > ' + bedtools_intersect_fname1
    subprocess.check_output(bedtools_intersect_cmd, shell=True)
    tp, tp_inverted, tp_del, tp_amp = assess_intersect(bedtools_intersect_fname1)

    bedtools_intersect_fname2 = pred_fname.replace('.bed', '_'+str(frac)+'_groundtruth_intersected_w_pred.bed')
    bedtools_intersect_cmd = 'bedtools intersect -a ' + gt_fname + ' -b ' + pred_fname + ' -wa -wb -f ' + str(frac) + ' -g ' + genome + ' -sorted > ' + bedtools_intersect_fname2
    subprocess.check_output(bedtools_intersect_cmd, shell=True)
    td, td_inverted, td_del, td_amp = assess_intersect(bedtools_intersect_fname2)

    if pp > 0:
        pr = tp / float(pp)
        pr_inverted = (tp + tp_inverted) / float(pp)
    else:
        pr = 0.0
        pr_inverted = 0.0

    if gtp > 0:
        rc = td / float(gtp)
        rc_inverted = (td + td_inverted) / float(gtp)
    else:
        rc = 0.0
        rc_inverted = 0.0

    if pp_amp > 0:
        pr_amp = tp_amp / float(pp_amp)
    else:
        pr_amp = 0.0

    if gtp_amp > 0:
        rc_amp = td_amp / float(gtp_amp)
    else:
        rc_amp = 0.0

    if pp_del > 0:
        pr_del = tp_del / float(pp_del)
    else:
        pr_del = 0.0

    if gtp_del > 0:
        rc_del = td_del / float(gtp_del)
    else:
        rc_del = 0.0

    return tp, pp, td, gtp, pr, rc, pr_inverted, rc_inverted, pr_amp, rc_amp, pr_del, rc_del, pp_amp, pp_del


def eval_nt(pred_fname, gt_fname, genome):
    """Evaluate with respect to nucleotides.

    Parameters
    ----------
    pred_fname : string
        Name of file with predicted CNV calls in BED format (chr, start, stop, amplification/deletion)
    gt_fname : string
        Name of file with ground truth CNV calls in BED format (chr, start, stop, amplification/deletion)
    genome : string
        Name of genome file.

    Returns
    -------
    tp : int
        Number of true nucleotides.
    pp : int
        Number of predicted nucleotides.
    gtp : int
        Number of ground truth nucleotides.
    pr : float
        Nucleotide precision defined as TP/(TP+FP).
    rc : float
        Nucleotide recall defined as TP/(TP+FN).
    pr_amp : float
        Nucleotide precision for amplifications only
    rc_amp : float
        Nucleotide recall for amplifications only
    pr_del : float
        Nucleotide precision for deletions only
    rc_del : float
        Nucleotide recall for deletions only
    """

    pred = pd.read_table(pred_fname,header=None,index_col=(3),names=['contig','start','stop','ampdel_status'])
    pred_amp = pred.query('ampdel_status == \'amplification\'')
    pred_del = pred.query('ampdel_status == \'deletion\'')
    pp = sum(pred['stop']-pred['start'])
    pp_amp = sum(pred_amp['stop']-pred_amp['start'])
    pp_del = sum(pred_del['stop']-pred_del['start'])

    gt = pd.read_table(gt_fname,header=None,index_col=(3),names=['contig','start','stop','ampdel_status'])
    gt_amp = gt.query('ampdel_status == \'amplification\'')
    gt_del = gt.query('ampdel_status == \'deletion\'')
    gtp = sum(gt['stop']-gt['start'])
    gtp_amp = sum(gt_amp['stop']-gt_amp['start'])
    gtp_del = sum(gt_del['stop']-gt_del['start'])

    bedtools_intersect_fname1 = pred_fname.replace('.bed', '_pred_intersected_w_groundtruth_nominoverlap_bpoverlap.bed')
    bedtools_intersect_cmd = 'bedtools intersect -a ' + pred_fname + ' -b ' + gt_fname + ' -wo -g ' + genome + ' -sorted > ' + bedtools_intersect_fname1
    subprocess.check_output(bedtools_intersect_cmd, shell=True)
    tp, tp_inverted, tp_del, tp_amp = assess_intersect_bp(bedtools_intersect_fname1)

    if pp > 0:
        pr = tp / float(pp)
        pr_inverted = (tp + tp_inverted) / float(pp)
    else:
        pr = 0.0
        pr_inverted = 0.0
    if gtp > 0:
        rc = tp / float(gtp)
        rc_inverted = (tp + tp_inverted) / float(gtp)
    else:
        rc = 0.0
        rc_inverted = 0.0

    if pp_amp > 0:
        pr_amp = tp_amp / float(pp_amp)
    else:
        pr_amp = 0.0
    if gtp_amp > 0:
        rc_amp = tp_amp / float(gtp_amp)
    else:
        rc_amp = 0.0

    if pp_del > 0:
        pr_del = tp_del / float(pp_del)
    else:
        pr_del = 0.0
    if gtp_del > 0:
        rc_del = tp_del / float(gtp_del)
    else:
        rc_del = 0.0

    return tp, pp, gtp, pr, rc, pr_inverted, rc_inverted, pr_amp, rc_amp, pr_del, rc_del


def eval_gene_based(pred_fname, gt_fname, frac):
    """Evaluate with a gene-centered view.

    Parameters
    ----------
    pred_fname : string
        Name of gene-based file with predicted CNV calls, as generated by gene_centered_analysis.py
    gt_fname : string
        Name of gene-based file with ground truth CNV calls, as generated by gene_centered_analysis.py
    frac : float
        fraction of overlap.

    Returns
    -------
    tp : int
        Number of genes with over n% overlap (both in ground truth and predicted CNVs).
    fp : int
        Number of genes with over n% overlap in predicted CNVs but not in ground thruth CNVs
    fn : int
        Number of genes with over n% overlap in ground truth CNVs but not in predicted CNVs
    pr : float
        Gene-based precision defined as TP/(TP+FP).
    rc : float
        Gene-based recall defined as TP/(TP+FN).
    pr_inverted : float
        Gene-based precision without considering the amp/del status during the comparison to ground truth
    rc_inverted : float
        Gene-based recall without considering the amp/del status during the comparison to ground truth
    pr_amp : float
        Gene-based precision - amplifications only
    rc_amp : float
        Gene-based recall - amplifications only
    pr_del : float
        Gene-based precision - deletions only
    rc_del : float
        Gene-based recall - deletions only
    """

    frac_thresh = frac * 100
    tp = 0
    tp2 = 0
    fp = 0
    fn = 0

    tp_inverted = 0
    tp2_inverted = 0
    fp_inverted = 0
    fn_inverted = 0

    tp_amp = 0
    tp2_amp = 0
    fp_amp = 0
    fn_amp = 0

    tp_del = 0
    tp2_del = 0
    fp_del = 0
    fn_del = 0

    pred_CNVs_gene_based = pd.read_table(pred_fname,header=None,index_col=(0,2),names=['gene_name','percentage','ampdel_status'])
    gt_CNVs_gene_based = pd.read_table(gt_fname,header=None,index_col=(0,2),names=['gene_name','percentage','ampdel_status'])

    pred_CNVs_gene_based_inverted = pd.read_table(pred_fname,header=None,index_col=(0),names=['gene_name','percentage','ampdel_status'])
    gt_CNVs_gene_based_inverted = pd.read_table(gt_fname,header=None,index_col=(0),names=['gene_name','percentage','ampdel_status'])

    assert(pd.Index.nunique(pred_CNVs_gene_based.index) == len(pred_CNVs_gene_based.index))
    assert(pd.Index.nunique(gt_CNVs_gene_based.index) == len(gt_CNVs_gene_based.index))

    subtable_pred_CNVs = pred_CNVs_gene_based.query('percentage >= @frac_thresh')
    subtable_gt_CNVs = gt_CNVs_gene_based.query('percentage >= @frac_thresh')
    subtable_pred_CNVs_inverted = pred_CNVs_gene_based_inverted.query('percentage >= @frac_thresh')
    subtable_gt_CNVs_inverted = gt_CNVs_gene_based_inverted.query('percentage >= @frac_thresh')

    for index_gt in subtable_gt_CNVs.index:
        if index_gt in subtable_pred_CNVs.index:
            tp += 1
            if 'deletion' in index_gt:
                tp_del += 1
            elif 'amplification' in index_gt:
                tp_amp += 1
            else:
                print 'WARN: neither deletion nor amplification in index_gt (tp).'
                assert(False)
        else:
            fn += 1
            if 'deletion' in index_gt:
                fn_del += 1
            elif 'amplification' in index_gt:
                fn_amp += 1
            else:
                print 'WARN: neither deletion nor amplification in index_gt (fn).'
                assert(False)

    for index_pred in subtable_pred_CNVs.index:
        if index_pred in subtable_gt_CNVs.index:
            tp2 += 1
            if 'deletion' in index_pred:
                tp2_del += 1
            elif 'amplification' in index_pred:
                tp2_amp += 1
            else:
                print 'WARN: neither deletion nor amplification in index_gt (tp2).'
                assert(False)
        else:
            fp += 1
            if 'deletion' in index_pred:
                fp_del += 1
            elif 'amplification' in index_pred:
                fp_amp += 1
            else:
                print 'WARN: neither deletion nor amplification in index_gt (fp).'
                assert(False)

    assert(tp == tp2)
    assert(tp_amp == tp2_amp)
    assert(tp_del == tp2_del)

    for index_gt in subtable_gt_CNVs_inverted.index:
        if index_gt in subtable_pred_CNVs_inverted.index:
            tp_inverted += 1
        else:
            fn_inverted += 1
    for index_pred in subtable_pred_CNVs_inverted.index:
        if index_pred in subtable_gt_CNVs_inverted.index:
            tp2_inverted += 1
        else:
            fp_inverted += 1
    tp_inverted = min(tp_inverted, tp2_inverted)

    if (tp+fp) > 0:
        pr = float(tp)/(tp+fp)
    else:
        pr = 0.0
    if (tp+fn) > 0:
        rc = float(tp)/(tp+fn)
    else:
        rc = 0.0

    if (tp_inverted+fp_inverted) > 0:
        pr_inverted = float(tp_inverted)/(tp_inverted+fp_inverted)
    else:
        pr_inverted = 0.0
    if (tp_inverted+fn_inverted) > 0:
        rc_inverted = float(tp_inverted)/(tp_inverted+fn_inverted)
    else:
        rc_inverted = 0.0

    if (tp_amp+fp_amp) > 0:
        pr_amp = float(tp_amp)/(tp_amp+fp_amp)
    else:
        pr_amp = 0.0
    if (tp_amp+fn_amp) > 0:
        rc_amp = float(tp_amp)/(tp_amp+fn_amp)
    else:
        rc_amp = 0.0

    if (tp_del+fp_del) > 0:
        pr_del = float(tp_del)/(tp_del+fp_del)
    else:
        pr_del = 0.0
    if (tp_del+fn_del) > 0:
        rc_del = float(tp_del)/(tp_del+fn_del)
    else:
        rc_del = 0.0

    print frac_thresh
    print tp, fp, fn, pr, rc, pr_inverted, rc_inverted, pr_amp, rc_amp, pr_del, rc_del

    return tp, fp, fn, pr, rc, pr_inverted, rc_inverted, pr_amp, rc_amp, pr_del, rc_del


def compare_cnv_calls(fname1, fname2, genome, evalgenebased, frac=0.8):
    """Compare two sets of CNV calls.

    Parameters
    ----------
    fname1 : string
        File name of BED file with CNV calls from experiment 1. Required columns: chr start stop {amplification|deletion}
    fname2 : string
        File name of file with CNV calls from experiment 2 / ground truth. Required columns: chr start stop {amplification|deletion}
    genome : string
        File name of genome file.
    evalgenebased : string
        y or n: perform gene-based evaluation
    frac : float
        Fraction of interval overlap.

    Returns
    -------
    jacc_idx : float
        Jaccard index.
    tp_int : int
        Number of true intervals.
    pp_int : int
        Number of predicted intervals.
    td_int : int
        Number of true discovered intervals.
    gtp_int : int
        Number of ground truth intervals.
    tp_nt : int
        Number of true nucleotides.
    pp_nt : int
        Number of predicted nucleotides.
    gtp_nt : int
        Number of ground truth nucleotides.
    pr_int : float
        Interval precision defined as TP/(TP+FP).
    rc_int : float
        Interval recall defined as TP/(TP+FN).
    pr_nt : float
        Nucleotide precision defined as TP/(TP+FP).
    rc_nt : float
        Nucleotide recall defined as TP/(TP+FN).
    tp_gene : int
        True Positives, gene-based evaluation
    fp_gene : int
        False Positives, gene-based evaluation
    fn_gene : int
        False Negatives, gene-based evaluation
    pr_gene : float
        Gene precision
    rc_gene : float
        Gene recall
    pr_inverted_int : float
        Interval precision without considering matching amp/del status between predicted CNVs and ground truth CNVs
    rc_inverted_int : float
        Interval recall without considering matching amp/del status between predicted CNVs and ground truth CNVs
    pr_inverted_nt : float
        Nucleotide precision without considering matching amp/del status between predicted CNVs and ground truth CNVs
    rc_inverted_nt : float
        Nucleotide recall without considering matching amp/del status between predicted CNVs and ground truth CNVs
    pr_inverted_gene : float
        Gene precision without considering matching amp/del status between predicted CNVs and ground truth CNVs
    rc_inverted_gene : float
        Gene recall without considering matching amp/del status between predicted CNVs and ground truth CNVs
    pr_amp_int : float
        Interval precision - amplifications only
    rc_amp_int : float
        Interval recall - amplifications only
    pr_del_int : float
        Interval precision - deletions only
    rc_del_int : float
        Interval recall - deletions only
    pr_amp_nt : float
        Nucleotide precision - amplifications only
    rc_amp_nt : float
        Nucleotide recall - amplifications only
    pr_del_nt : float
        Nucleotide precision - deletions only
    rc_del_nt : float
        Nucleotide recall - deletions only
    pr_amp_gene : float
        Gene precision - amplifications only
    rc_amp_gene : float
        Gene recall - amplifications only
    pr_del_gene : float
        Gene precision - deletions only
    rc_del_gene : float
        Gene recall - deletions only
    pp_amp: int
        Number of predicted amplifications
    pp_del: int
        Number of predicted deletions
    """

    bed_fname1 = fname1
    bed_fname2 = fname2
    jacc_idx = Jaccard_stats(bed_fname1, bed_fname2, genome)
    tp_int, pp_int, td_int, gtp_int, pr_int, rc_int, pr_inverted_int, rc_inverted_int, pr_amp_int, rc_amp_int, pr_del_int, rc_del_int, pp_amp, pp_del = eval_intervals(bed_fname1, bed_fname2, genome, frac)
    tp_nt, pp_nt, gtp_nt, pr_nt, rc_nt, pr_inverted_nt, rc_inverted_nt, pr_amp_nt, rc_amp_nt, pr_del_nt, rc_del_nt = eval_nt(bed_fname1, bed_fname2, genome)

    if evalgenebased == 'y':
        gene_based_fname1 = bed_fname1.replace('.bed','_gene_based.txt')
        gene_based_fname2 = bed_fname2.replace('.bed','_gene_based.txt')

        if os.path.exists(gene_based_fname1) and os.path.exists(gene_based_fname2):
            print gene_based_fname1, gene_based_fname2
            tp_gene, fp_gene, fn_gene, pr_gene, rc_gene, pr_inverted_gene, rc_inverted_gene, pr_amp_gene, rc_amp_gene, pr_del_gene, rc_del_gene = eval_gene_based(gene_based_fname1, gene_based_fname2, frac)
        else:
            tp_gene, fp_gene, fn_gene, pr_gene, rc_gene, pr_inverted_gene, rc_inverted_gene, pr_amp_gene, rc_amp_gene, pr_del_gene, rc_del_gene = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            print 'WARN: Either '+gene_based_fname1+' and/or '+gene_based_fname2+' does not exist. Not doing gene-based analysis.'
    else:
        tp_gene, fp_gene, fn_gene, pr_gene, rc_gene, pr_inverted_gene, rc_inverted_gene, pr_amp_gene, rc_amp_gene, pr_del_gene, rc_del_gene = float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN')

    return jacc_idx, tp_int, pp_int, td_int, gtp_int, tp_nt, pp_nt, gtp_nt, pr_int, rc_int, pr_nt, rc_nt, tp_gene, fp_gene, fn_gene, pr_gene, rc_gene, pr_inverted_int, rc_inverted_int, pr_inverted_nt, rc_inverted_nt, pr_inverted_gene, rc_inverted_gene, pr_amp_int, rc_amp_int, pr_del_int, rc_del_int, pr_amp_nt, rc_amp_nt, pr_del_nt, rc_del_nt, pr_amp_gene, rc_amp_gene, pr_del_gene, rc_del_gene, pp_amp, pp_del


def main(args):
    """Main function."""

    assert(args.evalgenebased == 'y' or args.evalgenebased == 'n')

    jacc_idx, tp_int, pp_int, td_int, gtp_int, tp_nt, pp_nt, gtp_nt, pr_int, rc_int, pr_nt, rc_nt, tp_gene, fp_gene, fn_gene, pr_gene, rc_gene, pr_inverted_int, rc_inverted_int, pr_inverted_nt, rc_inverted_nt, pr_inverted_gene, rc_inverted_gene, pr_amp_int, rc_amp_int, pr_del_int, rc_del_int, pr_amp_nt, rc_amp_nt, pr_del_nt, rc_del_nt, pr_amp_gene, rc_amp_gene, pr_del_gene, rc_del_gene, pp_amp, pp_del = compare_cnv_calls(args.fname1, args.fname2, args.genome, args.evalgenebased)
    sys.stdout.write('Jaccard index: {0:.3f}\n'.format(jacc_idx))
    sys.stdout.write('Interval level: pr {0:.3f}, rc {1:.3f}\n'.format(pr_int, rc_int))
    sys.stdout.write('Nucleotide level: pr {0:.3f}, rc {1:.3f}\n'.format(pr_nt, rc_nt))
    if args.evalgenebased == 'y':
        sys.stdout.write('Gene level: pr {0:.3f}, rc {1:.3f}\n'.format(pr_gene, rc_gene))

if __name__ == '__main__':
    main(parse_args(sys.argv))
