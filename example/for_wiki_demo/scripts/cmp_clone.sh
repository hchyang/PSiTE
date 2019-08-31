#!/bin/bash

#####################################################
# Author: Bingxin Lu
# Description: This script is used to compare the output of PyClone (predicted clonal compositions) against simulated clonal compositions
# Input: $1 -- Y (The sample is single tumor); N (The sample is multi-regional tumor)
# Output: 
#   pyclone.txt, which contains the predictions of PyClone;
#   pyclone_stats.txt, which contains a brief summary of PyClone precitions; 

# Assumption: 
#   The following files are available:
#     pyclone_output/tables/cluster.tsv; 
#     pyclone_output/tables/loci.tsv; 
#     output/nodes_vars.txt;
#     output/nodes_ccf.txt;
#     output/phylovar_snvs/*.snv; 
#     scripts/compare_clone.R
# Command: bash cmp_clone.sh
#####################################################


# Postprocess the output of PyClone
idir=pyclone_output/tables
wdir=comparison

# Combine the original two output tables into a single file
cluster=$idir/cluster.tsv
loci=$idir/loci.tsv
# columns in $patient.pyclone: 1:chr  2:startpos  3:sample_id   4:cluster_id   5:cellular_prevalence   6:cellular_prevalence_std   7:variant_allele_frequency 8:mean
# Column 8 are from cluster.tsv
output=$wdir/pyclone.txt
awk 'NR>1 && $3>1' $cluster | perl -e 'open CLUSTER,$ARGV[0];while(<CLUSTER>){@line=split;$key=$line[0].$line[1]; $h{$key}=$line[3];} open LOCI,$ARGV[1];while(<LOCI>){@line=split;$ key=$line[1].$line[2]; if($h{$key}){$ccf=$h{$key}; print join("\t",@line, $ccf),"\n"; } }' - $loci | perl -lane '$F[0]=~s/.*_//;$F[0]=~s/:/\t/;$F[1]=~s/.*_//;print join("\t",@F)' | sort -k1,1 -k2,2n > $output


# Draw plot of mutation assignment
# Only show mutations in both output
# Columns required for the plots: CCFs for each muts
# Assign CCF for each node
Rscript scripts/compare_clone.R

