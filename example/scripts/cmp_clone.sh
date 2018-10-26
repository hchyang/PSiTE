#!/bin/bash

#####################################################
# Author: Bingxin Lu
# Description: This script is used to compare the output of PyClone (predicted clonal compositions) against simulated clonal compositions
# Input: None
# Output: 
#   pyclone.txt, which contains the predictions of PyClone;
#   pyclone_stats.txt, which contains a brief summary of PyClone precitions; 
#   mutect_snv.txt, a TSV file containing SNVs predicted by Mutect; 
#   real_snv.txt, a one-column file containing the position of simulated SNVs; 
#   trunk_snv.pos, a one-column file containing the position of truncal SNVs
# Assumption: 
#   The following files are available -- pyclone_output/tables/cluster.tsv; pyclone_output/tables/loci.tsv; output/phylovar_snvs/tumor.snv; comparison/tumor_snv_exome.txt;  comparison/trunk_snv.pos
#   The script file scripts/compare_clone.R
# Command: bash cmp_clone.sh
#####################################################


# Postprocess the output of PyClone
idir=pyclone_output/tables
wdir=comparison

cluster=$idir/cluster.tsv
loci=$idir/loci.tsv
# columns in $patient.pyclone: 1:chr  2:startpos  3:sample_id   4:cluster_id   5:cellular_prevalence   6:cellular_prevalence_std   7:variant_allele_frequency 8:mean
# Column 8 are from cluster.tsv
output=$wdir/pyclone.txt
awk 'NR>1 && $3>1' $cluster | perl -e 'open CLUSTER,$ARGV[0];while(<CLUSTER>){@line=split;$key=$line[0].$line[1]; $h{$key}=$line[3];} open LOCI,$ARGV[1];while(<LOCI>){@line=split;$ key=$line[1].$line[2]; if($h{$key}){$ccf=$h{$key}; print join("\t",@line, $ccf),"\n"; } }' - $loci | perl -lane '$F[0]=~s/.*_//;$F[0]=~s/:/\t/;$F[1]=~s/.*_//;print join("\t",@F)' | sort -k1,1 -k2,2n > $output


# Find the statistics of input SNVs
# Write evaluation statistic to a file
fstat="$wdir/pyclone_stats.txt"
cat /dev/null > $fstat

fpos=$wdir/pyclone_input.pos
fsnv="output/phylovar_snvs/tumor.snv"
fsnv_exome="$wdir/tumor_snv_exome.txt"

less pyclone_output/input.tsv | tail -n+2 | cut -f1 | sed 's/chr//g' | sed 's/:/_/g' > $fpos
total=`less $fpos | wc -l`
echo "The number of input SNVs to PyClone: $total" >> $fstat

# All the truncal SNVs must be TPs
truncal=`grep -f $fpos $wdir/trunk_snv.pos | wc -l`
echo "  $truncal truncal SNVs" >> $fstat
ntruncal=`echo $total - $truncal | bc`

# The other mutations may be FPs
# Input mutations to PyClone that are not truncal
less $wdir/real_snv.txt | tail -n+2 | cut -f1,2 | sed 's/\t/_/g' > $wdir/real_snv.pos
ntruncal_tp=`grep -f $fpos $wdir/trunk_snv.pos  -v | grep -f - $wdir/real_snv.pos | wc -l`
echo "  $ntruncal_tp simulated non-truncal SNVs" >> $fstat
ntruncal_fp=`echo $ntruncal - $ntruncal_tp | bc`
echo "  $ntruncal_fp FP non-truncal SNVs" >> $fstat
echo -e "\n" >> $fstat


# Draw plot of mutation assignment
# Only show mutations in both output
# Columns required for the plots: CCFs for each muts
# Assign CCF for each node
Rscript scripts/compare_clone.R

# Remove intermediate files
rm $wdir/pyclone_input.pos
rm $wdir/real_snv.pos