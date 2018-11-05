#!/bin/bash

#####################################################
# Author: Bingxin Lu
# Description: This script is used to check the output of phylovar (simulated SNVs and CNVs). Please run this script to check the number of simulated SNVs and CNVs. 
# Input: None
# Output: 
#   phylovar_stats.txt, which contains simple summary of simulated SNVs and CNVs; 
#   tumor_snv_exome.txt, a TSV file containing simulated SNVs in the exome region; 
#   tumor_cnv_exome.txt, a TSV file containing simulated CNVs in the exome region
# Assumption: The following files are available -- input/S03723314_Covered_c3.bed; input/trunk8000snvs.txt; output/phylovar_snvs/tumor.snv; output/phylovar_cnvs/tumor.cnv
# Command: bash check_phylovar.sh
#####################################################

ftarget=input/S03723314_Covered_c3.bed
ftrunk=input/trunk8000snvs.txt
fsnv=output/phylovar_snvs/tumor.snv
fcnv=output/phylovar_cnvs/tumor.cnv

wdir="comparison"   # The directory to store the output of comparison
if [ ! -d $wdir ]; then
  mkdir $wdir 
fi 

fsnv_exome=$wdir/tumor_snv_exome.txt
fcnv_exome=$wdir/tumor_cnv_exome.txt

fcnv_bed=$wdir/tumor_cnv.bed
fcnv_exome_bed=$wdir/tumor_cnv_exome.bed

# Parse the simulated SNVs
# Find regions in $fsnv that are also in $ftarget
bedtools intersect -u -a $fsnv -b $ftarget > $fsnv_exome 
# Extract only the position
cut -f1,3 $fsnv_exome| sed 's/\t/_/g' > $wdir/tumor_snv_exome.pos
cut -f1,3,4,5 $fsnv |> $wdir/real_snv.pos
# trunk_snv.pos, which contains simulated truncal SNVs;
tail -n+2 $ftrunk | awk '$5 !~ /^+.*/ && $5 !~ /^-.*/' | cut -f1,4 | sed 's/\t/_/g' > $wdir/trunk_snv.pos


# Parse the simulated CNVs
# trunk_cnv.pos, which contains simulated truncal CNVs.   
tail -n+2 $ftrunk | awk '$5 ~ /^+.*/ || $5 ~ /^-.*/' | cut -f1,3,4,5  > $wdir/trunk_cnv.pos
awk 'BEGIN{OFS="\t"}NR>1{if($4 ~ /-/) print $1,$2,$3,"deletion"; else print $1,$2,$3,"amplification"}'  $fcnv | sed 's/^X/999/g' - | sort -k1n -k2n | sed 's/^999/X/g' - > $fcnv_bed
# Extract CNVs in exome region
# $ftarget contains a set of small regions; $fcnv contains a set of large regions; -F Minimum overlap required as a fraction of B. -u:  Write the original A entry once if any overlaps found in B.
bedtools intersect -u -a $fcnv -b $ftarget -F 0.5 > $fcnv_exome
awk 'BEGIN{OFS="\t"}{if($4 ~ /-/) print $1,$2,$3,"deletion"; else print $1,$2,$3,"amplification"}'  $fcnv_exome | sed 's/^X/999/g' - | sort -k1n -k2n | sed 's/^999/X/g' - > $fcnv_exome_bed
bedtools intersect -u -a $wdir/trunk_cnv.pos -b $ftarget -F 0.5 > $wdir/trunk_cnv_exome

# Write evaluation statistic to a file
fstat="$wdir/phylovar_stats.txt"
cat /dev/null > $fstat

## SNVs
total=`tail -n+2 $fsnv | wc -l`
echo "The number of SNVs simulated: $total" >> $fstat
truncal=`tail -n+2 $ftrunk | awk '$5 !~ /^+.*/ && $5 !~ /^-.*/' | wc -l`
echo "  $truncal truncal SNVs" >> $fstat
ntruncal=`echo $total - $truncal | bc`
echo "  $ntruncal non-truncal SNVs" >> $fstat
echo -e "\n" >> $fstat

real=`cat $wdir/tumor_snv_exome.pos | wc -l`
echo "The number of SNVs in exome region: $real" >> $fstat
truncal=`grep -f $wdir/tumor_snv_exome.pos $wdir/trunk_snv.pos | wc -l`
echo "  $truncal truncal SNVs" >> $fstat
ntruncal=`echo $real - $truncal | bc`
echo "  $ntruncal non-truncal SNVs" >> $fstat
echo -e "\n" >> $fstat


## CNVs
total=`tail -n+2 $fcnv | wc -l`
ndel=`grep del $fcnv_bed | wc -l`
namp=`grep amp $fcnv_bed | wc -l`
echo "The number of CNVs simulated: $total" >> $fstat
echo "  $ndel deletions" >> $fstat
echo -e "  $namp amplifications\n" >> $fstat
truncal=`tail -n+2 $ftrunk | awk '$5 ~ /^+.*/ || $5 ~ /^-.*/' | wc -l`

echo "The number of truncal CNVs simulated: $truncal" >> $fstat
ndelt=`awk '$5 ~ /^-.*/' $ftrunk | wc -l`
nampt=`awk '$5 ~ /^+.*/' $ftrunk | wc -l`
echo "  $ndelt deletions" >> $fstat
echo -e "  $nampt amplifications\n" >> $fstat

ntruncal=`echo $total - $truncal | bc`
echo "The number of non-truncal CNVs simulated: $ntruncal" >> $fstat
ndeln=`echo $ndel - $ndelt | bc`
nampn=`echo $namp - $nampt | bc`
echo "  $ndeln deletions" >> $fstat
echo -e "  $nampn amplifications\n" >> $fstat

total=`cat $fcnv_exome_bed | wc -l`
ndel=`grep del $fcnv_exome_bed | wc -l`
namp=`grep amp $fcnv_exome_bed | wc -l`
echo "The number of CNVs overlaping with target regions: $total" >> $fstat
echo "  $ndel deletions" >> $fstat
echo -e "  $namp amplifications\n" >> $fstat

truncal=`cat $wdir/trunk_cnv_exome | wc -l`
echo "The number of truncal CNVs overlaping with target regions: $truncal" >> $fstat
ndelt=`awk '$4 ~ /^-.*/' $wdir/trunk_cnv_exome | wc -l`
nampt=`awk '$4 ~ /^+.*/' $wdir/trunk_cnv_exome | wc -l`
echo "  $ndelt deletions" >> $fstat
echo -e "  $nampt amplifications\n" >> $fstat

ntruncal=`echo $total - $truncal | bc`
echo "The number of non-truncal CNVs overlaping with target regions: $ntruncal" >> $fstat
ndeln=`echo $ndel - $ndelt | bc`
nampn=`echo $namp - $nampt | bc`
echo "  $ndeln deletions" >> $fstat
echo -e "  $nampn amplifications\n" >> $fstat

# Remove intermediate files
rm $wdir/tumor_snv_exome.pos
rm $wdir/real_snv.pos
rm $wdir/trunk_snv.pos

rm $fcnv_exome_bed
rm $fcnv_bed
rm $wdir/trunk_cnv.pos
rm $wdir/trunk_cnv_exome
