#!/bin/bash

#####################################################
# Author: Bingxin Lu
# Description: This script is used to compare the output of Mutect (predicted SNVs) against simulated SNVs
# Input: $1 -- The folder storing mutect.PASS.vcf.gz; $2 -- The name of output file
# Output: 
#   The file specified by $2, which contains a brief summary of SNV precitions; 
#   tumor_snv_exome.txt, a TSV file containing simulated SNVs in the exome region; 
#   mutect_snv.txt, a TSV file containing SNVs predicted by Mutect; 
#   real_snv.txt, a TSV file containing the position of simulated SNVs; 
#   trunk_snv.pos, a one-column file containing the position of truncal SNVs
# Assumption: The following files are available -- input/S03723314_Covered_c3.bed; input/trunk8000snvs.txt; output/phylovar_snvs/tumor.snv; $1/mutect.PASS.vcf.gz
# Command: bash check_snv.sh $1 $2
#####################################################

dir_mutect=$1
fout=$2
# fout="snv_stats.txt"

ftarget="input/S03723314_Covered_c3.bed"
ftrunk="input/trunk8000snvs.txt"
fsnv="output/phylovar_snvs/tumor.snv"
res_mutect="$dir_mutect/mutect.PASS.vcf.gz"
# threshold for low-frequency SNVs
cutoff=0.1

wdir="comparison"   # The directory to store the output of comparison
if [ ! -d $wdir ]; then
  mkdir $wdir 
fi 

# Parse the simulated SNVs
fsnv_exome="$wdir/tumor_snv_exome.txt"
# Find regions in $fsnv that are also in $ftarget
bedtools intersect -u -a $fsnv -b $ftarget > $fsnv_exome
# Extract only the position
cut -f1,3 $fsnv_exome | sed 's/\t/_/g' > $wdir/tumor_snv_exome.pos
cut -f1,3,4,5 $fsnv > $wdir/real_snv.txt
tail -n+2 $ftrunk | awk '$5 !~ /^+.*/ && $5 !~ /^-.*/' | cut -f1,4 | sed 's/\t/_/g' > $wdir/trunk_snv.pos

# Parse the predicted SNVs
zcat $res_mutect | awk '!/^#/' | cut -f1,2 > $wdir/mutect_snv.pos
sed -i 's/\t/_/g' $wdir/mutect_snv.pos
zcat $res_mutect | awk '!/^#/' | cut -f1,2,4,5,10,11 | perl -lane '@s1=split/:/,$F[4]; @s2=split/:/,$F[5]; print join("\t", $F[0], $F[1], $F[2], $F[3], $s1[4], $s2[4]);' -  > $wdir/mutect_snv.txt

# Find FNs and FPs
fpos=$wdir/mutect_snv.pos
fwes=$wdir/tumor_snv_exome.pos
ffn=$wdir/snv.mutect.fn
ffp=$wdir/snv_wes.mutect.fp
diff <(sort $fpos) <(sort $fwes) | grep '^>' | sed 's/^>\ //' > $ffn
diff <(sort $fwes) <(sort $fpos) | grep '^>' | sed 's/^>\ //' > $ffp

# Write evaluation statistic to a file
fstat="$wdir/$fout"
cat /dev/null > $fstat

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

# Find truncal and non-truncal mutations in the predicted SNVs
pred=`cat $fpos | wc -l`
echo "Predicted SNVs by mutect: $pred" >> $fstat
truncal=`grep -f $fpos $wdir/trunk_snv.pos | wc -l`
echo "  $truncal truncal SNVs" >> $fstat
ntruncal=`echo $pred - $truncal | bc`
echo "  $ntruncal non-truncal SNVs" >> $fstat
echo -e "\n" >> $fstat


fn=`cat $ffn | wc -l`
echo "FNs in the predictions: $fn" >> $fstat
fp=`cat $ffp | wc -l`
echo "FPs in the predictions: $fp" >> $fstat
tp=`echo "$pred - $fp" | bc`
echo "TPs in the predictions: $tp" >> $fstat
pr=`echo "scale=2; $tp / $pred" | bc`
echo "Precison: $pr" >> $fstat
rc=`echo "scale=2; $tp / $real" | bc`
echo "Recall: $rc" >> $fstat
echo -e "\n" >> $fstat

# Find the frequency of FPs and FNs
# Get simulated (real) frequency of FNs
sed 's/_/\t/g' $ffn | grep -F -f - $wdir/real_snv.txt | sort -k4n > $ffn.freq
# Find FNs with low frequency
nlf=`awk -v cutoff=$cutoff '$4<cutoff' $ffn.freq | wc -l`
echo "The number of FNs with frequency lower than $cutoff: $nlf"  >> $fstat
total=`cat $ffn  | wc -l`
# echo "The number of total FNs: $total" >> $fstat
frac=`echo "scale=2; $nlf / $total" | bc`
echo "The fraction of FNs with frequency lower than $cutoff: $frac" >> $fstat
echo -e "\n" >> $fstat

# Get predicted frequency of FNs
sed 's/_/\t/g' $ffp > tmp
zcat $res_mutect | grep -F -f tmp - | cut -f1,2,4,5,10,11 | perl -lane '@s1=split/:/,$F[4]; @s2=split/:/,$F[5]; print join("\t", $F[0], $F[1], $F[2], $F[3], $s1[4], $s2[4]);' -  > $ffp.freq
# Find FPs with low frequency
nlf=`awk -v cutoff=$cutoff '$6>cutoff' $ffp.freq | wc -l`
echo "The number of FPs with frequency lower than $cutoff: $nlf" >> $fstat
total=`cat $ffp.freq | wc -l`
# echo "The number of total FNs: $total" >> $fstat
frac=`echo "scale=2; $nlf / $total" | bc`
echo "The fraction of FPs with frequency lower than $cutoff: $frac" >> $fstat

# Remove intermediate files
rm tmp
rm $fpos
rm $fwes
rm $ffn
rm $ffp
rm $ffn.freq
rm $ffp.freq
