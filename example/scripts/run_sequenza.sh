#!/bin/bash

#####################################################
# Author: Bingxin Lu
# Description: This script is used to run Sequenza
# Dependency: sequenza-utils.py (in Sequenza package, available at $PATH); samtools; scripts/runSequenza.R
# Input: 
#   input/human_g1k_v37_decoy.fasta
#   input/S03723314_Covered_c3.bed
#####################################################

home=$PWD

fname=sequenza_output
if [ ! -d $fname ]; then
        mkdir $fname
fi
cd $fname
wdir=$PWD

# Step 1: create a configuration file called config.yaml, with the folloing content
cat /dev/null > config.yaml
echo -e "exome_bed: $home/input/S03723314_Covered_c3.bed" >> config.yaml
echo -e "reference: $home/input/human_g1k_v37_decoy.fasta" >> config.yaml
echo -e "sequenza: $home/scripts/runSequenza.R" >> config.yaml 
echo -e "samples:" >> config.yaml
echo -e " sample:" >> config.yaml
indir=$home/mutect_output/
normal=`ls $indir/out/normal/normal.bwamem.dedup.realn.recal.bam`
tumor=`ls $indir/out/tumor/tumor.bwamem.dedup.realn.recal.bam`
echo -e "   normal: $normal" >> config.yaml
echo -e "   tumor: $tumor" >> config.yaml

# Step 6: Call snakemake to run sequenza with the following command
snakemake -n --latency-wait 120 --rerun-incomplete -j 100 -d $wdir/output --configfile $wdir/config.yaml -s $home/scripts/Snakefile_sequenza_wes
