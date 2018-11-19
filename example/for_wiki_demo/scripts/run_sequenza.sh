#!/bin/bash

#####################################################
# Author: Bingxin Lu
# Description: This script is used to run Sequenza
# Dependency: sequenza-utils.py (in Sequenza package, available at $PATH); samtools; scripts/runSequenza.R
# Input: 
#   $1 -- The absolute path to the directory which contains the output of Mutect;
#   $2 -- The absolute path to the directory which contains the output of Sequenza
# Assumption: 
#   The following data files are available at the home directory:
#     input/human_g1k_v37.fasta; 
#     input/S03723314_Covered_c3.bed; 
#     scripts/runSequenza.R; 
#     scripts/Snakefile_sequenza_wes;
#     Two BAM files for a pair of normal and tumor sample under directory $indir(out/normal/normal.bwamem.dedup.realn.recal.bam; out/tumor/tumor.bwamem.dedup.realn.recal.bam)
#####################################################

indir=$1
odir=$2
# indir=$PWD/mutect_output/sct1 
# odir=$PWD/sequenza_output/sct1

home=$PWD

if [ ! -d $odir ]; then
        mkdir -p $odir
fi
cd $odir

# Step 1: create a configuration file called config.yaml, with the folloing content
cat /dev/null > config.yaml
echo -e "exome_bed: $home/input/S03723314_Covered_c3.bed" >> config.yaml
echo -e "reference: $home/input/human_g1k_v37.fasta" >> config.yaml
echo -e "sequenza: $home/scripts/runSequenza.R" >> config.yaml 
echo -e "samples:" >> config.yaml
echo -e " sample:" >> config.yaml

normal=$indir/out/normal/normal.bwamem.dedup.realn.recal.bam
tumor=$indir/out/tumor/tumor.bwamem.dedup.realn.recal.bam
echo -e "   normal: $normal" >> config.yaml
echo -e "   tumor: $tumor" >> config.yaml

# Step 6: Call snakemake to run sequenza with the following command
snakemake --latency-wait 120 --rerun-incomplete -j 100 -d $odir/output --configfile $odir/config.yaml -s $home/scripts/Snakefile_sequenza_wes
# Submiting to a cluster node
# cmd="snakemake --latency-wait 120 --rerun-incomplete -j 100 -d $odir/output --configfile $odir/config.yaml -s $home/scripts/Snakefile_sequenza_wes"
# qsub -V -pe OpenMP 1 -l h_rss=2G,h_rt=480:00:00 -N sequenza -wd $odir/ -o $odir/ -j y -b y "$cmd"
