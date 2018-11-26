#!/usr/bin/env Rscript

####################################################
# Author: Bingxin Lu 
# Description: This script is used to compare the cluster assignments of two different clusterings (simulated truth vs PyClone).
# Input:
#   comparison/pyclone.txt -- The file containing the results of PyClone
#   output/nodes_vars.txt -- The file containing simulated variants in each node
#   output/nodes_ccf.txt  -- The file containing simulated CCF (Cancer Cell Fraction) in each node
#   output/phylovar_snvs/*.tsv  -- The files containing simulated SNVs in a tumor sample
#   input/trunk8000snvs.txt -- The file containing input truncal mutations
# Output:
#   comparison/pyclone_stats.txt, which contains a brief summary of input SNVs to PyClone
#   comparison/pycloneVsReal.pdf, which contains the alluvial plot that shows the assignment of mutations in different clusters between the simulated clonal populations and the predicted clonal populations. The y-axis represents the number of mutations. Each rectangle in the white bar at the two sides of the plot represents a cluster. The cluster are ordered by increasing CCF (Cancer Cell Fraction) from top to bottom. The color stripes correspond to mutations in different simulated clusters. CRI: corrected rand index. VI: variation of information. CRI (range from -1 to 1) and VI (non-negative) are two common indexes to measure the agreement between two clusterings. They were computed by method cluster.stats in R library fpc. The larger the CRI is, the more similar the two clustersings are. The smaller the VI is, the more similar the two clustersings are.
##################################################### 

library(ggplot2)
library(fpc)
library(dplyr)
library(tidyr)
library(readr)
library(magrittr)
library(stringr)
# install.packages('ggalluvial')
library(ggalluvial)

# The results of PyClone
fpyclone = "comparison/pyclone.txt"
# The simulated variants in each node
fnvar = "output/nodes_vars.txt"
# The simulated CCF (Cancer Cell Fraction) in each node
fnccf = "output/nodes_ccf.txt"

indir = "pyclone_output/"
outdir = "comparison"
snvdir = "output/phylovar_snvs"

fstat="comparison/pyclone_stats.txt"
line="A brief summary of input SNVs to PyClone"
write(line, fstat)

# Read input trunk variants
ftrunk="input/trunk8000snvs.txt"
trunks=read_tsv(ftrunk)
names(trunks)=c("chr", "hap", "start", "end", "var")
head(trunks)
trunks  %>% unite(mut, chr, end) -> trunks

# Check the input information of PyClone
# Find files end with .tsv in the input directory of PyClone
ftsvs=list.files(indir, pattern = "*tsv")
pyinput=data.frame()
for(i in 1:length(ftsvs)){
  input=read_tsv(file.path(indir, ftsvs[i]))
  pyinput=rbind(pyinput, input)
}
pyinput %>% tidyr::separate(mutation_id, c("chr", "pos")) %>% dplyr::mutate(chr=str_replace(chr, "chr", "")) %>% unite(mut, chr, pos)-> pyinput

# Read simulated SNVs
fsnvs=list.files(snvdir, pattern = "*snv")
rsnvs=data.frame()
# Summary for each sector
if(length(fsnvs)>1){
  for(i in 1:length(fsnvs)){
    fname = fsnvs[i]
    sample = sub(pattern = "(.*)\\..*$", replacement = "\\1", fname)
    line = paste0("   Sample ", sample, ":")
    write(line, fstat, append = T)

    # Find the input for this sample
    input1=read_tsv(file.path(indir, paste0(sample,".tsv")))
    input1 %>% tidyr::separate(mutation_id, c("chr", "pos")) %>% dplyr::mutate(chr=str_replace(chr, "chr", "")) %>% unite(mut, chr, pos)-> input1
    num_input1=nrow(input1)
    line=paste0("       ", num_input1, " input SNVs to PyClone")
    write(line, fstat, append = T)

    rsnv=read_tsv(file.path(snvdir, fsnvs[i]))
    names(rsnv)=c("chr", "start", "end", "var", "frequency", "rcount", "rfreq")
    rsnv  %>% unite(mut, chr, end) -> rsnv
    rsnvs = rbind(rsnvs, rsnv)

    sim_input = intersect(rsnv$mut, input1$mut)
    num_sinput1 = length(sim_input)
    line = paste0("       ", num_sinput1, " simulated SNVs")
    write(line, fstat, append = T)

    trunk_input = intersect(trunks$mut, sim_input)
    # All the truncal SNVs must be TPs
    truncal_tp = length(trunk_input)
    line = paste0("         ", truncal_tp, " truncal SNVs")
    write(line, fstat, append = T)

    ntruncal_tp = num_sinput1 - truncal_tp
    line= paste0("         ", ntruncal_tp, " non-truncal SNVs")
    write(line, fstat, append = T)

    # Find non-truncal mutations that are not simulated but in the input
    ntruncal_fp = num_input1 - truncal_tp - ntruncal_tp
    line= paste0("       ", ntruncal_fp, " false positive non-truncal SNVs")
    write(line, fstat, append = T)
  }
} else {   # Only one sample
  rsnvs=read_tsv(file.path(snvdir, fsnvs[i]))
  names(rsnvs)=c("chr", "start", "end", "var", "frequency")
  rsnvs  %>% unite(mut, chr, end) -> rsnvs
}

line = paste0("   Total: ")
write(line, fstat, append = T)

num_input = length(unique(pyinput$mut))
line=paste0("       ", num_input, " input SNVs to PyClone")
write(line, fstat, append = T)

sim_input = intersect(unique(rsnvs$mut), unique(pyinput$mut))
num_sinput = length(sim_input)
line = paste0("       ", num_sinput, " simulated SNVs")
write(line, fstat, append = T)

trunk_input = intersect(trunks$mut, sim_input)
# All the truncal SNVs must be TPs
truncal_tp = length(trunk_input)
line = paste0("         ", truncal_tp, " truncal SNVs")
write(line, fstat, append = T)

ntruncal_tp = num_sinput - truncal_tp
line= paste0("         ", ntruncal_tp, " simulated non-truncal SNVs")
write(line, fstat, append = T)

# Find non-truncal mutations that are not simulated but in the input
ntruncal_fp = num_input - truncal_tp - ntruncal_tp
line= paste0("       ", ntruncal_fp, " false positive non-truncal SNVs")
write(line, fstat, append = T)


nodes_vars = read_tsv(fnvar)
names(nodes_vars) = c("node", "chr", "start", "end", "var")
nodes_ccf = read_tsv(fnccf)
colnames(nodes_ccf)[1] = c("node")
nodes_vars_ccf = merge(nodes_vars, nodes_ccf, by = c("node"))
nodes_vars_ccf <- nodes_vars_ccf %>% select(-start) %>% mutate(node = gsub("node", 
  "", node)) %>% unite(mut, chr, end)

# columns: 1:chr 2:startpos 3:sample 4:cluster_id 5:cellular_prevalence
# 6:cellular_prevalence_std 7:variant_allele_frequency 8:mean Column 8 are from
# cluster.tsv
res_pyclone = read.table(fpyclone, header = F, sep = "\t")
names(res_pyclone) = c("chr", "pos", "sample", "cluster", "prevalence", "prevalence_std", 
  "vaf", "ccf")
res_pyclone <- res_pyclone %>% mutate(chr = gsub("chr", "", chr)) %>% unite(mut, 
  chr, pos)

res_cmb0 = merge(nodes_vars_ccf, res_pyclone, by = c("mut"))

# Write figures to files
fname = "pycloneVsReal.pdf"
outfile = file.path(outdir, fname)
pdf(outfile, width = 11, height = 11 * 0.618)

if(length(unique(res_cmb0$sample))==1) {res_cmb0$sample="tumor"}
samples = unique(as.character(res_cmb0$sample))
for (i in 1:length(samples)) {
  sampleid = samples[i]
  # print(sampleid) Read SNV infomation for the whole tumor or each sector
  # accordingly
  if (length(samples) > 1) {
    fsnv = file.path(snvdir, paste0("sct", i, ".snv"))    
  } else {
    fsnv = file.path(snvdir, "tumor.snv")
  }
  snvs = read_tsv(fsnv)
  if(ncol(snvs)==5){
    names(snvs) = c("chr", "start", "end", "var", "frequency")
  } else {
    names(snvs) = c("chr", "start", "end", "var", "frequency", "rcount", "rfreq")
  }  
  snvs <- snvs %>% select(-start, -var) %>% unite(mut, chr, end)
  
  res_cmb = merge(res_cmb0, snvs, by = c("mut"))
  res_all = res_cmb %>% dplyr::filter(sample == sampleid)
  
  res_all$Real_clustering = as.factor(res_all$node)
  res_all$PyClone_clustering = as.factor(res_all$cluster)
  levels(res_all$Real_clustering) = seq(1, length(levels(res_all$Real_clustering)))
  levels(res_all$PyClone_clustering) = seq(1, length(levels(res_all$PyClone_clustering)))
  res_all$Real_clustering = as.integer(res_all$Real_clustering)
  res_all$PyClone_clustering = as.integer(res_all$PyClone_clustering)
  stat = cluster.stats(d = NULL, clustering = res_all$PyClone_clustering, alt.clustering = res_all$Real_clustering, 
    compareonly = TRUE)
  cri = sprintf("%.3f", stat$corrected.rand)
  vi = sprintf("%.3f", stat$vi)
  
  # Plot the figure
  fontsize = 18
  axissize = 18
  titlesize = 22
  themes1 <- theme_bw() + theme(legend.text = element_text(size = fontsize), legend.title = element_text(size = fontsize), 
    axis.title = element_text(size = axissize), axis.text = element_text(size = axissize), 
    plot.title = element_text(size = titlesize, hjust = 0.5), legend.position = "none", 
    axis.line = element_line(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  
  if (length(samples) > 1) {
    prefix = paste("Mutation assignment", i)
  } else {
    prefix = "Mutation assignment"
  }
  
  main = paste0(prefix, " (CRI: ", cri, "; VI: ", vi, ")")
  p = ggplot(data = res_all, aes(axis1 = sprintf("%.3f", res_all[,sampleid]), axis2 = sprintf("%.3f", 
    ccf))) + scale_x_discrete(limits = c("Simulated", "PyClone"), expand = c(0, 
    0)) + geom_alluvium(aes(fill = as.factor(res_all[,sampleid]))) + geom_stratum(width = 1/10) + 
    geom_label(stat = "stratum", label.strata = TRUE) + themes1 + ggtitle(main)
  print(p)
}

dev.off()

