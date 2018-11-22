#!/usr/bin/env Rscript

####################################################
# Author: Bingxin Lu 
# Description: This script is used to compare the cluster assignments of two different clusterings (simulated truth vs PyClone).
# Input:
#   comparison/pyclone.txt -- The file containing the results of PyClone
#   output/nodes_vars.txt -- The file containing simulated variants in each node
#   output/nodes_ccf.txt  -- The file containing simulated CCF (Cancer Cell Fraction) in each node
# Output:
#   comparison/pycloneVsReal.pdf, which contains the alluvial plot that shows the assignment of mutations in different clusters between the simulated clonal populations and the predicted clonal populations. The y-axis represents the number of mutations. Each rectangle in the white bar at the two sides of the plot represents a cluster. The cluster are ordered by increasing CCF (Cancer Cell Fraction) from top to bottom. The color stripes correspond to mutations in different simulated clusters. CRI: corrected rand index. VI: variation of information. CRI (range from -1 to 1) and VI (non-negative) are two common indexes to measure the agreement between two clusterings. They were computed by method cluster.stats in R library fpc. The larger the CRI is, the more similar the two clustersings are. The smaller the VI is, the more similar the two clustersings are.
##################################################### 

library(ggplot2)
library(fpc)
library(dplyr)
library(tidyr)
library(readr)
library(magrittr)
# install.packages('ggalluvial')
library(ggalluvial)

# The results of PyClone
fpyclone = "comparison/pyclone.txt"
# The simulated variants in each node
fnvar = "output/nodes_vars.txt"
# The simulated CCF (Cancer Cell Fraction) in each node
fnccf = "output/nodes_ccf.txt"

outdir = "comparison"

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

samples = unique(as.character(res_cmb0$sample))
for (i in 1:length(samples)) {
  sampleid = samples[i]
  # print(sampleid) Read SNV infomation for the whole tumor or each sector
  # accordingly
  if (length(samples) > 1) {
    fsnv = paste0("output/phylovar_snvs/sct", i, ".snv")    
  } else {
    fsnv = "output/phylovar_snvs/tumor.snv"
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
  p = ggplot(data = res_all, aes(axis1 = sprintf("%.3f", tumor), axis2 = sprintf("%.3f", 
    ccf))) + scale_x_discrete(limits = c("Simulated", "PyClone"), expand = c(0, 
    0)) + geom_alluvium(aes(fill = as.factor(tumor))) + geom_stratum(width = 1/10) + 
    geom_label(stat = "stratum", label.strata = TRUE) + themes1 + ggtitle(main)
  print(p)
}

dev.off()

