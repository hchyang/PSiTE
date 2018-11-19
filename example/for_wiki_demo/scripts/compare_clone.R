#!/usr/bin/env Rscript

#####################################################
# Author: Bingxin Lu
# Description: This script is used to compare the cluster assignments of two different clusterings (simulated truth vs PyClone). 
# The alluvial plot shows the assignment of mutations in different clusters between the simulated clonal populations and the predicted clonal populations. The y-axis represents the number of mutations. Each rectangle in the white bar at the two sides of the plot represents a cluster. The cluster are ordered by increasing CCF (Cancer Cell Fraction) from top to bottom. On the left, there are 9 clusters in the simulated data. On the right, there are 2 clusters predicted by PyClone. The color stripes correspond to mutations in different simulated clusters. CRI: corrected rand index. VI: variation of information. CRI (range from -1 to 1) and VI (non-negative) are two common indexes to measure the agreement between two clusterings. They were computed by method cluster.stats in R library fpc. The larger the CRI is, the more similar the two clustersings are. The smaller the VI is, the more similar the two clustersings are.
#####################################################

library(ggplot2)
library(fpc)
library(dplyr)
library(tidyr)
library(readr)
library(magrittr)
# install.packages("ggalluvial")
library(ggalluvial)

fpyclone="comparison/pyclone.txt"
fnvar="output/nodes_vars.txt"
fnccf="output/nodes_ccf.txt"
outdir="comparison"

nodes_vars=read_tsv(fnvar)
names(nodes_vars)=c("node", "chr", "start", "end", "var")
nodes_ccf=read_tsv(fnccf)
colnames(nodes_ccf)[1]=c("node")
nodes_vars_ccf=merge(nodes_vars, nodes_ccf, by=c("node"))
nodes_vars_ccf %>% select(-start) %>% mutate(node=gsub("node","",node)) %>% unite(mut, chr, end)  -> nodes_vars_ccf

# columns: 1:chr   2:startpos   3:sample     4:cluster_id   5:cellular_prevalence   6:cellular_prevalence_std       7:variant_allele_frequency 8:mean
# Column 8 are from cluster.tsv
res_pyclone = read.table(fpyclone, header=F, sep="\t")
names(res_pyclone)=c("chr", "pos", "sample", "cluster", "prevalence", "prevalence_std", "vaf", "ccf")
res_pyclone %>% mutate(chr=gsub("chr","",chr)) %>% unite(mut, chr, pos) -> res_pyclone

res_cmb = merge(nodes_vars_ccf, res_pyclone, by=c("mut"))
# Write figures to files
fname="pycloneVsReal.pdf"
outfile=file.path(outdir, fname)
pdf(outfile, width=11, height=11*0.618)

samples = unique(as.character(res_cmb$sample))
for (i in 1:length(samples)){
  sampleid=samples[i]
  # print(sampleid)
  res_all = res_cmb %>% dplyr::filter(sample==sampleid)
  
  res_all$Real_clustering=as.factor(res_all$node)
  res_all$PyClone_clustering=as.factor(res_all$cluster)
  levels(res_all$Real_clustering)=seq(1,length(levels(res_all$Real_clustering)))
  levels(res_all$PyClone_clustering)=seq(1,length(levels(res_all$PyClone_clustering)))
  res_all$Real_clustering=as.integer(res_all$Real_clustering)
  res_all$PyClone_clustering=as.integer(res_all$PyClone_clustering)
  stat=cluster.stats(d = NULL, clustering=res_all$PyClone_clustering, alt.clustering = res_all$Real_clustering, compareonly = TRUE)
  cri=sprintf("%.3f", stat$corrected.rand)
  vi=sprintf("%.3f", stat$vi)
  
  # Plot the figure
  fontsize=18
  axissize=18
  titlesize=22
  themes1 <- theme_bw() + theme(legend.text = element_text(size=fontsize)
                                , legend.title = element_text(size=fontsize)
                                , axis.title=element_text(size=axissize), axis.text=element_text(size=axissize)
                                , plot.title = element_text(size=titlesize, hjust = 0.5)
                                , legend.position="none"
                                , axis.line=element_line()
                                , panel.grid.major = element_blank(), panel.grid.minor = element_blank()
  )
  
  
  prefix="mutation assignment"
  main=paste0(prefix, " (CRI: ", cri, "; VI: ", vi,")")
  p=ggplot(data = res_all,
           aes(axis1 = sprintf("%.3f",tumor), axis2 = sprintf("%.3f",ccf))) +
    scale_x_discrete(limits = c("Real", "PyClone"), expand=c(0, 0)) +
    geom_alluvium(aes(fill = as.factor(tumor))) +
    geom_stratum(width = 1/10) + geom_label(stat = "stratum", label.strata = TRUE) +
    themes1 + ggtitle(main)
  print(p)
}

dev.off()
