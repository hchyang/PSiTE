#!/usr/bin/env Rscript

#####################################################
# Author: Bingxin Lu
# Description: This script is used towrite an affiliation file which records the relationship between cells and samples, based on the output file of tumopp (sampled_ids.txt).
#####################################################


library(readr)

idir="input"
odir="input"


# Write the affiliation file
fname="affiliation_tumopp_random.txt"
fpath=file.path(odir, fname)
fout=file(fpath, 'w')
header="#sector purity depth prune_p cells"
purity=0.6
depth=100
prune_p=0.05
write_lines(header, fout)

# Read sector information
fpath=file.path(odir, "sampled_ids.txt")
sampled_ids=read_tsv(fpath, col_types = "cic")
names(sampled_ids)=c("sector", "count", "id")
# Count the number of cells for each sample
sample_size = unlist(lapply(sampled_ids$id, function(x) length(strsplit(x, ",")[[1]])))
print(sample_size)
for(i in 1:nrow(sampled_ids)){
  # i = 5
  sect_info=sampled_ids[i,]
  sect = sect_info$sector
  cells = sect_info$id
  line = paste(sect, purity, depth, prune_p, cells, sep=" ")
  write_lines(line, fout)
}
close(fout)

