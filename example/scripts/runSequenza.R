#!/usr/bin/env Rscript
library("sequenza")
args <- commandArgs(trailingOnly = T)
data.file <- args[1]

seqz <- sequenza.extract(data.file,gz=TRUE,chromosome.list=c(1:22,'X'))
cp.data <- sequenza.fit(seqz)

base_name=gsub('.binning.seqz.gz$','',basename(data.file))
output_dir=dirname(data.file)

sequenza.results(sequenza.extract=seqz,cp.table=cp.data,sample.id="output",out.dir=paste(output_dir,base_name,sep='/'))
