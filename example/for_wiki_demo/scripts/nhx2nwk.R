#!/usr/bin/env Rscript

####################################################
# Author: Bingxin Lu 
# Description: This script is used to convert a NHX tree to a NEWICK tree.
# Input:
#   input/tumopp_subtree.nhx -- The file containing a phylogenetic tree in NHX format
# Output:
#   input/tumopp_subtree.nwk, which contains a phylogenetic tree in NEWICK format
# Command:  Rscript nhx2nwk.R
##################################################### 

library(ggtree)
library(readr)
library(ape)

idir="input"
odir="input"

itree = read.nhx(file.path(idir, "tumopp_subtree.nhx"))
# Convert nhx to newick format, input to phylovar
ftree=file.path(odir,"tumopp_subtree.nwk")
write.tree(itree@phylo, ftree)

