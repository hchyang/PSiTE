#!/usr/bin/env python

####################################################
# Author: Bingxin Lu 
# Description: This script is used to convert a NHX tree to a NEWICK tree.
# Input:
#   input/tumopp_subtree.nhx -- The file containing a phylogenetic tree in NHX format
# Output:
#   input/tumopp_subtree.nwk, which contains a phylogenetic tree in NEWICK format
# Command:  python nhx2nwk.py
##################################################### 

# http://etetoolkit.org/docs/latest/tutorial/tutorial_trees.html#reading-newick-trees

from ete3 import Tree

ifile="input/tumopp_subtree.nhx"
ofile="input/tumopp_subtree.nwk"

t = Tree(ifile)
t.write(format=5, outfile=ofile)
