#!/usr/bin/env python3

#########################################################################
# Author: Hechuan Yang
# Created Time: 2017-01-13 17:06:00
# File Name: highlight_mutations.py
# Description: 
#########################################################################

import os
import sys
import argparse
import pickle
import logging

from simulate_somatic_vars import *

#import parse_newick

parse=argparse.ArgumentParser(description='Find the mutations on the tree and output a NHX tree with the branch with the muations highlighted')
parse.add_argument('-m','--mutations',required=True,help='Mutations to be selected')
parse.add_argument('-l','--log',required=True,help='Log file')
parse.add_argument('-t','--tree', help='The original tree with mutations')
args=parse.parse_args()
logging.basicConfig(filename=args.log, filemode='w', format='%(levelname)s: %(message)s', level=10)

mutations={}
with open(args.mutations) as input:
    for line in input:
        line=line.rstrip()
        fields=line.split('\t')
        mutations[str(fields[0])]=float(fields[1])

with open(args.tree,'rb') as input:
    mytree=pickle.load(input)
    mytree.highlight_snvs(mutations.values())
#    print(set(mutations.values()))
    print('{};'.format(mytree.tree2newick(with_lens=True,attrs=['C'])))

