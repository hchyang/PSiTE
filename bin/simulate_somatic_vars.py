#!/usr/bin/env python

#########################################################################
# Author: Hechuan
# Created Time: 2017-04-04 18:00:34
# File Name: simulate_somatic_vars.py
# Description: 
#########################################################################

import sys
import pickle
import numpy
import argparse
import logging
import trunk_vars
import tree
#import pdb; pdb.set_trace()

#handle the error below
#python | head == IOError: [Errno 32] Broken pipe 
from signal import signal, SIGPIPE, SIG_DFL 
signal(SIGPIPE,SIG_DFL) 

#@profile
def __main__():
    parse=argparse.ArgumentParser(description='Generate snvs on a coalescent tree in newick format')
    parse.add_argument('-t','--tree',required=True,help='a tree in newick format')
    default=300
    parse.add_argument('-r','--snv_rate',type=float,default=default,help='the muation rate of SNVs [{}]'.format(default))
    default=3
    parse.add_argument('-R','--cnv_rate',type=float,default=default,help='the muation rate of CNVs [{}]'.format(default))
    default=0.5
    parse.add_argument('-d','--del_prob',type=int,default=default,help='the probability of being deletion for a CNV mutation [{}]'.format(default))
    default=4000000
    parse.add_argument('-l','--cnv_length_lambda',type=int,default=default,help='the lambda of CNVs length [{}]'.format(default))
    default=20000000
    parse.add_argument('-L','--cnv_length_max',type=int,default=default,help='the maximium of CNVs length [{}]'.format(default))
    default=5
    parse.add_argument('-c','--copy_max',type=int,default=default,help='the maximium ADDITIONAL copy of a CNV [{}]'.format(default))
    default=2
    parse.add_argument('-p','--ploid',type=int,default=default,help='the ploid to simulate [{}]'.format(default))
    default=0
    parse.add_argument('-x','--prune',type=int,default=default,help='trim all its children for the branches with equal or less than this number of tips [{}]'.format(default))
    default=0.0
    parse.add_argument('-X','--prune_proportion',type=float,default=default,help='trim all its children for the branches with equal or less than this proportion of tips [{}]'.format(default))
    default=50
    parse.add_argument('-D','--depth',type=int,default=default,help='the mean depth for simulating coverage data [{}]'.format(default))
    default=None
    parse.add_argument('-s','--random_seed',type=int,help='the seed for random random number generator [{}]'.format(default))
    default='raw.cnvs'
    parse.add_argument('-V','--cnv',type=str,default=default,help='the output file to save CNVs [{}]'.format(default))
    default='raw.snvs'
    parse.add_argument('-S','--snv',type=str,default=default,help='the output file to save SNVs [{}]'.format(default))
    default='depth.profile'
    parse.add_argument('-P','--depth_profile',type=str,default=default,help='the file to save depth profile [{}]'.format(default))
    default='nodes.snvs'
    parse.add_argument('-n','--nodes_snvs',type=str,default=default,help='the file to save snvs in each nodes [{}]'.format(default))
    default='log.txt'
    parse.add_argument('-g','--log',type=str,default=default,help='the log file [{}]'.format(default))
    default=10
    parse.add_argument('--loglevel',type=int,default=default,choices=[10,20,30,40,50],help='the logging level [{}]'.format(default))
    parse.add_argument('--trunk_vars',type=str,help='the trunk variants file')
    default='tree.dat'
    parse.add_argument('--tree_data',type=str,default=default,help='the file to dump the tree data [{}]'.format(default))
    default=None
    parse.add_argument('--expands',type=str,default=default,help='the basename of the file to output the snv and segment data for EXPANDS [{}]'.format(default))
    default=100000000
    parse.add_argument('--length',type=int,default=default,help='the length of sequence to simulate [{}]'.format(default))
    args=parse.parse_args()

    logging.basicConfig(filename=args.log, filemode='w', format='%(levelname)s: %(message)s', level=args.loglevel)
    logging.info(' Command: %s',' '.join(sys.argv))
    if args.random_seed==None:
        seed=numpy.random.randint(4294967296) #2**32: Must be convertible to 32 bit unsigned integers.
    else:
        seed=args.random_seed
    logging.info(' Random seed: %s',seed)
    numpy.random.seed(seed)

    with open(args.tree) as input:
        for line in input:
            newick=line.rstrip()
            mytree=tree.newick2tree(newick)
#            print(mytree.tree2newick())
            leaves_number=mytree.leaves_number()
            if args.prune>0:
                mytree.prune(tips=args.prune)
            elif args.prune_proportion>0.0:
                trim=leaves_number*args.prune_proportion
                mytree.prune(tips=trim)
#trunk vars
            trunk_snvs={}
            trunk_dels={}
            trunk_cnvs={}
            if args.trunk_vars!=None:
                trunk_snvs,trunk_dels,trunk_cnvs=trunk_vars.classify_vars(args.trunk_vars,args.ploid,leaves_number,mytree)

            snvs_freq,cnvs,depth_profile,nodes_snvs,tree_with_snvs=mytree.snvs_freq_cnvs_profile(
                                                       ploid=args.ploid,
                                                       snv_rate=args.snv_rate,
                                                       cnv_rate=args.cnv_rate,
                                                       del_prob=args.del_prob,
                                                       cnv_length_lambda=args.cnv_length_lambda,
                                                       cnv_length_max=args.cnv_length_max,
                                                       copy_max=args.copy_max,
                                                       trunk_snvs=trunk_snvs,
                                                       trunk_dels=trunk_dels,
                                                       trunk_cnvs=trunk_cnvs,
                                                       length=args.length,
                                                       )
            cnv_file=open(args.cnv,'w')
            for cnv in cnvs:
                cnv_file.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(cnv['seg'],cnv['start'],cnv['end'],cnv['copy'],cnv['leaves_count'],cnv['pre_snvs']))

            snv_file=open(args.snv,'w')
            for pos,freq in snvs_freq:
                snv_file.write('{}\t{}\n'.format(pos,freq))

            depth_profile_file=open(args.depth_profile,'w')
            for seg in depth_profile:
                depth_profile_file.write('{}\t{}\t{}\n'.format(*seg))

            nodes_snvs_file=open(args.nodes_snvs,'w')
            for node in sorted(nodes_snvs.keys()):
                for snv in sorted(nodes_snvs[node]):
                    nodes_snvs_file.write('{}\t{}\n'.format(node,snv))

#TODO: Should we change pickle to json or yaml?
#http://stackoverflow.com/questions/4677012/python-cant-pickle-type-x-attribute-lookup-failed
            tree_data_file=open(args.tree_data,'wb')
            #pickle.dump(tree_with_snvs,tree_data_file)

            if args.expands != None:
#output for expands
                expands_snps_file=open(args.expands+'.snps','w')
                expands_snps_file.write('chr\tstartpos\tAF_Tumor\tPN_B\n')
                chroms=1
                for pos,freq in snvs_freq:
                    total_dp,b_allele_dp=simulate_sequence_coverage(args.depth,freq)
                    expands_snps_file.write('{}\t{}\t{}\t{}\n'.format(chroms,pos,b_allele_dp/total_dp,0))

#in the segment input for expands
#CN_Estimate - the copy number estimated for each segment (average value across all subpopulations in the sample)
                expands_segs_file=open(args.expands+'.segs','w')
                expands_segs_file.write("chr\tstartpos\tendpos\tCN_Estimate\n")
                for start,end,copy in depth_profile:
                    expands_segs_file.write('{}\t{}\t{}\t{}\n'.format(chroms,start,end,copy/leaves_number))

if __name__ == '__main__':
    __main__()
