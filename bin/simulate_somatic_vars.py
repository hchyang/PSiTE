#!/usr/bin/env python3

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

def check_seed(value):
    ivalue=int(value)
    if not 0<=ivalue<=4294967296: #2**32: Must be convertible to 32 bit unsigned integers.
         raise argparse.ArgumentTypeError("{} is an invalid value for random seed. It should be an interger between 0 and 4294967296.".format(value))
    return ivalue

def check_max_cnv_length(seq_length,max_cnv_length):
    if max_cnv_length>seq_length:
         raise argparse.ArgumentTypeError("The value of -L/--cnv_length_max ({}) should NOT be larger than --length ({}).".format(max_cnv_length,seq_length))

def cn_dist(copy_max=None,copy_parameter=None):
    '''
    When an amplification event happens, it will randomly pick a copy number according to a geometric like distribution,
    p(n+1)=w*p(n).
    This function returns the configure of the distribution of CNVs' copy number.
    '''
    cn_dist_cfg={}
    cn_dist_cfg['copy']=[x for x in range(1,copy_max+1)]
    w=1/sum([copy_parameter**x for x in range(copy_max)])
    scale=sum([w*copy_parameter**x for x in range(copy_max)])
    cn_dist_cfg['prob']=[(w*copy_parameter**x)/scale for x in range(copy_max)]
    print(cn_dist_cfg['copy'])
    print(cn_dist_cfg['prob'])
    return cn_dist_cfg

#@profile
def __main__():
    parse=argparse.ArgumentParser(description='Simulate SNVs/CNVs on a coalescent tree in newick format')
    parse.add_argument('-t','--tree',required=True,help='a tree in newick format')
    default=300
    parse.add_argument('-r','--snv_rate',type=float,default=default,help='the muation rate of SNVs [{}]'.format(default))
    default=3
    parse.add_argument('-R','--cnv_rate',type=float,default=default,help='the muation rate of CNVs [{}]'.format(default))
    default=0.5
    parse.add_argument('-d','--del_prob',type=int,default=default,help='the probability of being deletion for a CNV mutation [{}]'.format(default))
#https://en.wikipedia.org/wiki/Copy-number_variation
    default=120000
    parse.add_argument('-l','--cnv_length_beta',type=int,default=default,help='the mean of CNVs length [{}]'.format(default))
    default=20000000
    parse.add_argument('-L','--cnv_length_max',type=int,default=default,help='the maximium of CNVs length [{}]'.format(default))
    default=5
    parse.add_argument('-c','--copy_max',type=int,default=default,help='the maximium ADDITIONAL copy of a CNVs [{}]'.format(default))
    default=0.5
    parse.add_argument('-C','--copy_parameter',type=float,default=default,help="the p parameter of CNVs' copy number distribution [{}]".format(default))
    default=2
    parse.add_argument('-P','--ploidy',type=int,default=default,help='the ploidy to simulate [{}]'.format(default))
    default=1.0
    parse.add_argument('-p','--purity',type=float,default=default,help="the purity of tumor cells in the simulated sample [{}]".format(default))
    default=50
    parse.add_argument('-D','--depth',type=int,default=default,help='the mean depth for simulating coverage data [{}]'.format(default))
    default=0
    parse.add_argument('-x','--prune',type=int,default=default,help='trim all their children for the branches with equal or less than this number of tips [{}]'.format(default))
    default=0.0
    parse.add_argument('-X','--prune_proportion',type=float,default=default,help='trim all their children for the branches with equal or less than this proportion of tips [{}]'.format(default))
    default=None
    parse.add_argument('-s','--random_seed',type=check_seed,help='the seed for random number generator [{}]'.format(default))
    default='output.snvs'
    parse.add_argument('-S','--snv',type=str,default=default,help='the output file to save SNVs [{}]'.format(default))
    default='output.cnvs'
    parse.add_argument('-V','--cnv',type=str,default=default,help='the output file to save CNVs [{}]'.format(default))
    default='output.nodes.snvs'
    parse.add_argument('-n','--nodes_snvs',type=str,default=default,help='the file to save SNVs on each nodes [{}]'.format(default))
    default='log.txt'
    parse.add_argument('-g','--log',type=str,default=default,help='the log file [{}]'.format(default))
    default='DEBUG'
    parse.add_argument('-G','--loglevel',type=str,default=default,choices=['DEBUG','INFO'],help='the logging level [{}]'.format(default))
    default='output.cnv.profile'
    parse.add_argument('--cnv_profile',type=str,default=default,help='the file to save CNVs profile [{}]'.format(default))
    parse.add_argument('--trunk_vars',type=str,help='the trunk variants file')
    parse.add_argument('--snv_genotype',type=str,help='the file to save SNV genotypes for each sample')
    parse.add_argument('--ind_cnvs',type=str,help='the file to save CNVs for each sample individual')
    default=0
    parse.add_argument('--trunk_length',type=float,help='the trunk length [{}]'.format(default))
    default='tree.dat'
    parse.add_argument('--tree_data',type=str,default=default,help='the file to dump the tree data [{}]'.format(default))
    default=None
    parse.add_argument('--expands',type=str,default=default,help='the basename of the file to output the snv and segment data for EXPANDS [{}]'.format(default))
    default=100000000
    parse.add_argument('--length',type=int,default=default,help='the length of the sequence to simulate [{}]'.format(default))
    args=parse.parse_args()

    check_max_cnv_length(args.length,args.cnv_length_max)

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
            if args.trunk_length:
                mytree.lens=args.trunk_length
            leaves_number=mytree.leaves_counting()
            leaves_names=sorted(mytree.leaves_naming())
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
                trunk_snvs,trunk_dels,trunk_cnvs=trunk_vars.classify_vars(args.trunk_vars,args.ploidy,args.length,leaves_number,mytree)

            cn_dist_cfg=cn_dist(copy_max=args.copy_max,copy_parameter=args.copy_parameter)

            snvs_freq,cnvs,cnv_profile,nodes_snvs,tree_with_snvs,leaf_snv_alts,leaf_snv_refs,leaf_cnvs=mytree.snvs_freq_cnvs_profile(
                ploidy=args.ploidy,
                snv_rate=args.snv_rate,
                cnv_rate=args.cnv_rate,
                del_prob=args.del_prob,
                cnv_length_beta=args.cnv_length_beta,
                cnv_length_max=args.cnv_length_max,
                cn_dist_cfg=cn_dist_cfg,
                trunk_snvs=trunk_snvs,
                trunk_dels=trunk_dels,
                trunk_cnvs=trunk_cnvs,
                purity=args.purity,
                length=args.length,
                )

            if args.snv_genotype!=None:
                genotype_file=open(args.snv_genotype,'w')
                genotype_file.write('{}\t{}\n'.format('#Positon','\t'.join(leaves_names)))
                for snv in snvs_freq:
                    genotype_file.write('{}\t{}\n'.format(snv[0],
                        '\t'.join([str(leaf_snv_alts[leaf][snv[0]])+':'+str(leaf_snv_refs[leaf][snv[0]]) for leaf in leaves_names])))

            if args.ind_cnvs!=None:
                ind_cnvs_file=open(args.ind_cnvs,'w')
                ind_cnvs_file.write('#Sample\tParental\tStart\tEnd\tCopy\n')
                for leaf in sorted(leaf_cnvs.keys()):
                    for cnv in leaf_cnvs[leaf]:
                        ind_cnvs_file.write('{}\n'.format('\t'.join([str(x) for x in [leaf,cnv['parental'],cnv['start'],cnv['end'],cnv['copy']]])))

            cnv_file=open(args.cnv,'w')
            for cnv in cnvs:
                cnv_file.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(cnv['seg'],cnv['start'],cnv['end'],cnv['copy'],cnv['leaves_count'],cnv['pre_snvs']))

            snv_file=open(args.snv,'w')
            snv_file.write('#Position\tTrue_freq\tTotal_depth\tSimulated_freq\n')
            for pos,freq in snvs_freq:
                total_dp,b_allele_dp=tree.simulate_sequence_coverage(args.depth,freq)
                b_allele_freq=0
                if total_dp!=0:
                    b_allele_freq=b_allele_dp/total_dp
                snv_file.write('{}\t{}\t{}\t{}\n'.format(pos,freq,total_dp,b_allele_freq))

            cnv_profile_file=open(args.cnv_profile,'w')
            for seg in cnv_profile:
                cnv_profile_file.write('{}\t{}\t{}\n'.format(*seg))

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
                    total_dp,b_allele_dp=tree.simulate_sequence_coverage(args.depth,freq)
                    expands_snps_file.write('{}\t{}\t{}\t{}\n'.format(chroms,pos,b_allele_dp/total_dp,0))

#in the segment input for expands
#CN_Estimate - the copy number estimated for each segment (average value across all subpopulations in the sample)
                expands_segs_file=open(args.expands+'.segs','w')
                expands_segs_file.write("chr\tstartpos\tendpos\tCN_Estimate\n")
                for start,end,copy in cnv_profile:
                    expands_segs_file.write('{}\t{}\t{}\t{}\n'.format(chroms,start,end,copy/leaves_number))

if __name__ == '__main__':
    __main__()
