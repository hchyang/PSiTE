#!/usr/bin/env python3

#########################################################################
# Author: Hechuan
# Created Time: 2017-04-04 18:00:34
# File Name: fa2ngs.py
# Description: 
#########################################################################

import sys
import os
import argparse
import numpy
import logging
import pyfaidx
import subprocess
from csite.phylovar import check_seed,random_int

#handle the error below
#python | head == IOError: [Errno 32] Broken pipe 
from signal import signal, SIGPIPE, SIG_DFL 
signal(SIGPIPE,SIG_DFL) 

def main(progname=None):
    parse=argparse.ArgumentParser(
        description='a wrapper of simulating short reads from normal and tumor genome fasta',
        prog=progname if progname else sys.argv[0])
    parse.add_argument('-n','--normal',required=True,
        help='the directory of the normal fasta')
    parse.add_argument('-t','--tumor',required=True,
        help='the directory of the tumor fasta')
    parse.add_argument('-c','--chain',required=True,
        help='the directory of the tumor chain')
    default='art_reads'
    parse.add_argument('-o','--output',type=str,default=default,
        help='output directory [{}]'.format(default))
    default=50
    parse.add_argument('-d','--depth',type=float,default=default,
        help='the mean depth of tumor for ART to simulate short reads [{}]'.format(default))
    default=50
    parse.add_argument('-D','--normal_depth',type=float,default=default,
        help='the mean depth of normal for ART to simulate short reads [{}]'.format(default))
    default=0.5
    parse.add_argument('-p','--purity',type=float,default=default,
        help='the proportion of tumor cells in simulated sample [{}]'.format(default))
    default=None
    parse.add_argument('-s','--random_seed',type=check_seed,
        help='the seed for random number generator [{}]'.format(default))
    default='fa2ngs.log'
    parse.add_argument('-g','--log',type=str,default=default,
        help='the log file to save the settings of each command [{}]'.format(default))
    default='art_illumina --noALN --quiet --paired --len 100 --mflen 500 --sdev 20'
    parse.add_argument('--art',type=str,default=default,
        help='the parameters for ART program [{}]'.format(default))
    args=parse.parse_args()

###### logging and random seed setting
    logging.basicConfig(filename=args.log, 
        filemode='w',format='[%(asctime)s] %(levelname)s: %(message)s',
        datefmt='%m-%d %H:%M:%S',level='INFO')
    logging.info(' Command: %s',' '.join(sys.argv))
    if args.random_seed==None:
        seed=random_int()
    else:
        seed=args.random_seed
    logging.info(' Random seed: %s',seed)
    numpy.random.seed(seed)

#there must be two haplotype fasta in the normal dir
    assert os.path.isdir(args.normal),'{} is not exist or not a folder.'.format(args.normal)
    for hap in 0,1:
        assert os.path.isfile('{}/normal_parental_{}.fa'.format(args.normal,hap)), 'Can not find normal_parental_{}.fa under the normal directory: {}'.format(hap,args.normal)

#tumor directory and chain directory must exist.
#also file chain_dir/tip_node_sample.count.
    assert os.path.isdir(args.chain),'{} is not exist or not a folder.'.format(args.chain)
    assert os.path.isfile(args.chain+'/tip_node_sample.count'),'Can not find tip_node_sample.count under the chain directory: {}'.format(args.chain)
    assert os.path.isdir(args.tumor),'{} is not exist or not a folder.'.format(args.chain)

#compute coverage and run ART
#FIXME: cell number: float? int?
    normal_gsize=0
    for hap in 0,1:
        normal_gsize+=genomesize(fasta='{}/normal_parental_{}.fa'.format(args.normal,hap))
    total_seq_bases=normal_gsize/2*args.depth

    tip_node_leaves=tip_node_leaves_counting(f='{}/tip_node_sample.count'.format(args.chain))
    tumor_cells=sum(tip_node_leaves.values())
    total_cells=tumor_cells/args.purity
    normal_cells=total_cells-tumor_cells

    normal_dna=normal_gsize*normal_cells
    tumor_dna=0
    tip_node_gsize={}
    for tip_node,leaves in tip_node_leaves.items():
        assert os.path.isfile('{}/{}.genome.fa'.format(args.tumor,tip_node)), 'Can not find {}.genome.fa under the tumor directory: {}'.format(tip_node,args.tumor)
        tip_node_gsize[tip_node]=genomesize(fasta='{}/{}.genome.fa'.format(args.tumor,tip_node))
        tumor_dna+=tip_node_gsize[tip_node]*tip_node_leaves[tip_node]
    seq_per_base=total_seq_bases/(normal_dna+tumor_dna)

    tumor_dir=args.output+'/tumor'
    normal_dir=args.output+'/normal'
    os.mkdir(args.output)
    os.mkdir(tumor_dir)
    os.mkdir(normal_dir)
    art_params=args.art.split()

#create a reference meta file which can be used by wessim to simulate exome-seq data
    ref_meta=open('reference.meta','w')

#two normal cell haplotypes
    cmd_params=art_params[:]
    for hap in 0,1:
        fcov=str(normal_cells*seq_per_base)
        ref='{}/normal_parental_{}.fa'.format(args.normal,hap)
        prefix='{}/normal_parental_{}.'.format(tumor_dir,hap)
        final_cmd_params=cmd_params+['--fcov',fcov,'--in',ref,'--out',prefix,'--id','n_hap{}'.format(hap),'--rndSeed',str(random_int())]
        logging.info(' Command: %s',' '.join(final_cmd_params))
        subprocess.run(args=final_cmd_params,check=True)
        compress_fq(prefix=prefix)

        if args.normal_depth>0:
            fcov=str(args.normal_depth/2)
            prefix='{}/normal_parental_{}.'.format(normal_dir,hap)
            final_cmd_params=cmd_params+['--fcov',fcov,'--in',ref,'--out',prefix,'--id','n_hap{}'.format(hap),'--rndSeed',str(random_int())]
            logging.info(' Command: %s',' '.join(final_cmd_params))
            subprocess.run(args=final_cmd_params,check=True)
            compress_fq(prefix=prefix)

        fullname=os.path.abspath(ref)
        ref_meta.write('{}\t{}\n'.format(fullname,str(normal_cells/total_cells/2)))

#tumor cells haplotypes
    for tip_node in sorted(tip_node_leaves.keys()):
        fcov=str(tip_node_leaves[tip_node]*seq_per_base)
        ref='{}/{}.genome.fa'.format(args.tumor,tip_node)
        prefix='{}/{}.'.format(tumor_dir,tip_node)
        final_cmd_params=cmd_params+['--fcov',fcov,'--in',ref,'--out',prefix,'--id','t_{}'.format(tip_node),'--rndSeed',str(random_int())]
        logging.info(' Command: %s',' '.join(final_cmd_params))
        subprocess.run(args=final_cmd_params,check=True)
        compress_fq(prefix=prefix)

        fullname=os.path.abspath(ref)
        ref_meta.write('{}\t{}\n'.format(fullname,str(tip_node_leaves[tip_node]/total_cells)))

    ref_meta.close()

def compress_fq(prefix=None):
    suffixes=['fq','1.fq','2.fq']
    for suffix in suffixes:
        fq=prefix+suffix
        if os.path.isfile(fq):
            subprocess.run(args=['gzip',fq],check=True)

def tip_node_leaves_counting(f=None):
    '''
    Return a dictionay with structure: 
    {tip_node1:leaves_count1,tip_node2:leaves_count2,...}
    '''
    tip_node_leaves={}
    with open(f,'r') as input:
        for line in input:
            if not line.startswith('#'):
                tip_node,leaves=line.split()
                tip_node_leaves[tip_node]=int(leaves)
    return tip_node_leaves

def genomesize(fasta=None):
    '''
    Extract genome size from .fa file.
    '''
    fa=pyfaidx.Faidx(fasta)
    gsize=0
    for chroms in fa.index.keys():
        gsize+=fa.index[chroms].rlen
    return gsize

