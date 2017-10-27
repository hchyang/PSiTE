#!/usr/bin/env python3

#########################################################################
# Author: Hechuan
# Created Time: 2017-04-04 18:00:34
# File Name: allinone.py
# Description: 
#########################################################################

import sys
import os
import argparse
import numpy
import yaml
import logging
import subprocess
import pyfaidx
from csite.phylovar import check_prune,check_proportion,check_seed,random_int,check_config_file

#handle the error below
#python | head == IOError: [Errno 32] Broken pipe 
from signal import signal, SIGPIPE, SIG_DFL 
signal(SIGPIPE,SIG_DFL) 

def main(progname=None):
    parse=argparse.ArgumentParser(
        description='a wrapper of simulating short reads from genome with germline and somatic variants',
        prog=progname if progname else sys.argv[0])
    parse.add_argument('-r','--reference',required=True,
        help='a fasta file of the reference genome')
    parse.add_argument('-v','--vcf',required=True,
        help='a vcf file contains germline variants')
    parse.add_argument('-t','--tree',required=True,
        help='a newick file contains ONE tree')
    parse.add_argument('-c','--config',required=True,
        help='a yaml file contains configure for somatic variant simulation')
    parse.add_argument('-o','--output',required=True,
        help='output folder')
    default=None
    parse.add_argument('--trunk_vars',type=str,
        help='the trunk variants file supplied by user [{}]'.format(default))
    default=0
    parse.add_argument('--trunk_length',type=float,
        help='the length of the truncal branch [{}]'.format(default))
    default=0
    parse.add_argument('-x','--prune',type=check_prune,default=default,
        help='trim all the children of the nodes with equal or less than this number of tips [{}]'.format(default))
    default=0.0
    parse.add_argument('-X','--prune_proportion',type=check_proportion,default=default,
        help='trim all the children of the nodes with equal or less than this proportion of tips [{}]'.format(default))
    default=50
    parse.add_argument('-d','--depth',type=int,default=default,
        help='the mean depth for ART to simulate short reads [{}]'.format(default))
    default=0.5
    parse.add_argument('-p','--purity',type=int,default=default,
        help='the proportion of tumor cells in simulated sample [{}]'.format(default))
    default='art_illumina --paired --len 100 --mflen 500 --sdev 20'
    parse.add_argument('--art',type=str,default=default,
        help='the parameters for ART program [{}]'.format(default))
    default=None
    parse.add_argument('-s','--random_seed',type=check_seed,
        help='the seed for random number generator [{}]'.format(default))
    default='allinone.log'
    parse.add_argument('-g','--log',type=str,default=default,
        help='the log file to save the settings of each command [{}]'.format(default))
    args=parse.parse_args()
    with open(args.config,'r') as configfile:
        config=yaml.safe_load(configfile)
    check_config_file(config=config)

#get absolute paths for the input files
    reference=os.path.abspath(args.reference)
    vcf=os.path.abspath(args.vcf)
    tree=os.path.abspath(args.tree)
    config=os.path.abspath(args.config)
    if args.trunk_vars:
        trunk_vars=os.path.abspath(args.trunk_vars)
    outdir=args.output
    try:
        os.mkdir(outdir)
    except FileExistsError:
        exit('{} already exists. Try another directory to output! (-o/--output)'.format(outdir))
    os.chdir(outdir)

###### logging and random seed setting
    logging.basicConfig(filename=args.log, filemode='w',
        format='[%(asctime)s] %(levelname)s: %(message)s',
        datefmt='%m-%d %H:%M:%S',level='INFO')
    logging.info(' Command: %s',' '.join(sys.argv[0:1]+['allinone']+sys.argv[1:]))
    if args.random_seed==None:
        seed=random_int()
    else:
        seed=args.random_seed
    logging.info(' Random seed: %s',seed)
    numpy.random.seed(seed)

#vcf2fa
    normal_fa='normal_fa'
    cmd_params=[sys.argv[0],'vcf2fa',
                '--vcf',vcf,
                '--reference',reference,
                '--output',normal_fa]
    logging.info(' Command: %s',' '.join(cmd_params))
    subprocess.run(args=cmd_params,check=True)

#phylovar
    cmd_params=[sys.argv[0],'phylovar',
                '--tree',tree,
                '--config',config,
                '--random_seed',str(random_int()),
                '--chain','tumor_chain']
    if args.trunk_vars:
        cmd_params.extend(['--trunk_vars',trunk_vars])
    if args.trunk_length:
        cmd_params.extend(['--trunk_length',str(args.trunk_length)])
    if args.prune:
        cmd_params.extend(['--prune',str(args.prune)])
    elif args.prune_proportion:
        cmd_params.extend(['--prune_proportion',str(args.prune_proportion)])
    logging.info(' Command: %s',' '.join(cmd_params))
    subprocess.run(args=cmd_params,check=True)

#chain2fa
    tumor_fa='tumor_fa'
    cmd_params=[sys.argv[0],'chain2fa',
                '--chain','tumor_chain',
                '--reference','{dir}/normal_hap0.fa,{dir}/normal_hap1.fa'.format(dir=normal_fa),
                '--output',tumor_fa]
    logging.info(' Command: %s',' '.join(cmd_params))
    subprocess.run(args=cmd_params,check=True)

#compute coverage and run ART
#FIXME: cell number: float? int?
    normal_gsize=0
    normal_gsize+=genomesize(fasta='{}/normal_hap0.fa'.format(normal_fa))
    normal_gsize+=genomesize(fasta='{}/normal_hap1.fa'.format(normal_fa))
    total_seq_bases=normal_gsize/2*args.depth
    #print(total_seq_bases)

    tip_node_leaves=tip_node_leaves_counting(f='tumor_chain/tip_node_sample.count')
    tumor_cells=sum(tip_node_leaves.values())
    normal_cells=tumor_cells/args.purity-tumor_cells
    #print(normal_cells)

    normal_dna=normal_gsize*normal_cells
    tumor_dna=0
    tip_node_gsize={}
    for tip_node,leaves in tip_node_leaves.items():
        tip_node_gsize[tip_node]=genomesize(fasta='{}/{}.genome.fa'.format(tumor_fa,tip_node))
        tumor_dna+=tip_node_gsize[tip_node]*tip_node_leaves[tip_node]
    seq_per_base=total_seq_bases/(normal_dna+tumor_dna)

    art_folder='art_reads'
    os.mkdir(art_folder)
    art_params=args.art.split()

#two normal cell haplotypes
    cmd_params=art_params[:]
    cmd_params.extend(['--fcov',str(normal_cells*seq_per_base)])
    for hap in 0,1:
        ref='{}/normal_hap{}.fa'.format(normal_fa,hap)
        prefix='{}/normal_hap{}.'.format(art_folder,hap)
        final_cmd_params=cmd_params+['--in',ref,'--out',prefix,'--rndSeed',str(random_int())]
        logging.info(' Command: %s',' '.join(final_cmd_params))
        subprocess.run(args=final_cmd_params,check=True)

#tumor cells haplotypes
    for tip_node in sorted(tip_node_leaves.keys()):
        cmd_params=art_params[:]
        cmd_params.extend(['--fcov',str(tip_node_leaves[tip_node]*seq_per_base)])
        ref='{}/{}.genome.fa'.format(tumor_fa,tip_node)
        prefix='{}/{}.'.format(art_folder,tip_node)
        final_cmd_params=cmd_params+['--in',ref,'--out',prefix,'--rndSeed',str(random_int())]
        logging.info(' Command: %s',' '.join(final_cmd_params))
        subprocess.run(args=final_cmd_params,check=True)

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

