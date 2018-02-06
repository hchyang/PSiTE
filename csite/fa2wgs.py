#!/usr/bin/env python3

#########################################################################
# Author: Hechuan Yang
# Created Time: 2017-04-04 18:00:34
# File Name: fa2wgs.py
# Description: 
#########################################################################

import sys
import os
import glob
import argparse
import numpy
import math
import logging
import pyfaidx
import subprocess
import multiprocessing
from csite.phylovar import check_seed,random_int

#handle the error below
#python | head == IOError: [Errno 32] Broken pipe 
from signal import signal, SIGPIPE, SIG_DFL 
signal(SIGPIPE,SIG_DFL) 

def check_folder(directory=None):
    if not os.path.isdir(directory):
        raise argparse.ArgumentTypeError("'{}' doesn't exist or isn't a folder.".format(directory))
    return directory
    
def check_file(directory=None):
    if not os.path.isfile(directory):
        raise argparse.ArgumentTypeError("'{}' doesn't exist or isn't a file.".format(directory))
    return directory
    
def check_purity(value=None):
    fvalue=float(value)
    if not 0<fvalue<=1: 
        raise argparse.ArgumentTypeError("{} is an invalid value for --purity. ".format(value)+
            "It should be a float number in the range of (0,1].")
    return fvalue

def check_depth(value=None):
    fvalue=float(value)
    if fvalue<0: 
        raise argparse.ArgumentTypeError("{} is an invalid value for --depth/--normal_depth. ".format(value)+
            "It should be a non-negative float number.")
    return fvalue

def main(progname=None):
    parse=argparse.ArgumentParser(
        description='A wrapper of simulating WGS reads from normal and tumor genome fasta',
        prog=progname if progname else sys.argv[0])
    parse.add_argument('-n','--normal',type=check_folder,required=True,metavar='DIR',
        help='the directory of the normal fasta')
    parse.add_argument('-t','--tumor',type=check_folder,required=True,metavar='DIR',
        help='the directory of the tumor fasta')
    parse.add_argument('-m','--map',type=check_file,required=True,metavar='FILE',
        help='the map file containing the relationship between tip nodes and samples')
    default='art_reads'
    parse.add_argument('-o','--output',type=str,default=default,metavar='DIR',
        help='output directory [{}]'.format(default))
    default=50
    parse.add_argument('-d','--depth',type=check_depth,default=default,metavar='FLOAT',
        help='the mean depth of tumor sample for ART to simulate NGS reads [{}]'.format(default))
    default=0
    parse.add_argument('-D','--normal_depth',type=check_depth,default=default,metavar='FLOAT',
        help='the mean depth of normal sample for ART to simulate NGS reads [{}]'.format(default))
    default=0.5
    parse.add_argument('-p','--purity',type=check_purity,default=default,metavar='FLOAT',
        help='the proportion of tumor cells in simulated tumor sample [{}]'.format(default))
    default=None
    parse.add_argument('-s','--random_seed',type=check_seed,metavar='INT',
        help='the seed for random number generator [{}]'.format(default))
    default='fa2wgs.log'
    parse.add_argument('-g','--log',type=str,default=default,metavar='FILE',
        help='the log file to save the settings of each command [{}]'.format(default))
    default='art_illumina --noALN --quiet --paired --len 100 --mflen 500 --sdev 20'
    parse.add_argument('--art',type=str,default=default,metavar='STR',
        help='the parameters for ART program [{}]'.format(default))
    default=1
    parse.add_argument('--cores',type=int,default=default,metavar='INT',
        help='number of cores used to run the program [{}]'.format(default))
    parse.add_argument('--compress',action="store_true",
        help='compress the generated fastq files using gzip')
    parse.add_argument('--single',action="store_true",
        help="single cell mode. Output each tip node's NGS reads file seperately")
    args=parse.parse_args()

# logging and random seed setting
    logging.basicConfig(filename=args.log, 
        filemode='w',format='[%(asctime)s] %(levelname)s: %(message)s',
        datefmt='%m-%d %H:%M:%S',level='INFO')
    argv_copy=sys.argv[:]
    try:
        art_index=argv_copy.index('--art')
        argv_copy[art_index+1]="'{}'".format(argv_copy[art_index+1])
    except ValueError:
        pass
    argv_copy.insert(1,'fa2wgs')
    logging.info(' Command: %s',' '.join(argv_copy))
    if args.random_seed==None:
        seed=random_int()
    else:
        seed=args.random_seed
    logging.info(' Random seed: %s',seed)
    numpy.random.seed(seed)

#exit the program if you do not want to simulate any reads for normal or tumor samples
    if args.depth+args.normal_depth==0:
        sys.exit('Do nothing as the total depth (normal+tumor) is 0!')

#create index file (.fai) for each fasta
    tip_node_leaves=tip_node_leaves_counting(f=args.map)
    pool=multiprocessing.Pool(processes=args.cores)
    for parental in 0,1:
        fasta='{}/normal.parental_{}.fa'.format(args.normal,parental)
        assert os.path.isfile(fasta),\
            "Couldn't find {} under the normal directory: {}".format(fasta,args.normal)
        pool.apply_async(pyfaidx.Faidx,args=(fasta,))
        for tip_node in tip_node_leaves.keys():
            fasta='{}/{}.parental_{}.fa'.format(args.tumor,tip_node,parental)
            assert os.path.isfile(fasta),\
                "Couldn't find {} under the tumor directory: {}".format(fasta,args.tumor)
            pool.apply_async(pyfaidx.Faidx,args=(fasta,))
    pool.close()
    pool.join()

#compute coverage and run ART
#FIXME: cell number: float? int?
    normal_gsize=0
    for parental in 0,1:
        normal_gsize+=genomesize(fasta='{}/normal.parental_{}.fa'.format(args.normal,parental))
    tumor_seq_bases=normal_gsize/2*args.depth

    tumor_cells=sum(tip_node_leaves.values())
    total_cells=tumor_cells/args.purity
    normal_cells=total_cells-tumor_cells

    normal_dna=normal_gsize*normal_cells
    tumor_dna=0
    tip_node_gsize={}
    for tip_node,leaves in tip_node_leaves.items():
#The value of tip_node_gsize[tip_node] is a list of three elements:
#0)genomesize of parental 0
#1)genomesize of parental 1
#2)the sum of parental 0 and 1
        tip_node_gsize[tip_node]=[]
        for parental in 0,1:
            tip_node_gsize[tip_node].append(genomesize(fasta='{}/{}.parental_{}.fa'.format(args.tumor,tip_node,parental)))
        tip_node_gsize[tip_node].append(tip_node_gsize[tip_node][0]+tip_node_gsize[tip_node][1])
        tumor_dna+=tip_node_gsize[tip_node][2]*tip_node_leaves[tip_node]
    tumor_seq_per_base=tumor_seq_bases/(normal_dna+tumor_dna)

#create output folders
    if os.path.exists(args.output):
        if os.path.isdir(args.output):
            pass
        else:
            raise OutputExistsError("A file in the name of '{}' exists.\nDelete it or try another name as output folder.".format(args.output))
    else:
        os.mkdir(args.output,mode=0o755)
    tumor_dir=args.output+'/tumor'
    normal_dir=args.output+'/normal'

#collect simulation parameters first
    params_matrix=[]
    art_params=args.art

#simulation for normal sample
    if args.normal_depth>0:
        try:
            os.mkdir(normal_dir,mode=0o755)
        except FileExistsError as e:
            raise OutputExistsError("'{}' exists already! \nCan not use it as the output folder of normal NGS reads.".format(normal_dir)+
                '\nDelete it or use another folder as output folder.') from e
        for parental in 0,1:
            prefix='{}/normal.parental_{}.'.format(normal_dir,parental)
            fcov=args.normal_depth/2
            ref='{}/normal.parental_{}.fa'.format(args.normal,parental)
            sim_cfg={
                'gsize':normal_gsize/2,
                'base_cmd':art_params,
                'fcov':fcov,
                'in':ref,
                'out':prefix,
                'id':'nm_prt{}'.format(parental)}
            params_matrix.append(sim_cfg)

#simulation for tumor sample
    if args.depth>0:
        try:
            os.mkdir(tumor_dir,mode=0o755)
        except FileExistsError as e:
            raise OutputExistsError("'{}' exists already! \nCan not use it as the output folder of tumor NGS reads.".format(normal_dir)+
                '\nDelete it or use another folder as output folder.') from e
#create a reference meta file which can be used by wessim to simulate exome-seq data
        ref_meta=open('reference.meta','w')
#two normal cell haplotypes
        for parental in 0,1:
            prefix='{}/normal.parental_{}.'.format(tumor_dir,parental)
            fcov=normal_cells*tumor_seq_per_base
            ref='{}/normal.parental_{}.fa'.format(args.normal,parental)
            sim_cfg={
                'gsize':normal_gsize/2,
                'base_cmd':art_params,
                'fcov':fcov,
                'in':ref,
                'out':prefix,
                'id':'nm_prt{}'.format(parental)}
            params_matrix.append(sim_cfg)
            fullname=os.path.abspath(ref)
            ref_meta.write('{}\t{}\n'.format(fullname,str(normal_cells/total_cells/2)))

#tumor cells haplotypes
        for tip_node in sorted(tip_node_leaves.keys()):
            fcov=tip_node_leaves[tip_node]*tumor_seq_per_base
            for parental in 0,1:
                ref='{}/{}.parental_{}.fa'.format(args.tumor,tip_node,parental)
                prefix='{}/{}.parental_{}.'.format(tumor_dir,tip_node,parental)
                sim_cfg={
                    'gsize':tip_node_gsize[tip_node][parental],
                    'base_cmd':art_params,
                    'fcov':fcov,
                    'in':ref,
                    'out':prefix,
                    'id':'{}_prt{}'.format(tip_node,parental)}
                params_matrix.append(sim_cfg)
                fullname=os.path.abspath(ref)
                ref_meta.write('{}\t{}\n'.format(fullname,str(tip_node_leaves[tip_node]/total_cells*tip_node_gsize[tip_node][parental]/tip_node_gsize[tip_node][2])))
        ref_meta.close()

#generate fastq (and compress them) parallelly
#every thread will generate at most 2 percent of the total data you want to simulate
    sizeBlock=normal_gsize*(args.depth+args.normal_depth)*0.02
    final_params_matrix=[]
    for cfg in params_matrix:
        n=math.ceil(cfg['gsize']*cfg['fcov']/sizeBlock)
        if n==0:
            continue
        cfg['fcov']/=n
        for i in range(n):
            final_params_matrix.append(cfg.copy())
            final_params_matrix[-1]['out']=cfg['out']+'{:02d}.'.format(i)
            final_params_matrix[-1]['id']=cfg['id']+'_{:02d}-'.format(i)
            final_params_matrix[-1]['rndSeed']=str(random_int())
    pool=multiprocessing.Pool(processes=args.cores)
    for x in final_params_matrix:
        pool.apply_async(generate_fq,args=(x,args.compress))
    pool.close()
    pool.join()

#merge small fastq files into one for normal/tumor sample
    sample_fq_files=[]
    suffixes=['fq','1.fq','2.fq']
    if args.compress:
        suffixes=[x+'.gz' for x in suffixes]
    if args.normal_depth>0:
        for suffix in suffixes:
            prefix='{}/normal.parental_[01].[0-9][0-9].'.format(normal_dir)
            source=glob.glob(prefix+suffix)
            if len(source):
                target='{}/normal.{}'.format(normal_dir,suffix)
                source.sort()
                sample_fq_files.append([target,source])
    if args.depth>0:
        for suffix in suffixes:
            if args.single:
                for tip_node in ['normal']+sorted(tip_node_leaves.keys()):
                    prefix='{}/{}.parental_[01].[0-9][0-9].'.format(tumor_dir,tip_node)
                    source=glob.glob(prefix+suffix)
                    if len(source):
                        target='{}/{}.{}'.format(tumor_dir,tip_node,suffix)
                        source.sort()
                        sample_fq_files.append([target,source])
            else:
                prefix='{}/*.parental_[01].[0-9][0-9].'.format(tumor_dir)
                source=glob.glob(prefix+suffix)
                if len(source):
                    target='{}/tumor.{}'.format(tumor_dir,suffix)
                    source.sort()
                    sample_fq_files.append([target,source])
    pool=multiprocessing.Pool(processes=args.cores)
    for x in sample_fq_files:
        pool.apply_async(merge_fq,args=x)
    pool.close()
    pool.join()

def merge_fq(target=None,source=None):
    '''
    after generating fq in multiprocessing mode,
    there will be multiple fq files for each genome.
    I will merge them into one file for each genome.
    '''
    assert not os.path.isfile(target),"'{}' exists already!"
    with open(target,'a') as output:
        for f in source:
            subprocess.run(args=['cat',f],check=True,stdout=output)
            subprocess.run(args=['rm',f],check=True)

def generate_fq(params=None,compress=False):
    '''
    run art command to generate the fastq file, and call compress_fq if required.
    '''
    cmd_params=params['base_cmd'].split()+['--fcov',str(params['fcov']),
                                           '--in',params['in'],
                                           '--id',params['id'],
                                           '--out',params['out'],
                                           '--rndSeed',params['rndSeed']]
    subprocess.run(args=cmd_params,check=True)
    logging.info(' Command: {}'.format(' '.join(cmd_params)))
    if compress:
        compress_fq(prefix=params['out'])

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
                tip_node,leaves=line.split()[:2]
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

class OutputExistsError(Exception):
    pass
