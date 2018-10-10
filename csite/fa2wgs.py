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
import shutil
import gzip
import time
from csite.phylovar import check_seed,check_purity,random_int

#handle the error below
#python | head == IOError: [Errno 32] Broken pipe
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL)

def check_folder(directory=None):
    if not os.path.isdir(directory):
        raise argparse.ArgumentTypeError("'{}' doesn't exist or isn't a folder.".format(directory))
    return directory

def check_file(f=None):
    if not os.path.isfile(f):
        raise argparse.ArgumentTypeError("'{}' doesn't exist or isn't a file.".format(f))
    return f

def check_depth(value=None):
    fvalue=float(value)
    if fvalue<0:
        raise argparse.ArgumentTypeError("{} is an invalid value for read depth.".format(value)+
            "It should be a non-negative float number.")
    return fvalue


def main(progname=None):
    t0 = time.time()
    prog = progname if progname else sys.argv[0]
    parser=argparse.ArgumentParser(
        description='A wrapper of simulating WGS reads from normal and tumor genome fasta',
        prog=prog)
    group1 = parser.add_argument_group('Input arguments')
    group1.add_argument('-n','--normal',type=check_folder,required=True,metavar='DIR',
        help='the directory of the normal fasta')
    group1.add_argument('-t','--tumor',type=check_folder,required=True,metavar='DIR',
        help='the directory of the tumor fasta')
    group1.add_argument('-m','--map',type=check_folder,required=True,metavar='DIR',
        help='the directory of map files, which contains the relationship between tip nodes and samples')
    default=None
    group1.add_argument('-s','--sectors',type=check_file,default=default,metavar='FILE',
        help='the file containing purity and depth profile of each tumor sector. \
              After this setting, -d/-p will be ignored [{}]'.format(default))
    group2 = parser.add_argument_group('Arguments for simulation')
    default=50
    group2.add_argument('-d','--tumor_depth',type=check_depth,default=default,metavar='FLOAT',
        help='the mean depth of tumor sample for ART to simulate WGS reads [{}]'.format(default))
    default=0
    group2.add_argument('-D','--normal_depth',type=check_depth,default=default,metavar='FLOAT',
        help='the mean depth of normal sample for ART to simulate WGS reads [{}]'.format(default))
    default=0.6
    group2.add_argument('-p','--purity',type=check_purity,default=default,metavar='FLOAT',
        help='the proportion of tumor cells in simulated tumor sample [{}]'.format(default))
    default=None
    group2.add_argument('--random_seed',type=check_seed,metavar='INT',
        help='the seed for random number generator [{}]'.format(default))
    default=150
    group2.add_argument('--rlen',type=int,default=default,metavar='INT',
        help="the length of reads to simulate [{}]".format(default))
    default='art_illumina --noALN --quiet --paired --mflen 500 --sdev 20'
    group2.add_argument('--art',type=str,default=default,metavar='STR',
        help="the parameters for ART program ['{}']".format(default))
    default=1
    group2.add_argument('--cores',type=int,default=default,metavar='INT',
        help='number of cores used to run the program [{}]'.format(default))
    group2.add_argument('--separate',action="store_true",
        help="keep each tip node's WGS reads file separately")
    group2.add_argument('--single',action="store_true",
        help="single cell mode. "+\
        "After this setting,  -p will be ignored and the value of --tumor_depth is the depth of each tumor cell "+\
        "(not the total depth of tumor sample anymore).")
    group3 = parser.add_argument_group('Output arguments')
    default='art_reads'
    group3.add_argument('-o','--output',type=str,default=default,metavar='DIR',
        help='output directory [{}]'.format(default))
    default='fa2wgs.log'
    group3.add_argument('-g','--log',type=str,default=default,metavar='FILE',
        help='the log file to save the settings of each command [{}]'.format(default))
    args=parser.parse_args()

#always compress the simulated fastq files
    compress=True

#logging and random seed setting
    logging.basicConfig(filename=args.log,
        filemode='w',format='[%(asctime)s] %(levelname)s: %(message)s',
        datefmt='%m-%d %H:%M:%S',level='INFO')
    argv_copy=sys.argv[:]
    if '--art' in argv_copy:
        art_index=argv_copy.index('--art')
        argv_copy[art_index+1]="'{}'".format(argv_copy[art_index+1])
    argv_copy.insert(1,'fa2wgs')
    logging.info(' Command: %s',' '.join(argv_copy))
    if args.random_seed==None:
        seed=random_int()
    else:
        seed=args.random_seed
    logging.info(' Random seed: %s',seed)
    numpy.random.seed(seed)

#construct the sectors dictionary to store the meta information of all tumor sectors
    sectors={}
    if args.sectors!=None:
        sectors=read_sectors_file(f=args.sectors)
        for sector in sectors:
            mapfile=os.path.join(args.map,'{}.tipnode.map'.format(sector))
            assert os.path.isfile(mapfile),\
                "Couldn't find the map file ({}.tipnode.map) for sector '{}' ".format(sector,sector)+\
                "under the map directory ({}).".format(os.path.abspath(args.map))
    else:
        mapfiles=glob.glob(os.path.join(args.map,'*.tipnode.map'))
        infered_sectors=['.'.join(os.path.basename(x).split('.')[:-2]) for x in mapfiles]
        for sector in infered_sectors:
            sectors[sector]={'purity':args.purity,'depth':args.tumor_depth}
    for sector in sectors:
        mapfile=os.path.join(args.map,'{}.tipnode.map'.format(sector))
        sectors[sector]['composition']=tipnode_leaves_counting(f=mapfile)

#exit the program if you do NOT want to simulate any reads for normal and tumor samples
    if args.normal_depth==0:
        for sector in sectors:
            if sectors[sector]['depth']!=0:
                break
        else:
            sys.exit('Do nothing as the depth for each sample is 0!')

#single cell mode or bulk tumor mode
    if args.single:
        for sector in sectors:
            for tipnode,leaves_n in sectors[sector]['composition'].items():
                assert leaves_n==1,\
                    'In single mode, each tip node should represent only one cell.\n'+\
                    'But {} leaves are found underneath tipnode {} in one of your map files!'.format(leaves_n,tipnode)

#create index file (.fai) for each fasta
    pool=multiprocessing.Pool(processes=args.cores)
    tipnodes=set()
    for sector in sectors:
        tipnodes=tipnodes.union(set(sectors[sector]['composition'].keys()))
    results=[]
    for parental in 0,1:
        fasta=os.path.join(args.normal,'normal.parental_{}.fa'.format(parental))
        assert os.path.isfile(fasta),\
            "Couldn't find {} under the normal directory: {}".format(fasta,args.normal)
        results.append(pool.apply_async(build_fai,args=(fasta,)))
        for tipnode in tipnodes:
            fasta=os.path.join(args.tumor,'{}.parental_{}.fa'.format(tipnode,parental))
            assert os.path.isfile(fasta),\
                "Couldn't find {} under the tumor directory: {}".format(fasta,args.tumor)
            results.append(pool.apply_async(build_fai,args=(fasta,)))
    pool.close()
    pool.join()
    for result in results:
        result.get()

#create output folders
    if os.path.exists(args.output):
        if os.path.isdir(args.output):
            pass
        else:
            raise OutputExistsError("A FILE in the name of '{}' exists.\nDelete it or try another name as output folder.".format(args.output))
    else:
        os.mkdir(args.output,mode=0o755)
    normal_dir=os.path.join(args.output,'normal')
    if args.normal_depth>0:
        try:
            os.mkdir(normal_dir,mode=0o755)
        except FileExistsError as e:
            raise OutputExistsError("'{}' exists already! \nCan not use it as the output folder of normal WGS reads.".format(normal_dir)+
                '\nDelete it or use another folder as output folder.') from e

#collect simulation parameters first
    params_matrix=[]
    total_sim_bases=0
    art_params=args.art

#collect genome size for each genome
    normal_gsize=0
    for parental in 0,1:
        normal_gsize+=genomesize(fasta=os.path.join(args.normal,'normal.parental_{}.fa'.format(parental)))
    tipnode_gsize={}
    for tipnode in tipnodes:
#The value of tipnode_gsize[tipnode] is a list of three elements:
#0)genomesize of parental 0
#1)genomesize of parental 1
#2)the sum of parental 0 and 1
        tipnode_gsize[tipnode]=[]
        for parental in 0,1:
            tipnode_gsize[tipnode].append(genomesize(fasta=os.path.join(args.tumor,'{}.parental_{}.fa'.format(tipnode,parental))))
        tipnode_gsize[tipnode].append(tipnode_gsize[tipnode][0]+tipnode_gsize[tipnode][1])

#simulation for normal sample
    if args.normal_depth>0:
        for parental in 0,1:
            prefix=os.path.join(normal_dir,'normal.parental_{}.'.format(parental))
            fcov=args.normal_depth/2
            ref=os.path.join(args.normal,'normal.parental_{}.fa'.format(parental))
            sim_cfg={
                'gsize':normal_gsize/2,
                'base_cmd':art_params,
                'rlen':args.rlen,
                'fcov':fcov,
                'in':ref,
                'out':prefix,
                'id':'nm_prt{}'.format(parental)}
            params_matrix.append(sim_cfg)
            total_sim_bases+=sim_cfg['gsize']*sim_cfg['fcov']

#simulation for tumor sample
    for sector in sorted(sectors.keys()):
        if sectors[sector]['depth']>0:
#compute coverage and run ART
            sector_dir=os.path.join(args.output,sector)
            try:
                os.mkdir(sector_dir,mode=0o755)
            except FileExistsError as e:
                raise OutputExistsError("'{}' exists already! \nCan not use it as the output folder of tumor WGS reads.".format(sector_dir)+
                    '\nDelete it or use another folder as output folder.') from e

            tipnode_leaves=sectors[sector]['composition']
            sector_sim_bases=normal_gsize/2*sectors[sector]['depth']
            tumor_cells=sum(tipnode_leaves.values())
            total_cells=tumor_cells/sectors[sector]['purity']
            normal_cells=total_cells-tumor_cells
            normal_dna=normal_gsize*normal_cells
            tumor_dna=0
            for tipnode,leaves_n in tipnode_leaves.items():
                tumor_dna+=tipnode_gsize[tipnode][2]*leaves_n
            mean_depth_per_base=sector_sim_bases/(normal_dna+tumor_dna)

#two normal cell haplotypes
            if not args.single:
                for parental in 0,1:
                    prefix=os.path.join(sector_dir,'normal.parental_{}.'.format(parental))
                    fcov=normal_cells*mean_depth_per_base
                    ref=os.path.join(args.normal,'normal.parental_{}.fa'.format(parental))
                    sim_cfg={
                        'gsize':normal_gsize/2,
                        'base_cmd':art_params,
                        'rlen':args.rlen,
                        'fcov':fcov,
                        'in':ref,
                        'out':prefix,
                        'id':'nm_prt{}'.format(parental)}
                    params_matrix.append(sim_cfg)
                    total_sim_bases+=sim_cfg['gsize']*sim_cfg['fcov']

#tumor cells haplotypes
            for tipnode in sorted(tipnode_leaves.keys()):
                fcov=None
                if args.single:
                    fcov=sector_sim_bases/tipnode_gsize[tipnode][2]
                else:
                    fcov=tipnode_leaves[tipnode]*mean_depth_per_base
                for parental in 0,1:
                    ref=os.path.join(args.tumor,'{}.parental_{}.fa'.format(tipnode,parental))
                    prefix=os.path.join(sector_dir,'{}.parental_{}.'.format(tipnode,parental))
                    sim_cfg={
                        'gsize':tipnode_gsize[tipnode][parental],
                        'base_cmd':art_params,
                        'rlen':args.rlen,
                        'fcov':fcov,
                        'in':ref,
                        'out':prefix,
                        'id':'{}_prt{}'.format(tipnode,parental)}
                    params_matrix.append(sim_cfg)
                    total_sim_bases+=sim_cfg['gsize']*sim_cfg['fcov']

#generate fastq and compress them parallelly
#every thread will generate at most 2 percent of the total data you want to simulate
#In order to let users replicate the results (with same random seed) even using different number of cores,
#I use the fixed size of block to parallelize the program.
    assert total_sim_bases>0,'The genome sizes of all cells in the sample is 0!'
    sizeBlock=total_sim_bases*0.02
    final_params_matrix=[]
    for cfg in params_matrix:
        n=math.ceil(cfg['gsize']*cfg['fcov']/sizeBlock)
        if n==0:
            continue
        cfg['fcov']=round(cfg['fcov']/n,6)
        for i in range(n):
            final_params_matrix.append(cfg.copy())
            final_params_matrix[-1]['out']=cfg['out']+'{:03d}.'.format(i)
            final_params_matrix[-1]['id']=cfg['id']+'_{:03d}-'.format(i)
            final_params_matrix[-1]['rndSeed']=str(random_int())
    pool=multiprocessing.Pool(processes=args.cores)
    results=[]
    for x in final_params_matrix:
        results.append(pool.apply_async(generate_fq,args=(x,compress)))
    pool.close()
    pool.join()
    for result in results:
        result.get()

#merge small fastq files into one fastq for normal/tumor sample
    sample_fq_files=[]
    suffixes=['fq','1.fq','2.fq']
    if compress:
        suffixes=[x+'.gz' for x in suffixes]
    if args.normal_depth>0:
        for suffix in suffixes:
            prefix=os.path.join(normal_dir,'normal.parental_[01].[0-9][0-9][0-9].')
            source=glob.glob(prefix+suffix)
            if len(source):
                target=os.path.join(normal_dir,'normal.{}'.format(suffix))
                source.sort()
                sample_fq_files.append([target,source])
    for sector in sorted(sectors.keys()):
        if sectors[sector]['depth']>0:
            sector_dir=os.path.join(args.output,sector)
            tipnode_leaves=sectors[sector]['composition']
            for suffix in suffixes:
                if args.single or args.separate:
                    for tipnode in ['normal']+sorted(tipnode_leaves.keys()):
                        prefix=os.path.join(sector_dir,'{}.parental_[01].[0-9][0-9][0-9].'.format(tipnode))
                        source=glob.glob(prefix+suffix)
                        if len(source):
                            target=os.path.join(sector_dir,'{}.{}'.format(tipnode,suffix))
                            source.sort()
                            sample_fq_files.append([target,source])
                else:
                    prefix=os.path.join(sector_dir,'*.parental_[01].[0-9][0-9][0-9].')
                    source=glob.glob(prefix+suffix)
                    if len(source):
                        target=os.path.join(sector_dir,'{}.{}'.format(sector,suffix))
                        source.sort()
                        sample_fq_files.append([target,source])
    pool=multiprocessing.Pool(processes=args.cores)
    results=[]
    for x in sample_fq_files:
        results.append(pool.apply_async(merge_fq,args=x))
    pool.close()
    pool.join()
    for result in results:
        result.get()

    t1 = time.time()
    print ("Total time running {}: {} seconds".format
       (prog, str(t1-t0)))

def build_fai(fasta=None):
    '''
    In order to handle exceptions in child process--pyfaidx.Faidx,
    I must use the mothod result.get().
    But just using pyfaidx.Faidx in apply_async will induce error as there is no return value of pyfaidx.Faidx.
    So I just wrapper the function here and add a string as the return value.
    '''
    pyfaidx.Faidx(fasta)
    return 'Built index for {}'.format(fasta)

def merge_fq(target=None,source=None, remove=True):
    '''
    After generating short reads in multiprocessing mode,
    there will be multiple fq files for each genome.
    I will merge them into one file for each genome.
    '''
    assert not os.path.isfile(target),"'{}' exists already!".format(target)
    with open(target,'wb') as outfile:
        for f in source:
            with open(f,'rb') as infile:
                shutil.copyfileobj(infile, outfile)
            if remove:
                os.remove(f)

def generate_fq(params=None,compress=False):
    '''
    run art command to generate the fastq file, and call compress_fq if required.
    '''
    cmd_params=params['base_cmd'].split()+['--len',str(params['rlen']),
                                           '--fcov',str(params['fcov']),
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
            with open(fq, 'rb') as infile:
                with gzip.open(fq+'.gz', 'wb') as outfile:
                    shutil.copyfileobj(infile, outfile)
            os.remove(fq)

def read_sectors_file(f=None):
    '''
    Read the sectors file and build a dictionary.
    The sectors file should have 3 columns:
    1: sector id
    2: purity of the sector
    3: simulation depth of the sector
    The built dictionay is with structure:
    {sector1:{'purity':purity1,'depth':depth1},
     sector2:{'purity':purity2,'depth':depth2}, ...}
    '''
    sectors={}
    with open(f,'r') as input:
        header=next(input)
        if not header.startswith('#'):
            raise SectorsFileError("The first line should be a header line that starts with '#'!")
        header=header.lstrip('#')
        header=header.rstrip()
        header=header.split()
        if header!=['sector','purity','depth']:
            raise SectorsFileError('The format of your sectors file is not right!')
        for line in input:
            line=line.rstrip()
            sector,purity,depth=line.split()[:3]
            if sector not in sectors:
                sectors[sector]={}
            else:
                raise SectorsFileError('Found two records about sector {} in your --sectors file.'.format(sector))
            sectors[sector]['purity']=float(purity)
            sectors[sector]['depth']=float(depth)
            if not 0<sectors[sector]['purity']<=1:
                raise SectorsFileError("{} is an invalid value for sector purity.".format(sectors[sector]['purity'])+
                    "It should be a float number in the range of (0,1].")
            if sectors[sector]['depth']<0:
                raise SectorsFileError("{} is an invalid value for read depth.".format(sectors[sector]['depth'])+
                    "It should be a non-negative float number.")
    return sectors

def tipnode_leaves_counting(f=None):
    '''
    Return a dictionay with structure:
    {tipnode1:leaves_count1,tipnode2:leaves_count2,...}
    '''
    tipnode_leaves={}
    with open(f,'r') as input:
        for line in input:
            if not line.startswith('#'):
                line=line.rstrip()
                tipnode,leaves_n=line.split()[:2]
                tipnode_leaves[tipnode]=int(leaves_n)
    return tipnode_leaves

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

class SectorsFileError(Exception):
    pass
