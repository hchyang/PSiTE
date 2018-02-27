#!/usr/bin/env python3

#########################################################################
# Author: Hechuan Yang
# Created Time: 2017-04-04 18:00:34
# File Name: allinone.py
# Description: 
#########################################################################

import sys
import os
import shutil
import argparse
import numpy
import yaml
import logging
import subprocess
import pyfaidx
from csite.vcf2fa import check_sex,check_vcf,check_autosomes
from csite.phylovar import check_prune,check_proportion,check_seed,check_purity,random_int,check_config_file
from csite.fa2wgs import check_depth

#handle the error below
#python | head == IOError: [Errno 32] Broken pipe 
from signal import signal, SIGPIPE, SIG_DFL 
signal(SIGPIPE,SIG_DFL) 

def main(progname=None):
    parse=argparse.ArgumentParser(
        description='an all-in-one wrapper for NGS reads simulation for tumor samples',
        prog=progname if progname else sys.argv[0])
    parse.add_argument('-r','--reference',type=str,required=True,metavar='FILE',
        help='a fasta file of the reference genome')
    parse.add_argument('-v','--vcf',type=check_vcf,required=True,metavar='FILE',
        help='a vcf file contains germline variants')
    parse.add_argument('-t','--tree',type=str,required=True,metavar='FILE',
        help='a newick file contains ONE tree')
    parse.add_argument('-c','--config',type=str,required=True,metavar='FILE',
        help='a YAML file which contains the configuration of somatic variant simulation')
    parse.add_argument('-o','--output',type=str,required=True,metavar='DIR',
        help='output directory')
    parse.add_argument('-a','--autosomes',type=check_autosomes,required=True,metavar='STR',
        help='autosomes of the genome (e.g. 1,2,3,4,5 or 1..4,5)')
    default=None
    parse.add_argument('-s','--sex_chr',type=check_sex,default=default,metavar='STR',
        help='sex chromosomes of the genome (seperated by comma) [{}]'.format(default))
    default=0
    parse.add_argument('-x','--prune',type=check_prune,default=default,metavar='INT',
        help='trim all the children of the nodes with equal or less than this number of leaves [{}]'.format(default))
    default=0.0
    parse.add_argument('-X','--prune_proportion',type=check_proportion,default=default,metavar='FLOAT',
        help='trim all the children of the nodes with equal or less than this proportion of total leaves [{}]'.format(default))
    default=50
    parse.add_argument('-d','--depth',type=check_depth,default=default,metavar='FLOAT',
        help='the mean depth of tumor sample for ART to simulate NGS reads [{}]'.format(default))
    default=0
    parse.add_argument('-D','--normal_depth',type=check_depth,default=default,metavar='FLOAT',
        help='the mean depth of normal sample for ART to simulate NGS reads [{}]'.format(default))
    default=0.8
    parse.add_argument('-p','--purity',type=check_purity,default=default,metavar='FLOAT',
        help='the proportion of tumor cells in simulated tumor sample [{}]'.format(default))
    default=None
    parse.add_argument('--trunk_vars',type=str,default=default,metavar='FILE',
        help='a file containing truncal variants predefined by user [{}]'.format(default))
    default=0
    parse.add_argument('--trunk_length',type=float,default=default,metavar='FLOAT',
        help='the length of the trunk [{}]'.format(default))
    default='art_illumina --noALN --quiet --paired --len 100 --mflen 500 --sdev 20'
    parse.add_argument('--art',type=str,default=default,metavar='STR',
        help='the parameters for ART program [{}]'.format(default))
    default=None
    parse.add_argument('--random_seed',type=check_seed,default=default,metavar='INT',
        help='the seed for random number generator [{}]'.format(default))
    default='allinone.log'
    parse.add_argument('--log',type=str,default=default,metavar='FILE',
        help='the log file to save the settings of each command [{}]'.format(default))
    default=1
    parse.add_argument('--start',type=int,default=default,choices=[1,2,3,4],
        help='the serial number of the module from which to start [{}]'.format(default))
    default=1
    parse.add_argument('--cores',type=int,default=default,metavar='INT',
        help='number of cores used to run the program [{}]'.format(default))
    parse.add_argument('--compress',action="store_true",
        help='compress the generated fastq files using gzip')
    parse.add_argument('--seperate',action="store_true",
        help="keep each tip node's NGS reads file seperately")
    parse.add_argument('--single',action="store_true",
        help="single cell mode. "+\
        "After this setting, the value of --depth is the depth of each tumor cell "+\
        "(not the total depth of tumor sample anymore).")
    args=parse.parse_args()
    if (args.prune or args.prune_proportion) and args.single:
        raise argparse.ArgumentTypeError("Can not prune the tree in single cell mode!")
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
    if args.start==1:
        try:
            os.mkdir(outdir,mode=0o755)
        except FileExistsError as e:
            raise OutputExistsError("'{}' already exists. Try another directory to output! (-o/--output)".format(outdir)) from e
    else:
        assert os.path.isdir(outdir),"Couldn't start from step {}, ".format(args.start)+\
            "because I can not find the directory of previous results: '{}'.".format(outdir)
    os.chdir(outdir)

###### logging and random seed setting
    logging.basicConfig(filename=args.log if args.start==1 else args.log+'.start'+str(args.start), 
        filemode='w',format='[%(asctime)s] %(levelname)s: %(message)s',
        datefmt='%m-%d %H:%M:%S',level='INFO')
    argv_copy=sys.argv[:]
    try:
        art_index=argv_copy.index('--art')
        argv_copy[art_index+1]="'{}'".format(argv_copy[art_index+1])
    except ValueError:
        pass
    argv_copy.insert(1,'allinone')
    logging.info(' Command: %s',' '.join(argv_copy))
    if args.random_seed==None:
        seed=random_int()
    else:
        seed=args.random_seed
    logging.info(' Random seed: %s',seed)
    numpy.random.seed(seed)

#subfolders
    normal_fa='normal_fa'
    tumor_fa='tumor_fa'
    tumor_chain='tumor_chain'
    art_reads='art_reads'
#map file
    map_file='tipnode_samples.map'

#vcf2fa
    if args.start<2:
        cmd_params=[sys.argv[0],'vcf2fa',
                    '--vcf',vcf,
                    '--reference',reference,
                    '--output',normal_fa,
                    '--autosomes',args.autosomes]
        if args.sex_chr:
            cmd_params.extend(['--sex_chr',args.sex_chr])
        logging.info(' Command: %s',' '.join(cmd_params))
        subprocess.run(args=cmd_params,check=True)

#phylovar
#I place random_int() here as I do not want to skip it in any situation.
#Without this, you can not replicate the result with different --start setting.
    random_n=random_int()
    if args.start<3:
        if os.path.isdir(tumor_chain):
            shutil.rmtree(tumor_chain)
        elif os.path.isfile(tumor_chain):
            os.remove(tumor_chain)
        cmd_params=[sys.argv[0],'phylovar',
                    '--tree',tree,
                    '--config',config,
                    '--purity',str(args.purity),
                    '--random_seed',str(random_n),
                    '--map',map_file,
                    '--chain',tumor_chain]
        if args.sex_chr:
            cmd_params.extend(['--sex_chr',args.sex_chr])
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
    if args.start<4:
        if os.path.isdir(tumor_fa):
            shutil.rmtree(tumor_fa)
        elif os.path.isfile(tumor_fa):
            os.remove(tumor_fa)

        cmd_params=[sys.argv[0],'chain2fa',
                    '--chain',tumor_chain,
                    '--normal','{dir}/normal.parental_0.fa,{dir}/normal.parental_1.fa'.format(dir=normal_fa),
                    '--cores',str(args.cores),
                    '--output',tumor_fa]
        logging.info(' Command: %s',' '.join(cmd_params))
        subprocess.run(args=cmd_params,check=True)

#fa2wgs
    if os.path.isdir(art_reads):
        shutil.rmtree(art_reads)
    elif os.path.isfile(art_reads):
        os.remove(art_reads)
    cmd_params=[sys.argv[0],'fa2wgs',
                '--normal',normal_fa,
                '--tumor',tumor_fa,
                '--map',map_file,
                '--depth',str(args.depth),
                '--normal_depth',str(args.normal_depth),
                '--purity',str(args.purity),
                '--output',art_reads,
                '--random_seed',str(random_int()),
                '--cores',str(args.cores),
                '--art','{}'.format(args.art)]
    if args.compress:
        cmd_params.extend(['--compress'])
    if args.single:
        cmd_params.extend(['--single'])
    cmd_params_copy=cmd_params[:]
    art_index=cmd_params_copy.index('--art')
    cmd_params_copy[art_index+1]="'{}'".format(cmd_params_copy[art_index+1])
    logging.info(' Command: %s',' '.join(cmd_params_copy))
    subprocess.run(args=cmd_params,check=True)

class OutputExistsError(Exception):
    pass
