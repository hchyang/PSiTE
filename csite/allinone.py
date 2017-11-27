#!/usr/bin/env python3

#########################################################################
# Author: Hechuan
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
from csite.phylovar import check_prune,check_proportion,check_seed,random_int,check_config_file
from csite.vcf2fa import check_sex

#handle the error below
#python | head == IOError: [Errno 32] Broken pipe 
from signal import signal, SIGPIPE, SIG_DFL 
signal(SIGPIPE,SIG_DFL) 

def main(progname=None):
    parse=argparse.ArgumentParser(
        description='a wrapper of simulating short reads from genome with germline and somatic variants',
        prog=progname if progname else sys.argv[0])
    parse.add_argument('-r','--reference',type=str,required=True,
        help='a fasta file of the reference genome')
    parse.add_argument('-v','--vcf',type=str,required=True,
        help='a vcf file contains germline variants')
    parse.add_argument('-t','--tree',type=str,required=True,
        help='a newick file contains ONE tree')
    parse.add_argument('-c','--config',type=str,required=True,
        help='a yaml file contains configure for somatic variant simulation')
    parse.add_argument('-o','--output',type=str,required=True,
        help='output folder')
    parse.add_argument('-a','--autosomes',type=str,required=True,
        help='autosomes of the genome (seperated by comma)')
    default=None
    parse.add_argument('-s','--sex_chr',type=check_sex,default=default,
        help='sex chromosomes of the genome (seperated by comma) [{}]'.format(default))
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
    parse.add_argument('-d','--depth',type=float,default=default,
        help='the mean depth of tumor for ART to simulate short reads [{}]'.format(default))
    default=50
    parse.add_argument('-D','--normal_depth',type=float,default=default,
        help='the mean depth of normal for ART to simulate short reads [{}]'.format(default))
    default=0.5
    parse.add_argument('-p','--purity',type=float,default=default,
        help='the proportion of tumor cells in simulated sample [{}]'.format(default))
    default='art_illumina --noALN --quiet --paired --len 100 --mflen 500 --sdev 20'
    parse.add_argument('--art',type=str,default=default,
        help='the parameters for ART program [{}]'.format(default))
    default=None
    parse.add_argument('--random_seed',type=check_seed,
        help='the seed for random number generator [{}]'.format(default))
    default='allinone.log'
    parse.add_argument('--log',type=str,default=default,
        help='the log file to save the settings of each command [{}]'.format(default))
    default=1
    parse.add_argument('--start',type=int,default=default,choices=[1,2,3,4],
        help='the serial number of the module from which to start [{}]'.format(default))
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
    if args.start==1:
        try:
            os.mkdir(outdir)
        except FileExistsError:
            exit('{} already exists. Try another directory to output! (-o/--output)'.format(outdir))
    else:
        assert os.path.isdir(outdir),'Can not start from step {}, '.format(args.start)+'because I can not find previous results directory {}.'.format(outdir)
    os.chdir(outdir)

###### logging and random seed setting
    logging.basicConfig(filename=args.log if args.start==1 else args.log+'.start'+str(args.start), 
        filemode='w',format='[%(asctime)s] %(levelname)s: %(message)s',
        datefmt='%m-%d %H:%M:%S',level='INFO')
    logging.info(' Command: %s',' '.join(sys.argv[0:1]+['allinone']+sys.argv[1:]))
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

#vcf2fa
    if args.start<2:
        cmd_params=[sys.argv[0],'vcf2fa',
                    '--vcf',vcf,
                    '--reference',reference,
                    '--output',normal_fa,
                    '--autosomes',args.autosomes]
        if args.sex_chr:
            cmd_params.extend(['--sex_chr',','.join(args.sex_chr)])
        logging.info(' Command: %s',' '.join(cmd_params))
        subprocess.run(args=cmd_params,check=True)

#phylovar
#I place random_int() here as I do not want to skip it in all situation.
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
                    '--random_seed',str(random_n),
                    '--chain',tumor_chain,
                    '--depth',str(args.depth),
                    '--purity',str(args.purity)]
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
                    '--reference','{dir}/normal_parental_0.fa,{dir}/normal_parental_1.fa'.format(dir=normal_fa),
                    '--output',tumor_fa]
        logging.info(' Command: %s',' '.join(cmd_params))
        subprocess.run(args=cmd_params,check=True)

    if os.path.isdir(art_reads):
        shutil.rmtree(art_reads)
    elif os.path.isfile(art_reads):
        os.remove(art_reads)
    cmd_params=[sys.argv[0],'fa2ngs',
                '--normal',normal_fa,
                '--tumor',tumor_fa,
                '--chain',tumor_chain,
                '--depth',str(args.depth),
                '--normal_depth',str(args.normal_depth),
                '--purity',str(args.purity),
                '--random_seed',str(random_int()),
                '--output',art_reads,
                '--art',"{}".format(args.art)]
    logging.info(' Command: %s',' '.join(cmd_params[:-1]+["'{}'".format(args.art)]))
    subprocess.run(args=cmd_params,check=True)

