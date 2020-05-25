#!/usr/bin/env python3

#########################################################################
# Author: Hechuan Yang
# Created Time: 2017-04-04 18:00:34
# File Name: phylovar.py
# Description:
#########################################################################

import sys
import os
import re
import argparse
import numpy
import logging
import copy
import yaml
import time
import psite.trunk_vars
import psite.tree
from psite.vcf2fa import check_sex

#handle the error below
#python | head == IOError: [Errno 32] Broken pipe
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL)

WHOLET='tumor'
#I defined those two variables as global variables. As they will be used in function
#random_int and check_config_file, which are also used in allinone.py.
LARGEST=2**31
CFG_PARAMS={'snv_rate':float,
            'cnv_rate':float,
            'trunk_snv_rate':float,
            'trunk_cnv_rate':float,
            'del_prob':float,
            'tandem_prob':float,
            'cnv_length_beta':int,
            'cnv_length_max':int,
            'copy_parameter':float,
            'copy_max':int,
            'parental':str,
            'tstv':float,
            'length':int,
            }

def random_int():
    '''
    The random seed for numpy must be convertible to 32 bit unsigned integers.
    But in jave (used in fa2wes) the range is (-2**31,2**31-1).
    Let's use this to generate a integers can be used by both.
    '''
    return numpy.random.randint(LARGEST)

def check_seed(value=None):
    ivalue=int(value)
#2**32: Must be convertible to 32 bit unsigned integers.
#in jave the range is (-2**31,2**31-1)
    if not 0<=ivalue<LARGEST:
        raise argparse.ArgumentTypeError("{} is an invalid value for --random_seed. ".format(value)+
            "It should be an integer between 0 and {}.".format(LARGEST-1))
    return ivalue

def check_prune(value=None):
    fvalue=float(value)
    if not 0<=fvalue<=1:
        raise argparse.ArgumentTypeError("{} is an invalid value for --prune. ".format(value)+
            "It should be an float number in the range of [0,1].")
    return fvalue

def check_del_prob(value=None):
    fvalue=float(value)
    if not 0<=fvalue<=1:
        raise argparse.ArgumentTypeError("{} is an invalid value for --del_prob. ".format(value)+
            "It should be an float number in the range of [0,1].")
    return fvalue

def check_tandem_prob(value=None):
    fvalue=float(value)
    if not 0<=fvalue<=1:
        raise argparse.ArgumentTypeError("{} is an invalid value for --tandem_prob. ".format(value)+
            "It should be an float number in the range of [0,1].")
    return fvalue

def check_tstv(value=None):
    fvalue=float(value)
    if fvalue<=0:
        raise argparse.ArgumentTypeError("{} is an invalid value for transition/transversion ratio. ".format(value)+
            "It should be an float larger than 0.")
    return fvalue

def check_purity(value=None):
    fvalue=float(value)
    if not 0<fvalue<=1:
        raise argparse.ArgumentTypeError("{} is an invalid value for --purity. ".format(value)+
            "It should be a float number in the range of (0,1].")
    return fvalue

def check_depth(value=None):
    fvalue=float(value)
    if fvalue<=0:
        raise argparse.ArgumentTypeError("{} is an invalid value for --depth. ".format(value)+
            "It should be a positive float number.")
    return fvalue

def check_folder(directory=None):
    good_charactors=re.compile('^[0-9a-zA-Z/_\-.]+$')
    if not good_charactors.match(directory):
        raise argparse.ArgumentTypeError("'{}' is an invalid string for folder name. ".format(directory)+
            "Please only use number, alphabet and ._/- in the directory name.")
    if os.path.exists(directory):
        raise argparse.ArgumentTypeError("'{}' exists already. Delete it or use another name instead.".format(directory))
    return directory

def check_cnv_length_cfg(chroms=None,cnv_length_beta=None,cnv_length_max=None,chr_length=None):
    if cnv_length_max>chr_length:
         raise argparse.ArgumentTypeError("{}: The value of cnv_length_max ".format(chroms)+
            "({}) should NOT be larger than the length of chromosome ({}).".format(cnv_length_max,chr_length))
    if cnv_length_max<cnv_length_beta:
         raise argparse.ArgumentTypeError("{}: The value of cnv_length_max ".format(chroms)+
            "({}) should be larger than the cnv_length_beta ({}).".format(cnv_length_max,cnv_length_beta))

def cn_dist(copy_max=None,copy_parameter=None):
    '''
    When an amplification event happens, it will randomly pick a copy number according to a geometric like distribution,
    Pr(n+1)=p*Pr(n).
    This function returns the configure of the distribution of CNVs' copy number.
    '''
    cn_dist_cfg={}
    cn_dist_cfg['copy']=[x for x in range(1,copy_max+1)]
    w=1/sum([copy_parameter**x for x in range(copy_max)])
    scale=sum([w*copy_parameter**x for x in range(copy_max)])
    cn_dist_cfg['prob']=[(w*copy_parameter**x)/scale for x in range(copy_max)]
    return cn_dist_cfg

def read_cnvl_dist(cfg_f=None):
    with open(cfg_f) as input:
        header=next(input)
        if not header.startswith('#'):
            raise CnvDistFileError("The first line should be a header line that starts with '#'!")
        header=header.lstrip('#')
        header=header.rstrip()
        header=header.split()
        if header!=['low','high','prob']:
            raise CnvDistFileError('The format of your CNV distribution file is not right!')
        i=0
        sum_prob=0
        cnvl_dist={'index':[],'bins':[],'prob':[]}
        for line in input:
            line=line.rstrip()
            cols=line.split()
            if len(cols)!=3:
                raise CnvDistFileError('The format of your CNV distribution file is not right!')
            low=int(cols[0])
            high=int(cols[1])
            prob=float(cols[2])
            if not (low>0 and high>0 and high>low and 0<prob<=1):
                raise CnvDistFileError('Check the record below from your CNV distribution file:\n'
                    +'{}\n'.format(line))
            sum_prob+=prob
            cnvl_dist['index'].append(i)
            cnvl_dist['bins'].append([low,high])
            cnvl_dist['prob'].append(prob)
            i+=1
        if sum_prob!=1:
            raise CnvDistFileError('The sum of all the probability in your CNV distribution file is not 1!')
        return cnvl_dist

def tstv_dist(tstv=None):
    '''
    For a SNV:
    ts+tv=1 in which tv=tv1+tv2 and ts/tv=tstv
    This function returns the configure of the distribution of different form of SNV mutation.
    '''
    tstv_dist_cfg={}
    tv=1/(1+tstv)
    ts=1-tv
    tv1=tv/2
    tv2=1-ts-tv1
    tstv_dist_cfg['form']=[0,1,2]
    tstv_dist_cfg['prob']=[ts,tv1,tv2]
    return tstv_dist_cfg

def check_config_file(config=None):
    '''
    Validate the settings in config file.
    1. There must be only two section in config file: genome and chromosomes.
    2. Every parameter in CFG_PARAMS must have a pair of key:value in genome section.
    3. Every chromosome at least has the key 'length'.
    4. The total length of all chromosomes must equal the length of genome.
    5. If the (trunk) SNV rate of each chromosome has been specified, the sum of them should be equal the (trunk) SNV rate of genome.
    6. If the (trunk) SNV rate of some chromosomes (not all) have been specified, the sum of them should be <= the (trunk) SNV rate of genome.
    7. The CNV rate should satisfy the same criteria as SNV rate.
    '''
    if not (isinstance(config,dict) and set(config)==set(['genome','chromosomes'])):
        raise ConfigFileError('Check your config file. The format is not correct.')
    if not isinstance(config['genome'],dict):
        raise ConfigFileError('Check your config file. The format in genome section is not correct.')
    lack=set(CFG_PARAMS)-set(config['genome'])
    if lack!=set():
        raise ConfigFileError("'{}' are required in the genome section of config file.".format(','.join([str(x) for x in lack])))
    over=set(config['genome'])-set(CFG_PARAMS)
    if over!=set():
        raise ConfigFileError("'{}' are not acceptable parameters in config file.".format(','.join([str(x) for x in over])))

    for parameter in CFG_PARAMS:
        assert isinstance(config['genome'][parameter],CFG_PARAMS[parameter]),\
            "'{}' in genome section of config file should be a {}!".format(parameter,CFG_PARAMS[parameter].__name__)

    if not isinstance(config['chromosomes'],list):
        raise ConfigFileError('Check your config file. The format in chromosomes section is not correct.')
    total_chroms_length=0
    total_snv_rate=0
    total_cnv_rate=0
    total_trunk_snv_rate=0
    total_trunk_cnv_rate=0
    any_chr_missing_snv_rate=False
    any_chr_missing_cnv_rate=False
    any_chr_missing_trunk_snv_rate=False
    any_chr_missing_trunk_cnv_rate=False
    for chroms in config['chromosomes']:
        if not isinstance(chroms,dict) or len(chroms)>1:
            raise ConfigFileError('Check your config file. The format in chromosomes section is not correct.')
        for chroms_n,chroms_cfg in chroms.items():
            assert isinstance(chroms_n,str),'The name of chromosome {} in config file should be a str!'.format(chroms_n)
            over=set(chroms_cfg)-set(CFG_PARAMS)
            if over!=set():
                raise ConfigFileError("'{}' are not acceptable parameters ".format(','.join([str(x) for x in over]))+
                    'in the section of chromosome {} in your config file!'.format(chroms))
            for parameter in chroms_cfg.keys():
                assert isinstance(chroms_cfg[parameter],CFG_PARAMS[parameter]),\
                    "'{}' of chromosome {} in config file should be a {}!".format(parameter,chroms_n,CFG_PARAMS[parameter].__name__)
            if 'length' not in chroms_cfg:
                raise ConfigFileError("Couldn't find the length of chromosome:{}.".format(chroms_n))
            total_chroms_length+=chroms_cfg['length']
            if 'snv_rate' in chroms_cfg:
                total_snv_rate+=chroms_cfg['snv_rate']
            else:
                any_chr_missing_snv_rate=True
            if 'cnv_rate' in chroms_cfg:
                total_cnv_rate+=chroms_cfg['cnv_rate']
            else:
                any_chr_missing_cnv_rate=True
            if 'trunk_snv_rate' in chroms_cfg:
                total_trunk_snv_rate+=chroms_cfg['trunk_snv_rate']
            else:
                any_chr_missing_trunk_snv_rate=True
            if 'trunk_cnv_rate' in chroms_cfg:
                total_trunk_cnv_rate+=chroms_cfg['trunk_cnv_rate']
            else:
                any_chr_missing_trunk_cnv_rate=True

    if config['genome']['length']!=total_chroms_length:
        raise ConfigFileError('In your config file, the length of genome is {},'.format(str(config['genome']['length']))+
            'But the total length of all chromosomes are {}'.format(str(total_chroms_length)))
    if any_chr_missing_snv_rate:
        if total_snv_rate>config['genome']['snv_rate']:
            raise ConfigFileError('Check your config file, the sum of snv_rate in chromosmomes section is larger than the rate in genome section!')
    else:
        if total_snv_rate!=config['genome']['snv_rate']:
            raise ConfigFileError('Check your config file, the sum of snv_rate in chromosmomes section is not equal the rate in genome section!')
    if any_chr_missing_cnv_rate:
        if total_cnv_rate>config['genome']['cnv_rate']:
            raise ConfigFileError('Check your config file, the sum of cnv_rate in chromosmomes section is larger than the rate in genome section!')
    else:
        if total_cnv_rate!=config['genome']['cnv_rate']:
            raise ConfigFileError('Check your config file, the sum of cnv_rate in chromosmomes section is not equal the rate in genome section!')
    if any_chr_missing_trunk_snv_rate:
        if total_trunk_snv_rate>config['genome']['trunk_snv_rate']:
            raise ConfigFileError('Check your config file, the sum of trunk_snv_rate in chromosmomes section is larger than the rate in genome section!')
    else:
        if total_trunk_snv_rate!=config['genome']['trunk_snv_rate']:
            raise ConfigFileError('Check your config file, the sum of trunk_snv_rate in chromosmomes section is not equal the rate in genome section!')
    if any_chr_missing_trunk_cnv_rate:
        if total_trunk_cnv_rate>config['genome']['trunk_cnv_rate']:
            raise ConfigFileError('Check your config file, the sum of trunk_cnv_rate in chromosmomes section is larger than the rate in genome section!')
    else:
        if total_trunk_cnv_rate!=config['genome']['trunk_cnv_rate']:
            raise ConfigFileError('Check your config file, the sum of trunk_cnv_rate in chromosmomes section is not equal the rate in genome section!')

class ConfigFileError(Exception):
    pass

class AffiliationFileError(Exception):
    pass

class CnvDistFileError(Exception):
    pass

def read_affiliation(affiliation_f=None):
    '''
    Check the format of affiliation file and dump the data into the sectors dictionary.
    There should be 3 columns in the affiliation file.
    1. sector id
    2. purity
    3. depth
    4. prune proportion of the sector
    5. tumor cells in the sector
    '''
    sectors={}
    with open(affiliation_f) as input:
        header=next(input)
        if not header.startswith('#'):
            raise AffiliationFileError('The format of your affiliation file is not right!')
        header=header.lstrip('#')
        header=header.rstrip()
        header=header.split()
        if header!=['sector','purity','depth','prune_p','cells']:
            raise AffiliationFileError('The format of your affiliation file is not right!')
        for line in input:
            line=line.rstrip()
            cols=line.split()
            if len(cols)!=5:
                raise AffiliationFileError("The format of your affiliation file is not right!")
            sector=cols[0]
            purity=float(cols[1])
            depth=cols[2]
            prune_p=float(cols[3])
            if sector==WHOLET:
                raise AffiliationFileError(
                    "Please do not use '{}' as a sector name. In PSiTE, I use it as the name of whole tumor sample.".format(WHOLET))
            if not 0<=purity<=1:
                raise AffiliationFileError(
                    "The purity {} for sector {} is not valid in your affiliation file.\n".format(purity,sector)+\
                    "It should be a float number in the range of [0,1]")
            if depth=='-':
                depth=None
            else:
                try:
                    depth=float(depth)
                except ValueError:
                    raise AffiliationFileError(
                        "'{}' is an invalide value for depth in your affiliation file.".format(depth))
            if not 0<=prune_p<=1:
                raise AffiliationFileError(
                    "The prune proportion {} for sector {} is not valid in your affiliation file.\n".format(prune_p,sector)+\
                    "It should be a float number in the range of [0,1]")
            cells=[]
            for i in cols[4].split(','):
                n=i.count('..')
                if n==0:
                    cells.append(i)
                elif n==1:
                    m=re.search('^(.*?)([0-9]+)\.\.(\\1)?([0-9]+)$',i)
                    if m:
                        prefix=m.group(1)
                        start=int(m.group(2))
                        end=int(m.group(4))
                        if start>=end:
                            raise AffiliationFileError("The string '{}' is not valid in your affiliation file.".format(i))
                        else:
                            cells.extend([prefix+str(x) for x in range(start,end+1)])
                    else:
                        raise AffiliationFileError("The string '{}' is not valid in your affiliation file.".format(i))
                else:
                    raise AffiliationFileError("The string '{}' is not valid in your affiliation file.".format(i))
            if sector in sectors:
                if prune_p!=sectors[sector]['prune_p']:
                    raise AffiliationFileError("Found two different prune proportions for sector {} in affiliation file:\n{} and {}"\
                        .format(prune_p,sectors[sector]['prune_p']))
                else:
                    sectors[sector]['members'].extend(cells)
            else:
                sectors[sector]={'purity':purity,'depth':depth,'prune_p':prune_p,'members':cells}
    for sector in sectors:
        sectors[sector]['members']=set(sectors[sector]['members'])
    return sectors

#use kernprof -l -v script.py to profile
# @profile
def main(progname=None):
    t0 = time.time()
    prog=progname if progname else sys.argv[0]
    parser=argparse.ArgumentParser(
        description='Simulate SNVs/CNVs on a phylogenetic tree in newick format',
        prog=prog)
    group1=parser.add_argument_group('Input arguments')
    group1.add_argument('-t','--tree',required=True,metavar='FILE',
        help='a file containing !!!ONE!!! tree in newick format')
    default=None
    group1.add_argument('--trunk_vars',type=str,default=default,metavar='FILE',
        help='a file containing truncal variants predefined by user [{}]'.format(default))
    default=None
    group1.add_argument('--config',type=str,default=default,metavar='FILE',
        help='a YAML file which contains the configuration of somatic variant simulation. '+
            '-n/-r/-R/-d/-l/-L/-c/-C/-p/--tstv/--length/--trunk_snv_rate/--trunk_cnv_rate will be ignored. [{}]'.format(default))
    default=None
    group1.add_argument('--affiliation',type=str,default=default,metavar='FILE',
        help='a file containing sector affiliation of the cells in the sample [{}]'.format(default))
    default=None
    group1.add_argument('--cnvl_dist',type=str,default=default,metavar='FILE',
        help="a file containing the distribution profile of CNVs' length [{}]".format(default))
    group2=parser.add_argument_group('Simulation arguments (can be set in config YAML)')
    default='1'
    group2.add_argument('-n','--name',type=str,default=default,metavar='STR',
        help='the name of the sequence to be simulated [{}]'.format(default))
    default=300
    group2.add_argument('-r','--snv_rate',type=float,default=default,metavar='FLOAT',
        help='the muation rate of SNVs [{}]'.format(default))
    default=3
    group2.add_argument('-R','--cnv_rate',type=float,default=default,metavar='FLOAT',
        help='the muation rate of CNVs [{}]'.format(default))
    default=None
    group2.add_argument('--trunk_snv_rate',type=float,default=default,metavar='FLOAT',
        help='the muation rate of SNVs on trunk. It will be the same as --snv_rate if not being specified [{}]'.format(default))
    default=None
    group2.add_argument('--trunk_cnv_rate',type=float,default=default,metavar='FLOAT',
        help='the muation rate of CNVs on trunk. It will be the same as --cnv_rate if not being specified [{}]'.format(default))
    default=0.5
    group2.add_argument('-d','--del_prob',type=float,default=default,metavar='FLOAT',
        help='the probability of being deletion for a CNV mutation [{}]'.format(default))
    default=1.0
    group2.add_argument('--tandem_prob',type=float,default=default,metavar='FLOAT',
        help='the probability of being tandem repeat for an amplification mutation [{}]'.format(default))
#https://en.wikipedia.org/wiki/Copy-number_variation
    default=20000000
    group2.add_argument('-l','--cnv_length_beta',type=int,default=default,metavar='INT',
        help='the mean of CNVs length [{}]'.format(default))
    default=40000000
    group2.add_argument('-L','--cnv_length_max',type=int,default=default,metavar='INT',
        help='the maximium of CNVs length [{}]'.format(default))
    default=0.5
    group2.add_argument('-c','--copy_parameter',type=float,default=default,metavar='FLOAT',
        help="the p parameter of CNVs' copy number distribution [{}]".format(default))
    default=5
    group2.add_argument('-C','--copy_max',type=int,default=default,metavar='INT',
        help='the maximium ADDITIONAL copy of a CNVs [{}]'.format(default))
    default='01'
    group2.add_argument('-p','--parental',type=str,default=default,metavar='STR',
        help='the parental to simulate [{}]'.format(default))
    default=2.0
    group2.add_argument('--tstv',type=check_tstv,default=default,metavar='FLOAT',
        help='the ratio of ts/tv of SNV [{}]'.format(default))
    default=100000000
    group2.add_argument('--length',type=int,default=default,metavar='INT',
        help='the length of the sequence to simulate [{}]'.format(default))
    group3=parser.add_argument_group('Other simulation arguments (can NOT be set in config YAML)')
    default=0.05
    group3.add_argument('-x','--prune',type=check_prune,default=default,metavar='FLOAT',
        help='trim all the children of the nodes with less than this proportion of total leaves [{}]'.format(default))
    default=None
    group3.add_argument('-s','--sex_chr',type=check_sex,default=default,metavar='STR',
        help='sex chromosomes of the genome (separated by comma) [{}]'.format(default))
    default=None
    group3.add_argument('--random_seed',type=check_seed,metavar='INT',
        help='the seed for random number generator (an integer between 0 and 2**31-1) [{}]'.format(default))
    default=0
    group3.add_argument('--trunk_length',type=float,default=default,metavar='FLOAT',
        help='the length of the trunk [{}]'.format(default))
    default=0.6
    group3.add_argument('--purity',type=check_purity,default=default,metavar='FLOAT',
        help='the proportion of tumor cells in simulated tumor sample [{}]'.format(default))
    default=None
    group3.add_argument('--depth',type=check_depth,default=default,metavar='FLOAT',
        help='the sequencing depth of the whole tumor sample for read count simulation [{}]'.format(default))
#Actually, the depth here and the depth in the affiliation file is the not the mean depth of the whole genome.
#It's impossible to get that without calculating the size of all the genomes in sample. 
#Here we set all the diploid part of the genome with this depth. And the depth of other part will be caculated 
#according this.
    default=150
    group3.add_argument('--rlen',type=int,default=default,metavar='INT',
        help='the read length for simulating the read count for each segment of the genome [{}]'.format(default))
    group4=parser.add_argument_group('Output arguments')
    default='phylovar_snvs'
    group4.add_argument('-S','--snv',type=str,default=default,metavar='DIR',
        help='the output directory to save SNVs files [{}]'.format(default))
    default='phylovar_cnvs'
    group4.add_argument('-V','--cnv',type=str,default=default,metavar='DIR',
        help='the output directory to save CNVs files [{}]'.format(default))
    default='phylovar.log'
    group4.add_argument('-g','--log',type=str,default=default,metavar='FILE',
        help='the log file [{}]'.format(default))
    default='INFO'
    group4.add_argument('-G','--loglevel',type=str,default=default,choices=['DEBUG','INFO'],
        help='the logging level [{}]'.format(default))
    default=None
    group4.add_argument('--nhx',type=str,default=default,metavar='FILE',
        help='the output file in NHX format to save the pruned tree with all variants [{}]'.format(default))
    default=None
    group4.add_argument('--NHX',type=str,default=default,metavar='FILE',
        help='the output file in NHX format to save the original tree with all variants [{}]'.format(default))
    default=None
    group4.add_argument('--nodes_vars',type=str,default=default,metavar='FILE',
        help='the output file to save SNVs/CNVs on each node [{}]'.format(default))
    default=None
    group4.add_argument('--nodes_ccf',type=str,default=default,metavar='FILE',
        help='the output file to save CCF (Cancer Cell Fraction) of each node in each sector [{}]'.format(default))
    default=None
    group4.add_argument('--cnv_profile',type=str,default=default,metavar='DIR',
        help='the output directory to save the files of CNV profile of each sector [{}]'.format(default))
    default=None
    group4.add_argument('--cnv_rc',type=str,default=default,metavar='DIR',
        help='the output directory to save the files of simulated CNV read depth of each sector [{}]'.format(default))
    default=None
    group4.add_argument('--snv_genotype',type=str,default=default,metavar='FILE',
        help='the file to save SNV genotypes for each cell [{}]'.format(default))
    default=None
    group4.add_argument('--ind_cnvs',type=str,default=default,metavar='FILE',
        help='the file to save CNVs for each cell individual [{}]'.format(default))
#    default=None
#    parser.add_argument('--haplotype_copy',type=str,default=default,metavar='FILE',
#        help='the file to save haplotype copy for each SNV')
#    default=None
#    parser.add_argument('--expands',type=str,default=default,metavar='FILE',
#        help='the basename of the file to output the snv and segment data for EXPANDS [{}]'.format(default))
    default=None
    group4.add_argument('--map',type=check_folder,default=default,metavar='DIR',
        help='directory to output the map file for each sector, which contain the relationship between tip nodes and original samples [{}]'.format(default))
    default=None
    group4.add_argument('--chain',type=check_folder,default=default,metavar='DIR',
        help='directory to output chain files for each sample [{}]'.format(default))
    args=parser.parse_args()

###### figure out the simulation setting for each chroms
#1. The setting in configure YAML file will override the setting in command line.
#2. In the configure file, the setting for individual chr will override the setting of genome.
    if args.trunk_snv_rate==None:
        args.trunk_snv_rate=args.snv_rate
    if args.trunk_cnv_rate==None:
        args.trunk_cnv_rate=args.cnv_rate
    final_chroms_cfg={}
    final_chroms_cfg[args.name]={}
    final_chroms_cfg['order']=[args.name]
    for parameter in CFG_PARAMS:
        final_chroms_cfg[args.name][parameter]=getattr(args,parameter)
    max_ploidy=len(final_chroms_cfg[args.name]['parental'])

    if args.config:
        config={}
        final_chroms_cfg={}
        with open(args.config,'r') as configfile:
            config=yaml.safe_load(configfile)
        check_config_file(config=config)
        final_chroms_cfg['order']=[list(x.keys())[0] for x in config['chromosomes']]
        max_ploidy=len(config['genome']['parental'])
        undefined_snv_rate=config['genome']['snv_rate']
        undefined_snv_rate_length=config['genome']['length']
        undefined_cnv_rate=config['genome']['cnv_rate']
        undefined_cnv_rate_length=config['genome']['length']
        undefined_trunk_snv_rate=config['genome']['trunk_snv_rate']
        undefined_trunk_snv_rate_length=config['genome']['length']
        undefined_trunk_cnv_rate=config['genome']['trunk_cnv_rate']
        undefined_trunk_cnv_rate_length=config['genome']['length']
        for chroms in config['chromosomes']:
            for chroms_cfg in chroms.values():
                if 'snv_rate' in chroms_cfg:
                    undefined_snv_rate-=chroms_cfg['snv_rate']
                    undefined_snv_rate_length-=chroms_cfg['length']
                if 'cnv_rate' in chroms_cfg:
                    undefined_cnv_rate-=chroms_cfg['cnv_rate']
                    undefined_cnv_rate_length-=chroms_cfg['length']
                if 'trunk_snv_rate' in chroms_cfg:
                    undefined_trunk_snv_rate-=chroms_cfg['trunk_snv_rate']
                    undefined_trunk_snv_rate_length-=chroms_cfg['length']
                if 'trunk_cnv_rate' in chroms_cfg:
                    undefined_trunk_cnv_rate-=chroms_cfg['trunk_cnv_rate']
                    undefined_trunk_cnv_rate_length-=chroms_cfg['length']
        for chroms in config['chromosomes']:
            for chroms_n,chroms_cfg in chroms.items():
                final_chroms_cfg[chroms_n]={}
                final_chroms_cfg[chroms_n]['snv_rate']=chroms_cfg.get('snv_rate',chroms_cfg['length']/undefined_snv_rate_length*undefined_snv_rate)
                final_chroms_cfg[chroms_n]['cnv_rate']=chroms_cfg.get('cnv_rate',chroms_cfg['length']/undefined_cnv_rate_length*undefined_cnv_rate)
                final_chroms_cfg[chroms_n]['trunk_snv_rate']=chroms_cfg.get('trunk_snv_rate',chroms_cfg['length']/undefined_trunk_snv_rate_length*undefined_trunk_snv_rate)
                final_chroms_cfg[chroms_n]['trunk_cnv_rate']=chroms_cfg.get('trunk_cnv_rate',chroms_cfg['length']/undefined_trunk_cnv_rate_length*undefined_trunk_cnv_rate)
                for parameter in set(CFG_PARAMS)-set(['snv_rate','cnv_rate','trunk_snv_rate','trunk_cnv_rate']):
                    final_chroms_cfg[chroms_n][parameter]=chroms_cfg.get(parameter,config['genome'][parameter])
                if 'parental' in chroms_cfg and len(chroms_cfg['parental'])>max_ploidy:
                    max_ploidy=len(chroms_cfg['parental'])

###### logging and random seed setting
    logging.basicConfig(filename=args.log, filemode='w',
        format='[%(asctime)s] %(levelname)s: %(message)s',
        datefmt='%m-%d %H:%M:%S',level=args.loglevel)
    argv_copy=sys.argv[:]
    argv_copy.insert(1,'phylovar')
    logging.info(' Command: %s',' '.join(argv_copy))
    if args.random_seed==None:
        seed=random_int()
    else:
        seed=args.random_seed
    logging.info(' Random seed: %s',seed)
    numpy.random.seed(seed)

###### build tree from newick string
    newick=''
    with open(args.tree) as input:
        for line in input:
            newick+=line.rstrip()
    mytree=psite.tree.newick2tree(newick)
    if args.trunk_length:
        mytree.lens=args.trunk_length

###### original_tree
    original_tree=copy.deepcopy(mytree)

###### In order to get the mytree.tipnode_leaves, we will prune the tree in all situation.
    leaves_number=mytree.leaves_counting()
    leaves_names=mytree.leaves_naming()
    sectors={}
    if args.affiliation:
        sectors=read_affiliation(args.affiliation)
        for sector in sectors:
            invalid=set(sectors[sector]['members'])-set(leaves_names)
            if invalid:
                raise AffiliationFileError("Can not find the cells below on your tree:\n{}".format(','.join([str(x) for x in invalid])))
    sectors[WHOLET]={'purity':args.purity,'depth':args.depth,'prune_p':args.prune,'members':set(mytree.leaves_naming())}
    for sector in sectors:
        sectors[sector]['prune_n']=sectors[sector]['prune_p']*len(sectors[sector]['members'])
    logging.info(' Start pruning ...')
    mytree.prune(sectors=sectors)
    mytree.collect_sectors_nodes(sectors=sectors)

##### get the ccf of each node in each sector
    if args.nodes_ccf:
        sectors_size={}
        for sector,info in sectors.items():
            sectors_size[sector]=len(info['members'])
        nodes_ccf={}
        mytree.nodes_ccf(sectors_size=sectors_size,nodes_ccf=nodes_ccf)
        with open(args.nodes_ccf,'w') as output:
            ordered_sectors=sorted(set(sectors_size)-set([WHOLET]))
            ordered_sectors.append(WHOLET)
            output.write('#node\t{}\n'.format('\t'.join(ordered_sectors)))
            for node in sorted(nodes_ccf):
                ordered_ccf=[nodes_ccf[node][sector] for sector in ordered_sectors]
                output.write('{}\t{}\n'.format(node,'\t'.join([str(x) for x in ordered_ccf])))

######
    tipnode_leaves=mytree.tipnode_leaves
    tipnode_list=list(tipnode_leaves.keys())
    tipnode_list.sort()
    leaf_tipnode={}
    leaves_names=[]
    for tipnode,leaves in tipnode_leaves.items():
        leaves_names.extend(leaves)
        for leaf in leaves:
            leaf_tipnode[leaf]=tipnode
    leaves_names.sort()
    logging.info(' There are %s leaves on your input tree.',len(leaves_names))
    logging.info(' After pruning, there are %s tip nodes on the tree.',len(tipnode_list))

###### output the map of tip_node(after pruning):leaf
    if args.chain!=None:
        os.mkdir(args.chain,mode=0o755)
    if args.map!=None:
        os.mkdir(args.map,mode=0o755)
        for sector in sectors:
            with open(os.path.join(args.map,'{}.tipnode.map'.format(sector)),'w') as tipnode_samples_map_f:
                tipnode_samples_map_f.write('#tip_node\tcell_count\tcells\n')
                for tip_node in tipnode_list:
                    focal_members=sectors[sector]['members'].intersection(set(tipnode_leaves[tip_node]))
                    if len(focal_members):
                        tipnode_samples_map_f.write('{}\t{}\t'.format(tip_node,len(focal_members)))
                        tipnode_samples_map_f.write(','.join(sorted(focal_members)))
                        tipnode_samples_map_f.write('\n')

###### add trunk vars if supplied
    trunk_snvs={}
    trunk_cnvs={}
    if args.trunk_vars!=None:
        trunk_snvs,trunk_cnvs=psite.trunk_vars.classify_vars(
            args.trunk_vars,final_chroms_cfg,leaves_number,mytree)

###### open all required output file and output the headers 
    sectors_snvs_dir=args.snv
    os.mkdir(sectors_snvs_dir,mode=0o755)
    sectors_cnvs_dir=args.cnv
    os.mkdir(sectors_cnvs_dir,mode=0o755)
    for sector,info in sectors.items():
        info['snv_file']=open(os.path.join(sectors_snvs_dir,'{}.snv'.format(sector)),'w')
        info['snv_file'].write('#chr\tstart\tend\tform\tfrequency')
        if info['depth']!=None:
            info['snv_file'].write('\trcount\trfreq\n')
        else:
            info['snv_file'].write('\n')
        info['cnv_file']=open(os.path.join(sectors_cnvs_dir,'{}.cnv'.format(sector)),'w')
        info['cnv_file'].write('#chr\tstart\tend\tcopy\tcarrier\n')

    if args.cnv_profile!=None:
        sectors_cnv_prof_dir=args.cnv_profile
        os.mkdir(sectors_cnv_prof_dir,mode=0o755)
        for sector,info in sectors.items():
            info['cnv_profile_file']=open(os.path.join(sectors_cnv_prof_dir,'{}.cnv_prof'.format(sector)),'w')
            info['cnv_profile_file'].write('#chr\tstart\tend\tparental0_cn\tparental1_cn\ttotal_cn\n')

    if args.cnv_rc!=None:
        sectors_cnv_rc_dir=args.cnv_rc
        os.mkdir(sectors_cnv_rc_dir,mode=0o755)
        for sector,info in sectors.items():
            if info['depth']!=None:
                info['cnv_rc_file']=open(os.path.join(sectors_cnv_rc_dir,'{}.cnv_rc'.format(sector)),'w')
                info['cnv_rc_file'].write('#chr\tstart\tend\tparental0_rc\tparental1_rc\ttotal_rc\n')

    if args.snv_genotype!=None:
        genotype_file=open(args.snv_genotype,'w')
        genotype_file.write('#chr\tstart\tend\tform\t{}\n'.format('\t'.join(tipnode_list)))

    if args.ind_cnvs!=None:
        ind_cnvs_file=open(args.ind_cnvs,'w')
        ind_cnvs_file.write('#cell\tparental\tchr\tstart\tend\tcopy\n')

#    if args.haplotype_copy!=None:
#        parental_copy_file=open(args.haplotype_copy,'w')
#        parental_copy_file.write('#chr\tpos\t{}\n'.format('\t'.join(['haplotype'+str(x) for x in range(max_ploidy)])))

#    if args.expands != None:
#        expands_snps_file=open(args.expands+'.snps','w')
#        expands_snps_file.write('chr\tstartpos\tAF_Tumor\tPN_B\n')
#        expands_segs_file=open(args.expands+'.segs','w')
#        expands_segs_file.write("chr\tstartpos\tendpos\tCN_Estimate\n")

###### cnvl_dist
    cnvl_dist=None
    if args.cnvl_dist:
        cnvl_dist=read_cnvl_dist(args.cnvl_dist)

###### simulate variants for each chroms
    sex_chrs=set()
    if args.sex_chr:
        sex_chrs=set(args.sex_chr.split(','))
    all_nodes_vars={}
    for chroms in final_chroms_cfg['order']:
        chroms_cfg=final_chroms_cfg[chroms]
        if cnvl_dist==None:
            check_cnv_length_cfg(chroms=chroms,cnv_length_beta=chroms_cfg['cnv_length_beta'],
                cnv_length_max=chroms_cfg['cnv_length_max'],chr_length=chroms_cfg['length'])
        cn_dist_cfg=cn_dist(copy_max=chroms_cfg['copy_max'],copy_parameter=chroms_cfg['copy_parameter'])
        tstv_dist_cfg=tstv_dist(tstv=chroms_cfg['tstv'])
        logging.info(' Start the simulation for chromosome: %s',chroms)
#I need the normal_dosage to adjust the frequency of snv under under different purity
        for sector,info in sectors.items():
            if chroms in sex_chrs and len(sex_chrs)==2:
                n=1
            else:
                n=2
            tumor_cells=len(info['members'])
            total_cells=round(tumor_cells/info['purity'])
            normal_cells=total_cells-tumor_cells
            info['standard_total_dosage']=total_cells*n
            info['normal_dosage']=normal_cells*n

        (nodes_vars,tipnode_snv_alts,tipnode_snv_refs,tipnode_cnvs,
            )=mytree.snvs_freq_cnvs_profile(
                parental=chroms_cfg['parental'],
                snv_rate=chroms_cfg['snv_rate'],
                cnv_rate=chroms_cfg['cnv_rate'],
                trunk_snv_rate=chroms_cfg['trunk_snv_rate'],
                trunk_cnv_rate=chroms_cfg['trunk_cnv_rate'],
                del_prob=chroms_cfg['del_prob'],
                tandem_prob=chroms_cfg['tandem_prob'],
                cnv_length_beta=chroms_cfg['cnv_length_beta'],
                cnv_length_max=chroms_cfg['cnv_length_max'],
                cn_dist_cfg=cn_dist_cfg,
                tstv_dist_cfg=tstv_dist_cfg,
                trunk_snvs=trunk_snvs.get(chroms,{}),
                trunk_cnvs=trunk_cnvs.get(chroms,{}),
                length=chroms_cfg['length'],
                chain=args.chain,
                chroms=chroms,
                sectors=sectors,
                wholeT=WHOLET,
                cnvl_dist=cnvl_dist,
            )
        snvs_alt_total=sectors[WHOLET]['snvs_alt_total']
        cnvs=sectors[WHOLET]['cnvs']
        if args.nhx or args.NHX or args.nodes_vars:
            all_nodes_vars=psite.tree.merge_two_dict_set(dict1=all_nodes_vars,dict2=nodes_vars)

        if args.snv_genotype!=None:
            for pos,mutation,alt,total in snvs_alt_total:
                genotype_file.write('{}\t{}\t{}\t{}\t{}\n'.format(chroms,pos,pos+1,mutation,
                    '\t'.join([str(tipnode_snv_alts[tipnode][pos])+':'+str(tipnode_snv_refs[tipnode][pos]) for tipnode in tipnode_list])))

        if args.ind_cnvs!=None:
            for tipnode in tipnode_list:
                for cnv in tipnode_cnvs[tipnode]:
                    cnv_copy='+{}'.format(cnv['copy']) if cnv['copy']>0 else str(cnv['copy'])
                    ind_cnvs_file.write('{}\n'.format('\t'.join([str(x) for x in [tipnode,cnv['parental'],chroms,cnv['start'],cnv['end'],cnv_copy]])))

#        if args.haplotype_copy!=None:
#            for snv in hap_local_copy_for_all_snvs:
#                parental_copy_file.write('{}\t{}\n'.format(chroms,'\t'.join([str(x) for x in snv])))

        for sector,info in sectors.items():
            for pos,mutation,alt,total in info['snvs_alt_total']:
                freq=alt/total
                info['snv_file'].write('{}\t{}\t{}\t{}\t{}'.format(chroms,pos,pos+1,mutation,round(freq,4)))
                if info['depth']!=None:
                    expected_total_dp=info['depth']*total/info['standard_total_dosage']
                    total_dp,b_allele_dp=psite.tree.simulate_sequence_coverage(expected_total_dp,freq)
                    if total_dp!=0:
                        rfreq=round(b_allele_dp/total_dp,4)
                    else:
                        rfreq='-'
                    info['snv_file'].write('\t{}:{}\t{}\n'.format(b_allele_dp,total_dp,rfreq))
                else:
                    info['snv_file'].write('\n')
            for cnv in info['cnvs']:
                cnv_copy='+{}'.format(cnv['copy']) if cnv['copy']>0 else str(cnv['copy'])
                info['cnv_file'].write('{}\t{}\t{}\t{}\t{}\n'.format(chroms,cnv['start'],cnv['end'],cnv_copy,cnv['leaves_count']))

        if chroms in sex_chrs and len(sex_chrs)==2: # haploid sex chromosomes
            for sector,info in sectors.items():
                for seg in info['cnv_profile']:
#cnv_profile means the local copy of each segment across the cell population of the sample (normal+tumor)
                    seg[2]=seg[2]+info['normal_dosage']
                    seg[3]=seg[3]+0
                    seg[4]=seg[4]+info['normal_dosage']
        else:
            for sector,info in sectors.items():
                for seg in info['cnv_profile']:
                    seg[2]=seg[2]+round(info['normal_dosage']/2)
                    seg[3]=seg[3]+round(info['normal_dosage']/2)
                    seg[4]=seg[4]+info['normal_dosage']

        if args.cnv_profile!=None:
            for sector,info in sectors.items():
                for seg in info['cnv_profile']:
#cnv_profile means the local copy of each segment across the cell population of the sample (normal+tumor)
                    info['cnv_profile_file'].write('{}\n'.format('\t'.join([str(x) for x in [chroms]+seg])))

        if args.cnv_rc!=None:
            read_length=args.rlen
            for sector,info in sectors.items():
                if info['depth']!=None:
                    temp=[]
                    for seg in info['cnv_profile']:
                        if temp:
                            if seg[2]==temp[-1][2] and seg[3]==temp[-1][3] and seg[4]==temp[-1][4] and seg[0]==temp[-1][1]:
                                temp[-1][1]=seg[1]
                            else:
                                temp.append(seg)
                        else:
                            temp.append(seg)
                    for seg in temp:
                        start,end,parental0_cn,parental1_cn,total_cn=seg
                        expected_total_dp=info['depth']*total_cn/info['standard_total_dosage']
                        total_rc,parental0_rc,parental1_rc=psite.tree.simulate_cnv_rc(
                            mean_coverage=expected_total_dp,
                            parental0_cn=parental0_cn,
                            parental1_cn=parental1_cn,
                            seg_length=end-start,
                            read_length=read_length)
                        info['cnv_rc_file'].write('{}\n'.format('\t'.join([str(x) for x in (chroms,start,end,parental0_rc,parental1_rc,total_rc)])))

##output for expands
#        if args.expands != None:
#            for pos,mutation,freq in snvs_freq:
#                total_dp,b_allele_dp=psite.tree.simulate_sequence_coverage(args.depth,freq)
#                expands_snps_file.write('{}\t{}\t{}\t{}\n'.format(chroms,pos,b_allele_dp/total_dp,0))
#
##in the segment input for expands
##CN_Estimate - the copy number estimated for each segment (average value across all subpopulations in the sample)
#            for start,end,copy in cnv_profile:
#                expands_segs_file.write('{}\t{}\t{}\t{}\n'.format(chroms,start,end,copy/leaves_number))

###### close all opened files
    for sector,info in sectors.items():
        info['snv_file'].close()
        info['cnv_file'].close()

    if args.cnv_profile!=None:
        for sector,info in sectors.items():
            info['cnv_profile_file'].close()

    if args.cnv_rc!=None:
        for sector,info in sectors.items():
            if info['depth']!=None:
                info['cnv_rc_file'].close()

    if args.snv_genotype!=None:
        genotype_file.close()

    if args.ind_cnvs!=None:
        ind_cnvs_file.close()

#    if args.haplotype_copy!=None:
#        parental_copy_file.close()

#    if args.expands != None:
#        expands_snps_file.close()
#        expands_segs_file.close()

    if args.nhx or args.NHX:
        nodes_nSNV={}
        nodes_nAMP={}
        nodes_nDEL={}
        for node,variants in all_nodes_vars.items():
            nSNV,nAMP,nDEL=0,0,0
            for var in variants:
                form=var.split('#')[4]
                if form.startswith('+'):
                    nAMP+=1
                elif form.startswith('-'):
                    nDEL+=1
                else:
                    nSNV+=1
            nodes_nSNV[node]=nSNV
            nodes_nAMP[node]=nAMP
            nodes_nDEL[node]=nDEL

        if args.nhx:
            mytree.attach_info(attr='vars',info=all_nodes_vars)
            mytree.attach_info(attr='nSNV',info=nodes_nSNV,null=0)
            mytree.attach_info(attr='nAMP',info=nodes_nAMP,null=0)
            mytree.attach_info(attr='nDEL',info=nodes_nDEL,null=0)
            with open(args.nhx,'w') as tree_data_file:
                tree_data_file.write('{};\n'.format(mytree.tree2nhx(with_lens=True,attrs=['nodeid','vars','nSNV','nAMP','nDEL'])))
        if args.NHX:
            original_tree.attach_info(attr='vars',info=all_nodes_vars)
            original_tree.attach_info(attr='nSNV',info=nodes_nSNV,null=0)
            original_tree.attach_info(attr='nAMP',info=nodes_nAMP,null=0)
            original_tree.attach_info(attr='nDEL',info=nodes_nDEL,null=0)
            with open(args.NHX,'w') as tree_data_file:
                tree_data_file.write('{};\n'.format(original_tree.tree2nhx(with_lens=True,attrs=['nodeid','vars','nSNV','nAMP','nDEL'])))

#output SNVs/CNVs on each node
    if args.nodes_vars:
        with open(args.nodes_vars,'w') as nodes_vars_file:
            nodes_vars_file.write('#node\tchr\thap\tstart\tend\tvar\n')
            for node in sorted(all_nodes_vars.keys(),key=lambda x: int(x[4:])):
                vars_list=[x.split('#') for x in all_nodes_vars[node]]
                vars_list=sorted(vars_list,key=lambda x:(x[0],int(x[2]),int(x[3])))
                for var in vars_list:
                    nodes_vars_file.write('{}\t{}\n'.format(node,'\t'.join(var)))
    t1 = time.time()
    print ("Total time running {}: {} seconds".format
      (prog, str(t1-t0)))
