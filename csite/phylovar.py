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
import csite.trunk_vars
import csite.tree
from csite.vcf2fa import check_sex
#import pickle

#handle the error below
#python | head == IOError: [Errno 32] Broken pipe 
from signal import signal, SIGPIPE, SIG_DFL 
signal(SIGPIPE,SIG_DFL) 

#TODO: check whether parental is parental. sometimes it's mean haplotype...
#TODO: SNV true_freq
#TODO: rewrite the description of the output of SNVs 

#I defined those two parameters as global variables. As they will be used in function
#random_int and check_config_file, which are also used in allinone.py.
largest=2**32
cfg_params={'snv_rate':float,
            'cnv_rate':float,
            'del_prob':float,
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
    Let's use this to generate a integers can be used.
    '''
    return numpy.random.randint(largest) 

def check_seed(value=None):
    ivalue=int(value)
#2**32: Must be convertible to 32 bit unsigned integers.
    if not 0<=ivalue<=largest: 
        raise argparse.ArgumentTypeError("{} is an invalid value for --random_seed. ".format(value)+
            "It should be an integer between 0 and {}.".format(largest))
    return ivalue

def check_prune(value=None):
    ivalue=int(float(value))
    if not ivalue>0: 
        raise argparse.ArgumentTypeError("{} is an invalid value for --prune. ".format(value)+
            "It should be an positive integer and less than the number of leaves of the tree.")
    return ivalue

def check_proportion(value=None):
    fvalue=float(value)
    if not 0<fvalue<1: 
        raise argparse.ArgumentTypeError("{} is an invalid value for --prune_proportion. ".format(value)+
            "It should be an float between 0 and 1.")
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

def check_folder(directory=None):
    good_charactors=re.compile('^[0-9a-zA-Z/_\-]+$') 
    if not good_charactors.match(directory):
        raise argparse.ArgumentTypeError("'{}' is an invalid string for --chain. ".format(directory)+
            "Please only use number, alphabet and _/- in the directory name.")
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
    Check the validation of your config file.
    1. There must be only two section in config file: genome and chromosomes.
    2. Every parameter in cfg_params must have a pair of key:value in genome section.
    3. Every chromosome at least has the key 'length'.
    4. The total length of all chromosomes must equal the length of genome.
    5. If the SNV rate of each chromosome has been specified, the sum of them should be equal the SNV rate of genome.
    6. If the SNV rate of some chromosomes (not all) have been specified, the sum of them should be <= the SNV rate of genome.
    7. The CNV rate should satisfy the same criteria as SNV rate.
    '''
    if not (isinstance(config,dict) and set(config)==set(['genome','chromosomes'])):
        raise ConfigFileError('Check your config file. The format is not correct.')
    if not isinstance(config['genome'],dict):
        raise ConfigFileError('Check your config file. The format in genome section is not correct.')
    lack=set(cfg_params)-set(config['genome'])
    if lack!=set():
        raise ConfigFileError("'{}' are required in the genome section of config file.".format(','.join([str(x) for x in lack])))
    over=set(config['genome'])-set(cfg_params)
    if over!=set():
        raise ConfigFileError("'{}' are not acceptable parameters in config file.".format(','.join([str(x) for x in over])))

    for parameter in cfg_params:
        assert isinstance(config['genome'][parameter],cfg_params[parameter]),\
            "'{}' in genome section of config file should be a {}!".format(parameter,cfg_params[parameter].__name__)

    if not isinstance(config['chromosomes'],list):
        raise ConfigFileError('Check your config file. The format in chromosomes section is not correct.')
    total_chroms_length=0
    total_snv_rate=0
    missing_snv_rate=0
    total_cnv_rate=0
    missing_cnv_rate=0
    for chroms in config['chromosomes']: 
        if not isinstance(chroms,dict) or len(chroms)>1:
            raise ConfigFileError('Check your config file. The format in chromosomes section is not correct.')
        for chroms_n,chroms_cfg in chroms.items():
            assert isinstance(chroms_n,str),'The name of chromosome {} in config file should be a str!'.format(chroms_n)
            over=set(chroms_cfg)-set(cfg_params)
            if over!=set():
                raise ConfigFileError("'{}' are not acceptable parameters ".format(','.join([str(x) for x in over]))+
                    'in the section of chromosome {} in your config file!'.format(chroms))
            for parameter in chroms_cfg.keys():
                assert isinstance(chroms_cfg[parameter],cfg_params[parameter]),\
                    "'{}' in of chromosome {} in config file should be a {}!".format(parameter,chroms_n,cfg_params[parameter].__name__)
            if 'length' not in chroms_cfg:
                raise ConfigFileError("Couldn't find the length of chromosome:{}.".format(chroms_n))
            total_chroms_length+=chroms_cfg['length']
            if 'snv_rate' in chroms_cfg:
                total_snv_rate+=chroms_cfg['snv_rate']
            else:
                missing_snv_rate=1
            if 'cnv_rate' in chroms_cfg:
                total_cnv_rate+=chroms_cfg['cnv_rate']
            else:
                missing_cnv_rate=1

    if config['genome']['length']!=total_chroms_length:
        raise ConfigFileError('In your config file, the length of genome is {},'.format(str(config['genome']['length']))+
            'But the total length of all chromosomes are {}'.format(str(total_chroms_length)))
    if missing_snv_rate==1:
        if total_snv_rate>config['genome']['snv_rate']:
            raise ConfigFileError('Check your config file, the sum of snv_rate in chromosmomes section is larger than the rate in genome section!')
    else:
        if total_snv_rate!=config['genome']['snv_rate']:
            raise ConfigFileError('Check your config file, the sum of snv_rate in chromosmomes section is not equal the rate in genome section!')
    if missing_cnv_rate==1:
        if total_cnv_rate>config['genome']['cnv_rate']:
            raise ConfigFileError('Check your config file, the sum of cnv_rate in chromosmomes section is larger than the rate in genome section!')
    else:
        if total_cnv_rate!=config['genome']['cnv_rate']:
            raise ConfigFileError('Check your config file, the sum of cnv_rate in chromosmomes section is not equal the rate in genome section!')

class ConfigFileError(Exception):
    pass

#use kernprof -l -v script.py to profile
#@profile
def main(progname=None):
    parser=argparse.ArgumentParser(
        description='Simulate SNVs/CNVs on a coalescent tree in newick format',
        prog=progname if progname else sys.argv[0])
    group1=parser.add_argument_group('Input Files')#, 'input files')
    group1.add_argument('-t','--tree',required=True,metavar='FILE',
        help='a file containing ONE tree in newick format')
    default=None
    group1.add_argument('--trunk_vars',type=str,default=default,metavar='FILE',
        help='a file containing truncal variants predefined by user [{}]'.format(default))
    default=None
    group1.add_argument('--config',type=str,default=default,metavar='FILE',
        help='a YAML file which contains the configuration of somatic variant simulation. '+
            '-n/-r/-R/-d/-l/-L/-c/-C/-p/--tstv/--length will be ignored. [{}]'.format(default))
    group2=parser.add_argument_group('Simulation Parameters (can be set in config YAML)')
    default='1'
    group2.add_argument('-n','--name',type=str,default=default,metavar='STR',
        help='the name of the sequence to be simulated [{}]'.format(default))
    default=300
    group2.add_argument('-r','--snv_rate',type=float,default=default,metavar='FLOAT',
        help='the muation rate of SNVs [{}]'.format(default))
    default=3
    group2.add_argument('-R','--cnv_rate',type=float,default=default,metavar='FLOAT',
        help='the muation rate of CNVs [{}]'.format(default))
    default=0.5
    group2.add_argument('-d','--del_prob',type=float,default=default,metavar='FLOAT',
        help='the probability of being deletion for a CNV mutation [{}]'.format(default))
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
    group3=parser.add_argument_group('Other Simulation Parameters (can NOT be set in config YAML)')
    default=0
    group3.add_argument('-x','--prune',type=check_prune,default=default,metavar='INT',
        help='trim all the children of the nodes with equal or less than this number of leaves [{}]'.format(default))
    default=0.0
    group3.add_argument('-X','--prune_proportion',type=check_proportion,default=default,metavar='FLOAT',
        help='trim all the children of the nodes with equal or less than this proportion of total leaves [{}]'.format(default))
    default=None
    group3.add_argument('-s','--random_seed',type=check_seed,metavar='INT',
        help='the seed for random number generator [{}]'.format(default))
    default=0
    group3.add_argument('--trunk_length',type=float,default=default,metavar='FLOAT',
        help='the length of the trunk [{}]'.format(default))
    default=0.8
    group3.add_argument('--purity',type=check_purity,default=default,metavar='FLOAT',
        help='the proportion of tumor cells in simulated tumor sample [{}]'.format(default))
    default=None
    group3.add_argument('--sex_chr',type=check_sex,default=default,metavar='STR',
        help='sex chromosomes of the genome (separated by comma) [{}]'.format(default))
    group4=parser.add_argument_group('Output Related Parameters')
    default='phylovar.snvs'
    group4.add_argument('-S','--snv',type=str,default=default,metavar='FILE',
        help='the output file to save SNVs [{}]'.format(default))
    default='phylovar.cnvs'
    group4.add_argument('-V','--cnv',type=str,default=default,metavar='FILE',
        help='the output file to save CNVs [{}]'.format(default))
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
    group4.add_argument('--cnv_profile',type=str,default=default,metavar='FILE',
        help='the file to save CNVs profile [{}]'.format(default))
    group4.add_argument('--snv_genotype',type=str,metavar='FILE',
        help='the file to save SNV genotypes for each cell')
    group4.add_argument('--ind_cnvs',type=str,metavar='FILE',
        help='the file to save CNVs for each cell individual')
#    parser.add_argument('--haplotype_copy',type=str,
#        help='the file to save haplotype copy for each SNV')
#    default=None
#    parser.add_argument('--expands',type=str,default=default,
#        help='the basename of the file to output the snv and segment data for EXPANDS [{}]'.format(default))
    default=None
    group4.add_argument('--map',type=str,default=default,metavar='FILE',
        help='the map file to save the relationship between tip nodes and original samples [{}]'.format(default))
    default=None
    group4.add_argument('--chain',type=check_folder,default=default,metavar='DIR',
        help='directory to output chain files for each sample [{}]'.format(default))
    args=parser.parse_args()

###### figure out the simulation setting for each chroms
#1. The setting in configure YAML file will override the setting in command line.
#2. In the configure file, the setting for individual chr will override the setting of genome.
    final_chroms_cfg={} 
    final_chroms_cfg[args.name]={}
    final_chroms_cfg['order']=[args.name]
    for parameter in cfg_params:
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
        for chroms in config['chromosomes']: 
            for chroms_cfg in chroms.values():
                if 'snv_rate' in chroms_cfg:
                    undefined_snv_rate-=chroms_cfg['snv_rate']
                    undefined_snv_rate_length-=chroms_cfg['length']
                if 'cnv_rate' in chroms_cfg:
                    undefined_cnv_rate-=chroms_cfg['cnv_rate']
                    undefined_cnv_rate_length-=chroms_cfg['length']
        for chroms in config['chromosomes']: 
            for chroms_n,chroms_cfg in chroms.items():
                final_chroms_cfg[chroms_n]={}
                final_chroms_cfg[chroms_n]['snv_rate']=chroms_cfg.get('snv_rate',chroms_cfg['length']/undefined_snv_rate_length*undefined_snv_rate)
                final_chroms_cfg[chroms_n]['cnv_rate']=chroms_cfg.get('cnv_rate',chroms_cfg['length']/undefined_cnv_rate_length*undefined_cnv_rate)
                for parameter in set(cfg_params)-set(['snv_rate','cnv_rate']):
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
#TODO: We should do make sure the newick string is valided before processing it.
    newick=''
    with open(args.tree) as input:
        for line in input:
            newick+=line.rstrip()
    mytree=csite.tree.newick2tree(newick)
    if args.trunk_length:
        mytree.lens=args.trunk_length
#################original_tree
    original_tree=copy.deepcopy(mytree)
    leaves_number=mytree.leaves_counting()

####### prune the tree if required
    if args.prune>0 and args.prune_proportion>0:
        raise argparse.ArgumentTypeError("Use either --prune or --prune_proportion. Do not use both!")
    elif args.prune>0:
        if not args.prune<leaves_number:
            raise argparse.ArgumentTypeError("There are only {} leaves on the tree. Can not prune {} leaves.".format(
                leaves_number,args.prune))
        elif args.prune>=1:
            mytree.prune(tips=args.prune)
    elif args.prune_proportion>0:
        trim=leaves_number*args.prune_proportion
        if trim>=1:
            mytree.prune(tips=trim)
    else:
        mytree.prune(tips=0.5)
    tipnode_leaves=mytree.tipnode_leaves
    leaf_tipnode={}
    leaves_names=[]
    for tipnode,names in tipnode_leaves.items():
        leaves_names.extend(names)
        for name in names:
            leaf_tipnode[name]=tipnode
    leaves_names.sort()
    logging.info(' There are %s leaves on your input tree.',len(leaves_names))
    if args.prune>0 or args.prune_proportion>0:
        logging.info(' After pruning, there are %s tip nodes on the tree.',len(tipnode_leaves))

###### output the map of tip_node(after pruning):leaf
    if args.chain:
        os.mkdir(args.chain,mode=0o755)
    if args.map:
        with open(args.map,'w') as tipnode_samples_map_f:
            tipnode_samples_map_f.write('#tip_node\tcell_count\tcells\n')
            for tip_node in sorted(tipnode_leaves.keys()):
                tipnode_samples_map_f.write('{}\t{}\t'.format(tip_node,len(tipnode_leaves[tip_node])))
                tipnode_samples_map_f.write(','.join(sorted(tipnode_leaves[tip_node])))
                tipnode_samples_map_f.write('\n')



###### add trunk vars if supplied
    trunk_snvs={}
    trunk_cnvs={}
    if args.trunk_vars!=None:
        trunk_snvs,trunk_cnvs=csite.trunk_vars.classify_vars(
            args.trunk_vars,final_chroms_cfg,leaves_number,mytree)

###### open all required output file and output the headers 
    cnv_file=open(args.cnv,'w')
    cnv_file.write('#chr\tstart\tend\tcopy\tcarrier\n')
    snv_file=open(args.snv,'w')
    snv_file.write('#chr\tstart\tend\tform\tfrequency\n')

    if args.cnv_profile!=None:
        cnv_profile_file=open(args.cnv_profile,'w')
        cnv_profile_file.write('#chr\tstart\tend\tlocal_cp\n')

    if args.snv_genotype!=None:
        genotype_file=open(args.snv_genotype,'w')
        genotype_file.write('#chr\tpos\t{}\n'.format('\t'.join(leaves_names)))
    
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

###### simulate variants for each chroms
    sex_chrs=set()
    if args.sex_chr:
        sex_chrs=set(args.sex_chr.split(','))
    all_nodes_vars={}
    for chroms in final_chroms_cfg['order']:
        chroms_cfg=final_chroms_cfg[chroms]
        check_cnv_length_cfg(chroms=chroms,cnv_length_beta=chroms_cfg['cnv_length_beta'],
            cnv_length_max=chroms_cfg['cnv_length_max'],chr_length=chroms_cfg['length'])
        cn_dist_cfg=cn_dist(copy_max=chroms_cfg['copy_max'],copy_parameter=chroms_cfg['copy_parameter'])
        tstv_dist_cfg=tstv_dist(tstv=chroms_cfg['tstv'])
        logging.info(' Start the simulation for chromosome: %s',chroms)
#I need use normal_dosage to adjust the frequency of snv under under different purity
        normal_dosage=leaves_number*2*(1-args.purity)/args.purity
        if chroms in sex_chrs and len(sex_chrs)==2:
            normal_dosage=leaves_number*1*(1-args.purity)/args.purity
        
        (snvs_freq,cnvs,cnv_profile,nodes_vars,
            tipnode_snv_alts,tipnode_snv_refs,tipnode_cnvs,
            )=mytree.snvs_freq_cnvs_profile(
                parental=chroms_cfg['parental'],
                snv_rate=chroms_cfg['snv_rate'],
                cnv_rate=chroms_cfg['cnv_rate'],
                del_prob=chroms_cfg['del_prob'],
                cnv_length_beta=chroms_cfg['cnv_length_beta'],
                cnv_length_max=chroms_cfg['cnv_length_max'],
                cn_dist_cfg=cn_dist_cfg,
                tstv_dist_cfg=tstv_dist_cfg,
                trunk_snvs=trunk_snvs.get(chroms,{}),
                trunk_cnvs=trunk_cnvs.get(chroms,{}),
                length=chroms_cfg['length'],
                normal_dosage=normal_dosage,
                chain=args.chain,
                chroms=chroms,
            )
        all_snvs_pos=sorted(x[0] for x in snvs_freq)
        all_nodes_vars=csite.tree.merge_two_dict_set(dict1=all_nodes_vars,dict2=nodes_vars)

        if args.snv_genotype!=None:
            for pos in all_snvs_pos:
                genotype_file.write('{}\t{}\t{}\n'.format(chroms,pos,
                    '\t'.join([str(tipnode_snv_alts[leaf_tipnode[leaf]][pos])+':'+str(tipnode_snv_refs[leaf_tipnode[leaf]][pos]) for leaf in leaves_names])))

        if args.ind_cnvs!=None:
            for leaf in leaves_names:
                for cnv in tipnode_cnvs[leaf_tipnode[leaf]]:
                    cnv_copy='+{}'.format(cnv['copy']) if cnv['copy']>0 else str(cnv['copy'])
                    ind_cnvs_file.write('{}\n'.format('\t'.join([str(x) for x in [leaf,cnv['parental'],chroms,cnv['start'],cnv['end'],cnv_copy]])))

#        if args.haplotype_copy!=None:
#            for snv in hap_local_copy_for_all_snvs:
#                parental_copy_file.write('{}\t{}\n'.format(chroms,'\t'.join([str(x) for x in snv])))

        for cnv in cnvs:
            cnv_copy='+{}'.format(cnv['copy']) if cnv['copy']>0 else str(cnv['copy'])
            cnv_file.write('{}\t{}\t{}\t{}\t{}\n'.format(chroms,cnv['start'],cnv['end'],cnv_copy,cnv['leaves_count']))

        for pos,mutation,freq in snvs_freq:
            snv_file.write('{}\t{}\t{}\t{}\t{}\n'.format(chroms,pos,pos+1,mutation,freq))

        if args.cnv_profile!=None:
            for seg in cnv_profile:
                cnv_profile_file.write('{}\n'.format('\t'.join([str(x) for x in [chroms]+seg])))

##output for expands
#        if args.expands != None:
#            for pos,mutation,freq in snvs_freq:
#                total_dp,b_allele_dp=csite.tree.simulate_sequence_coverage(args.depth,freq)
#                expands_snps_file.write('{}\t{}\t{}\t{}\n'.format(chroms,pos,b_allele_dp/total_dp,0))
#
##in the segment input for expands
##CN_Estimate - the copy number estimated for each segment (average value across all subpopulations in the sample)
#            for start,end,copy in cnv_profile:
#                expands_segs_file.write('{}\t{}\t{}\t{}\n'.format(chroms,start,end,copy/leaves_number))

###### close all opened files
    cnv_file.close()
    snv_file.close()

    if args.cnv_profile!=None:
        cnv_profile_file.close()

    if args.snv_genotype!=None:
        genotype_file.close()
    
    if args.ind_cnvs!=None:
        ind_cnvs_file.close()

#    if args.haplotype_copy!=None:
#        parental_copy_file.close()

#    if args.expands != None:
#        expands_snps_file.close()
#        expands_segs_file.close()

#TODO: Should we change pickle to json or yaml?
#http://stackoverflow.com/questions/4677012/python-cant-pickle-type-x-attribute-lookup-failed
#FIXME: right now, it does not consider the deletion effect on pre_snvs.
#FIXME: output in .nhx format
    if args.nhx!=None:
        mytree.attach_info(attr='vars',info=all_nodes_vars)
        with open(args.nhx,'w') as tree_data_file:
            tree_data_file.write('{};\n'.format(mytree.tree2nhx(with_lens=True,attrs=['nodeid','vars'])))

    if args.NHX!=None:
        original_tree.attach_info(attr='vars',info=all_nodes_vars)
        with open(args.NHX,'w') as tree_data_file:
            tree_data_file.write('{};\n'.format(original_tree.tree2nhx(with_lens=True,attrs=['nodeid','vars'])))
