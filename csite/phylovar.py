#!/usr/bin/env python3

#########################################################################
# Author: Hechuan
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
import csite.trunk_vars
import csite.tree
import yaml
#import pickle

#handle the error below
#python | head == IOError: [Errno 32] Broken pipe 
from signal import signal, SIGPIPE, SIG_DFL 
signal(SIGPIPE,SIG_DFL) 

#TODO: SNV true_freq
#rewrite the description of the output of SNVs 

largest=2**32
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
            "It should be an interger between 0 and {}.".format(largest))
    return ivalue

def check_prune(value=None):
    ivalue=int(float(value))
    if ivalue<2: 
        raise argparse.ArgumentTypeError("{} is an invalid value for --prune. ".format(value)+
            "It should be an interger >=2 and <=the number of leaves of the tree.")
    return ivalue

def check_proportion(value=None):
    fvalue=float(value)
    if not 0<=fvalue<=1: 
        raise argparse.ArgumentTypeError("{} is an invalid value for --prune_proportion. ".format(value)+
            "It should be an float between 0 and 1.")
    return fvalue

def check_tstv(value=None):
    fvalue=float(value)
    if fvalue<=0: 
        raise argparse.ArgumentTypeError("{} is an invalid value for transition/transversion ratio. ".format(value)+
            "It should be an float larger than 0.")
    return fvalue

def check_folder(directory=None):
    good_charactors=re.compile('^[0-9a-zA-Z/_\-]+$') 
    if not good_charactors.match(directory):
        raise argparse.ArgumentTypeError("{} is an invalid string for --chain. ".format(directory)+
            "Please only use number, alphabet and _/- in the directory name.")
    if os.path.exists(directory):
        raise argparse.ArgumentTypeError("{} is already exist. Delete it or use another name instead.".format(directory))
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

def check_config_file(config=None,cfg_params=None):
    '''
    Check the validation of your config file.
    '''
#TODO: we should collect all invalid keys and print all of them in one time
    if config and not isinstance(config,dict):
        raise ConfigFileError('Check your config file. The format is not correct.')
    for key in set(config.keys())-set(['genome','chromosomes']):
        raise ConfigFileError('{} is not acceptable '.format(key)+
            'as the the 1st level key in your configure file!')
    if 'genome' in config and config['genome']:
        if not isinstance(config['genome'],dict):
            raise ConfigFileError('Check your config file. The format in genome section is not correct.')
        for key in set(config['genome'].keys())-set(cfg_params):
            raise ConfigFileError('{} is not an acceptable key '.format(key)+
                'in the genome section of your configure file!')
    if 'chromosomes' in config and config['chromosomes']:
        if not isinstance(config['chromosomes'],list):
            raise ConfigFileError('Check your config file. The format in chromosomes section is not correct.')
        for chroms in config['chromosomes']: 
            if not isinstance(chroms,dict) or len(chroms)>1:
                raise ConfigFileError('Check your config file. The format in chromosomes section is not correct.')
            for chroms_n,chroms_cfg in chroms.items():
                for key in set(chroms_cfg.keys())-set(cfg_params):
                    raise ConfigFileError('{} is not an acceptable key '.format(key)+
                        'in the section of chromosome {} in your configure file!'.format(chroms))

class ConfigFileError(Exception):
    pass

#use kernprof -l -v script.py to profile
#@profile
def main(progname=None):
    parse=argparse.ArgumentParser(
        description='Simulate SNVs/CNVs on a coalescent tree in newick format',
        prog=progname if progname else sys.argv[0])
    parse.add_argument('-t','--tree',required=True,
        help='a file contains ONE tree in newick format')
    default='1'
    parse.add_argument('-n','--name',type=str,default=default,
        help='the name of the sequence to be simulated [{}]'.format(default))
    default=300
    parse.add_argument('-r','--snv_rate',type=float,default=default,
        help='the muation rate of SNVs [{}]'.format(default))
    default=3
    parse.add_argument('-R','--cnv_rate',type=float,default=default,
        help='the muation rate of CNVs [{}]'.format(default))
    default=0.5
    parse.add_argument('-d','--del_prob',type=int,default=default,
        help='the probability of being deletion for a CNV mutation [{}]'.format(default))
#https://en.wikipedia.org/wiki/Copy-number_variation
    default=20000000
    parse.add_argument('-l','--cnv_length_beta',type=int,default=default,
        help='the mean of CNVs length [{}]'.format(default))
    default=40000000
    parse.add_argument('-L','--cnv_length_max',type=int,default=default,
        help='the maximium of CNVs length [{}]'.format(default))
    default=0.5
    parse.add_argument('-c','--copy_parameter',type=float,default=default,
        help="the p parameter of CNVs' copy number distribution [{}]".format(default))
    default=5
    parse.add_argument('-C','--copy_max',type=int,default=default,
        help='the maximium ADDITIONAL copy of a CNVs [{}]'.format(default))
    default='01'
    parse.add_argument('-p','--parental',type=str,default=default,
        help='the parental to simulate [{}]'.format(default))
    default=1.0
    parse.add_argument('-P','--purity',type=float,default=default,
        help="the purity of tumor cells in the simulated sample [{}]".format(default))
    default=50
    parse.add_argument('-D','--depth',type=int,default=default,
        help='the mean depth for simulating coverage data [{}]'.format(default))
    default=0
    parse.add_argument('-x','--prune',type=check_prune,default=default,
        help='trim all the children of the nodes with equal or less than this number of tips [{}]'.format(default))
    default=0.0
    parse.add_argument('-X','--prune_proportion',type=check_proportion,default=default,
        help='trim all the children of the nodes with equal or less than this proportion of tips [{}]'.format(default))
    default=None
    parse.add_argument('-s','--random_seed',type=check_seed,
        help='the seed for random number generator [{}]'.format(default))
    default='output.snvs'
    parse.add_argument('-S','--snv',type=str,default=default,
        help='the output file to save SNVs [{}]'.format(default))
    default='output.cnvs'
    parse.add_argument('-V','--cnv',type=str,default=default,
        help='the output file to save CNVs [{}]'.format(default))
    default='output.nodes_vars'
    parse.add_argument('-N','--nodes_vars',type=str,default=default,
        help='the output file to save SNVs/CNVs on each node [{}]'.format(default))
    default='output.named_tree.nhx'
    parse.add_argument('-T','--named_tree',type=str,default=default,
        help='the output file in NHX format to save the tree with all nodes named [{}]'.format(default))
    default='phylovar.log'
    parse.add_argument('-g','--log',type=str,default=default,
        help='the log file [{}]'.format(default))
    default='INFO'
    parse.add_argument('-G','--loglevel',type=str,default=default,choices=['DEBUG','INFO'],
        help='the logging level [{}]'.format(default))
    default='output.cnv.profile'
    parse.add_argument('--cnv_profile',type=str,default=default,
        help='the file to save CNVs profile [{}]'.format(default))
    parse.add_argument('--snv_genotype',type=str,
        help='the file to save SNV genotypes for each cell')
    parse.add_argument('--ind_cnvs',type=str,
        help='the file to save CNVs for each cell individual')
    parse.add_argument('--parental_copy',type=str,
        help='the file to save parental copy for each SNV')
    parse.add_argument('--trunk_vars',type=str,
        help='the trunk variants file supplied by user')
    default=0
    parse.add_argument('--trunk_length',type=float,
        help='the length of the truncal branch [{}]'.format(default))
    default=None
    parse.add_argument('--expands',type=str,default=default,
        help='the basename of the file to output the snv and segment data for EXPANDS [{}]'.format(default))
    default=2.0
    parse.add_argument('--tstv',type=check_tstv,default=default,
        help='the ratio of ts/tv of SNV [{}]'.format(default))
    default=100000000
    parse.add_argument('--length',type=int,default=default,
        help='the length of the sequence to simulate [{}]'.format(default))
    default=None
    parse.add_argument('--chain',type=check_folder,default=default,
        help='directory to output chain files for each sample [{}]'.format(default))
    default=None
    parse.add_argument('--config',type=str,default=default,
        help='configure file contains the setting of the somatic simulation in YAML format [{}]'.format(default))
    args=parse.parse_args()

###### check config file
#only those params in cfg_params are acceptable in configure file
    cfg_params=('snv_rate','cnv_rate','del_prob','cnv_length_beta','cnv_length_max',
        'copy_parameter','copy_max','parental','tstv','length')
    config={}
    if args.config:
        with open(args.config,'r') as configfile:
            config=yaml.safe_load(configfile)
        check_config_file(config=config,cfg_params=cfg_params)

###### figure out the simulation setting for each chroms
#1. The setting in configure YAML file will override the setting in command line.
#2. In the configure file, the setting for individual chr will override the setting of genome.
    genome_cfg={}
    if 'genome' in config and config['genome']:
        for parameter in cfg_params:
            genome_cfg[parameter]=config['genome'].get(parameter,getattr(args,parameter))
    else:
        for parameter in cfg_params:
            genome_cfg[parameter]=getattr(args,parameter)

    final_chroms_cfg={} 
    max_ploidy=0
    if 'chromosomes' in config and config['chromosomes']:
        for i in config['chromosomes']:
            for chroms,chroms_cfg in i.items():
                final_chroms_cfg[chroms]={}
                for parameter in cfg_params:
                    final_chroms_cfg[chroms][parameter]=chroms_cfg.get(parameter,genome_cfg[parameter])
                if max_ploidy<len(final_chroms_cfg[chroms]['parental']):
                    max_ploidy=len(final_chroms_cfg[chroms]['parental'])
    else:
        final_chroms_cfg[args.name]={}
        for parameter in cfg_params:
            final_chroms_cfg[args.name][parameter]=genome_cfg[parameter]
        max_ploidy=len(final_chroms_cfg[args.name]['parental'])

###### logging and random seed setting
    logging.basicConfig(filename=args.log, filemode='w',
        format='[%(asctime)s] %(levelname)s: %(message)s',
        datefmt='%m-%d %H:%M:%S',level=args.loglevel)
    logging.info(' Command: %s',' '.join(sys.argv))
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
#TODO: We should do make sure the newick string is valided before processing it.
    mytree=csite.tree.newick2tree(newick)

    if args.trunk_length:
        mytree.lens=args.trunk_length
    leaves_number=mytree.leaves_counting()
    leaves_names=sorted(mytree.leaves_naming())

####### prune the tree if required
    if args.prune>0:
        if args.prune>leaves_number:
            raise argparse.ArgumentTypeError("There are only {} leaves on the tree. It's impossible to prune {} leaves.".format(
                leaves_number,args.prune))
        if args.prune>=2:
            mytree.prune(tips=args.prune)
    elif args.prune_proportion>0.0:
        trim=leaves_number*args.prune_proportion
        if trim>=2:
            mytree.prune(tips=trim)

###### output the map of tip_node:leaf
    if args.chain:
        tip_leaves=mytree.tip_node_leaves()
        os.mkdir(args.chain,mode=0o755)
        with open(args.chain+'/tip_node_sample.map','w') as tip_leaves_f:
            tip_leaves_f.write('#tip_node\tsample\n')
            for tip_node in sorted(tip_leaves.keys()):
                for leaf in sorted(tip_leaves[tip_node]):
                    tip_leaves_f.write('node{}\t{}\n'.format(tip_node,leaf))
        with open(args.chain+'/tip_node_sample.count','w') as tip_leaves_count_f:
            tip_leaves_count_f.write('#tip_node\tsample_count\n')
            for tip_node in sorted(tip_leaves.keys()):
                tip_leaves_count_f.write('node{}\t{}\n'.format(tip_node,len(tip_leaves[tip_node])))

###### add trunk vars if supplied
#FIXME: Can not use args.length here!!!!!!!!!!!
#FIXME: Can not use args.ploidy here!!!!!!!!!!!
    trunk_snvs={}
    trunk_cnvs={}
    if args.trunk_vars!=None:
        trunk_snvs,trunk_cnvs=csite.trunk_vars.classify_vars(
            args.trunk_vars,args.ploidy,args.length,leaves_number,mytree)

#TODO: check the output of nodes_snvs/nodes_vars and parental_copy
###### open all required output file and output the headers 
#TODO: some file's headers are missing
    cnv_file=open(args.cnv,'w')
    snv_file=open(args.snv,'w')
    snv_file.write('#position\ttrue_freq\ttotal_depth\tsimulated_freq\n')
    cnv_profile_file=open(args.cnv_profile,'w')
    nodes_snvs_file=open(args.nodes_vars,'w')

    if args.snv_genotype!=None:
        genotype_file=open(args.snv_genotype,'w')
        genotype_file.write('{}\t{}\n'.format('#positon','\t'.join(leaves_names)))
    
    if args.ind_cnvs!=None:
        ind_cnvs_file=open(args.ind_cnvs,'w')
        ind_cnvs_file.write('#cell\tparental\tstart\tend\tcopy\n')

#TODO: check the format of this file. parental? haplotype?
    if args.parental_copy!=None:
        parental_copy_file=open(args.parental_copy,'w')
        parental_copy_file.write('#position\t{}\n'.format('\t'.join(['haplotype'+str(x) for x in range(max_ploidy)])))

    if args.expands != None:
        expands_snps_file=open(args.expands+'.snps','w')
        expands_snps_file.write('chr\tstartpos\tAF_Tumor\tPN_B\n')
        expands_segs_file=open(args.expands+'.segs','w')
        expands_segs_file.write("chr\tstartpos\tendpos\tCN_Estimate\n")

###### simulate variants for each chroms
    for chroms,chroms_cfg in final_chroms_cfg.items():
        check_cnv_length_cfg(chroms=chroms,cnv_length_beta=chroms_cfg['cnv_length_beta'],
            cnv_length_max=chroms_cfg['cnv_length_max'],chr_length=chroms_cfg['length'])
        cn_dist_cfg=cn_dist(copy_max=chroms_cfg['copy_max'],copy_parameter=chroms_cfg['copy_parameter'])
        tstv_dist_cfg=tstv_dist(tstv=chroms_cfg['tstv'])
        
        (snvs_freq,cnvs,cnv_profile,nodes_snvs,tree_with_snvs,
            leaf_snv_alts,leaf_snv_refs,leaf_cnvs,
            hap_local_copy_for_all_snvs,
            )=mytree.snvs_freq_cnvs_profile(
                parental=chroms_cfg['parental'],
                snv_rate=chroms_cfg['snv_rate'],
                cnv_rate=chroms_cfg['cnv_rate'],
                del_prob=chroms_cfg['del_prob'],
                cnv_length_beta=chroms_cfg['cnv_length_beta'],
                cnv_length_max=chroms_cfg['cnv_length_max'],
                cn_dist_cfg=cn_dist_cfg,
                tstv_dist_cfg=tstv_dist_cfg,
                trunk_snvs=trunk_snvs,
                trunk_cnvs=trunk_cnvs,
                purity=args.purity,
                length=chroms_cfg['length'],
                chain=args.chain,
                chroms=chroms,
            )

        if args.snv_genotype!=None:
            for snv in snvs_freq:
                genotype_file.write('{}\t{}\n'.format(snv[0],
                    '\t'.join([str(leaf_snv_alts[leaf][snv[0]])+':'+str(leaf_snv_refs[leaf][snv[0]]) for leaf in leaves_names])))

        if args.ind_cnvs!=None:
            for leaf in sorted(leaf_cnvs.keys()):
                for cnv in leaf_cnvs[leaf]:
                    ind_cnvs_file.write('{}\n'.format('\t'.join([str(x) for x in [leaf,cnv['parental'],cnv['start'],cnv['end'],cnv['copy']]])))

        if args.parental_copy!=None:
            for snv in hap_local_copy_for_all_snvs:
                parental_copy_file.write('\t'.join([str(x) for x in snv])+'\n')

        for cnv in cnvs:
            cnv_file.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(cnv['seg'],cnv['start'],cnv['end'],cnv['copy'],cnv['leaves_count'],cnv['pre_snvs']))

        for pos,freq in snvs_freq:
            total_dp,b_allele_dp=csite.tree.simulate_sequence_coverage(args.depth,freq)
            b_allele_freq=0
            if total_dp!=0:
                b_allele_freq=b_allele_dp/total_dp
            snv_file.write('{}\t{}\t{}\t{}\n'.format(pos,freq,total_dp,b_allele_freq))

        for seg in cnv_profile:
            cnv_profile_file.write('{}\t{}\t{}\n'.format(*seg))

#FIXME: output SNVs/CNVs, not only SNVs
        for node in sorted(nodes_snvs.keys()):
            for snv in sorted(nodes_snvs[node]):
                nodes_snvs_file.write('{}\t{}\n'.format(node,snv))

#output for expands
        if args.expands != None:
            for pos,freq in snvs_freq:
                total_dp,b_allele_dp=csite.tree.simulate_sequence_coverage(args.depth,freq)
                expands_snps_file.write('{}\t{}\t{}\t{}\n'.format(chroms,pos,b_allele_dp/total_dp,0))

#in the segment input for expands
#CN_Estimate - the copy number estimated for each segment (average value across all subpopulations in the sample)
            for start,end,copy in cnv_profile:
                expands_segs_file.write('{}\t{}\t{}\t{}\n'.format(chroms,start,end,copy/leaves_number))

#TODO: Should we change pickle to json or yaml?
#http://stackoverflow.com/questions/4677012/python-cant-pickle-type-x-attribute-lookup-failed
        #with open(args.named_tree,'wb') as tree_data_file:
        #    pickle.dump(tree_with_snvs,tree_data_file)

###### close all opened files
    cnv_file.close()
    snv_file.close()
    cnv_profile_file.close()
    nodes_snvs_file.close()

    if args.snv_genotype!=None:
        genotype_file.close()
    
    if args.ind_cnvs!=None:
        ind_cnvs_file.close()

    if args.parental_copy!=None:
        parental_copy_file.close()

    if args.expands != None:
        expands_snps_file.close()
        expands_segs_file.close()
