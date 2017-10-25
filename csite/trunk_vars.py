#!/usr/bin/env python3

#########################################################################
# Author: Hechuan Yang
# Created Time: 2017-04-20 10:10:59
# File Name: trunk_vars.py
# Description: 
#########################################################################

import logging
import copy as cp

#def classify_vars(vars_file,ploid,seq_length,leaves_number,tree):
def classify_vars(vars_file,chroms_cfg,leaves_number,tree):
    '''
    There should be 4 columns for each varians in the input file,
    chr:      the chromosome of the var
    hap:      which halotype of the chr the var locates in (0 based)
    start:    the start of the var
    end:      the end of the var
    variant:  an string. 0/1/2: SNV, -1: deletion, +int: amplification
    P.S. start and end are 0 based. And the region of each var is like in bed: [start,end).
    '''
    snvs={}
    amps={}
    dels={}
    cnvs={}

#classify vars into different categories
    with open(vars_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            chrom,hap,start,end,var=line.split()
            hap=int(hap)
            start=int(start)
            end=int(end)
            if chrom not in chroms_cfg:
                raise TrunkVarError('The chr of the variant below is not in the genome:\n{}'.format(line))
            if not 0<=hap<len(chroms_cfg[chrom]['parental']):
                raise TrunkVarError('The haplotype of the variant below is out of range:\n{}'.format(line))
            if not 0<=start<chroms_cfg[chrom]['length'] or not 0<=end<chroms_cfg[chrom]['length']: 
                raise TrunkVarError('The coordinant of the variant below is out of range:\n{}'.format(line))

            if var.startswith('+') or var.startswith('-'):
                copy=int(var)
                if chrom not in cnvs:
                    cnvs[chrom]={}
                if hap not in cnvs[chrom]:
                    cnvs[chrom][hap]=[]
#construct cnv
                cnvs[chrom][hap].append({'seg': [0,chroms_cfg[chrom]['length']],
                                         'start': start,
                                         'end': end,
                                         'copy': copy,
                                         'leaves_count': leaves_number,
                                         'pre_snvs': [],
                                         'new_copies': [],
                                        })

                if copy==-1:
                    cnvs[chrom][hap][-1]['type']='DEL'
                    if chrom not in dels:
                        dels[chrom]={}
                    if hap not in dels[chrom]:
                        dels[chrom][hap]=[]
                    dels[chrom][hap].append([start,end,copy])
                elif copy>0:
                    cnvs[chrom][hap][-1]['type']='AMP'
                    if chrom not in amps:
                        amps[chrom]={}
                    if hap not in amps[chrom]:
                        amps[chrom][hap]=[]
                    amps[chrom][hap].append([start,end,copy])
                    for i in range(copy):
                        segment=cp.deepcopy(tree)
                        cnvs[chrom][hap][-1]['new_copies'].append(segment)
                else:
#right now, copy must be -1 or a positive integer
                    raise TrunkVarError('The fourth column of the variant below is invalid:\n{}'.format(line))
            else:
                if end-start!=1:
                    raise TrunkVarError('The coordinant of the SNV below is not correct:\n{}'.format(line))
                form=int(var)
                if chrom not in snvs:
                    snvs[chrom]={}
                if hap not in snvs[chrom]:
                    snvs[chrom][hap]=[]
                snvs[chrom][hap].append({'type':'SNV',
                                         'start':start,
                                         'end':end,
                                         'mutation':form,
                                        })

    check_vars(snvs,amps,dels)
    
    logging.debug('trunk SNVs:%s',snvs)
    logging.debug('trunk AMPs:%s',amps)
    logging.debug('trunk DELs:%s',dels)
    logging.debug('trunk CNVs:%s',cnvs)
    return snvs,cnvs

def check_vars(snvs,amps,dels):
    '''
    Check: whether any snv/amp overlap with deletion.
    '''
    for chrom in sorted(dels.keys()):
        snvs_chrom=snvs.get(chrom,{})
        amps_chrom=amps.get(chrom,{})
        for hap in sorted(dels[chrom].keys()):
            snvs_chrom_hap=snvs_chrom.get(hap,[])
            amps_chrom_hap=amps_chrom.get(hap,[])
            for deletion in dels[chrom][hap]:
                for snv in snvs_chrom_hap:
                    if deletion[0]<=snv['start']<deletion[1]:
                        raise TrunkVarError('Those variants below are conflict:\n'+
                            '{}\n'.format('\t'.join([str(x) for x in [chrom,hap,snv['start'],snv['end'],snv['mutation']]]))+
                            '{}\n'.format('\t'.join([str(x) for x in [chrom,hap,deletion[0],deletion[1],'-1']])))
                for amp in amps_chrom_hap:
                    if deletion[0]<=amp[0]<deletion[1] or deletion[0]<amp[1]<=deletion[1]:
                        raise TrunkVarError('Those variants below are conflict:\n'+
                            '{}\n'.format('\t'.join([str(x) for x in [chrom,hap,amp[0],amp[1],'+'+str(amp[2])]]))+
                            '{}\n'.format('\t'.join([str(x) for x in [chrom,hap,deletion[0],deletion[1],'-1']])))

class TrunkVarError(Exception):
    pass

