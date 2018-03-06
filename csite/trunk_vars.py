#!/usr/bin/env python3

#########################################################################
# Author: Hechuan Yang
# Created Time: 2017-04-20 10:10:59
# File Name: trunk_vars.py
# Description: 
#########################################################################

import logging
import copy as cp

def classify_vars(vars_file,chroms_cfg,leaves_number,tree):
    '''
    There should be 4 columns for each varians in the input file,
    chr:      the chromosome of the var
    hap:      which halotype of the chr the var locates in (0 based)
    start:    the start of the var
    end:      the end of the var
    var:      an string. 0/1/2: SNV, -1: deletion, +int: amplification
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
            chroms,hap,start,end,var=line.split()
            hap=int(hap)
            start=int(start)
            end=int(end)
            if chroms not in chroms_cfg['order']:
                raise TrunkVarError('The chr of the variant below is not in the genome:\n{}'.format(line))
            if not 0<=hap<len(chroms_cfg[chroms]['parental']):
                raise TrunkVarError('The haplotype of the variant below is out of range:\n{}'.format(line))
            if not 0<=start<chroms_cfg[chroms]['length'] or not 0<=end<chroms_cfg[chroms]['length']: 
                raise TrunkVarError('The coordinant of the variant below is out of range:\n{}'.format(line))

            if var.startswith('+') or var.startswith('-'):
                copy=int(var)
                if chroms not in cnvs:
                    cnvs[chroms]={}
                if hap not in cnvs[chroms]:
                    cnvs[chroms][hap]=[]
#construct cnv
                cnvs[chroms][hap].append({'seg': [0,chroms_cfg[chroms]['length']],
                                         'start': start,
                                         'end': end,
                                         'copy': copy,
                                         'leaves_count': leaves_number,
                                         'pre_snvs': [],
                                         'new_copies': [],
                                        })

                if copy==-1:
                    cnvs[chroms][hap][-1]['type']='DEL'
                    if chroms not in dels:
                        dels[chroms]={}
                    if hap not in dels[chroms]:
                        dels[chroms][hap]=[]
                    dels[chroms][hap].append([start,end,copy])
                elif copy>0:
                    cnvs[chroms][hap][-1]['type']='AMP'
                    if chroms not in amps:
                        amps[chroms]={}
                    if hap not in amps[chroms]:
                        amps[chroms][hap]=[]
                    amps[chroms][hap].append([start,end,copy])
                    for i in range(copy):
                        segment=cp.deepcopy(tree)
                        cnvs[chroms][hap][-1]['new_copies'].append(segment)
                else:
#right now, copy must be -1 or a positive integer
                    raise TrunkVarError('The fourth column of the variant below is invalid:\n{}'.format(line))
            else:
                if end-start!=1:
                    raise TrunkVarError('The coordinant of the SNV below is not correct:\n{}'.format(line))
                if var not in ('0','1','2'):
                    raise TrunkVarError('The mutation form of the SNV below is not correct:\n{}'.format(line))
                form=int(var)
                if chroms not in snvs:
                    snvs[chroms]={}
                if hap not in snvs[chroms]:
                    snvs[chroms][hap]=[]
                snvs[chroms][hap].append({'type':'SNV',
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
    for chroms in sorted(dels.keys()):
        snvs_chroms=snvs.get(chroms,{})
        amps_chroms=amps.get(chroms,{})
        for hap in sorted(dels[chroms].keys()):
            snvs_chrom_hap=snvs_chroms.get(hap,[])
            amps_chrom_hap=amps_chroms.get(hap,[])
            for deletion in dels[chroms][hap]:
                for snv in snvs_chrom_hap:
                    if deletion[0]<=snv['start']<deletion[1]:
                        raise TrunkVarError('Those variants below are conflict:\n'+
                            '{}\n'.format('\t'.join([str(x) for x in [chroms,hap,snv['start'],snv['end'],snv['mutation']]]))+
                            '{}\n'.format('\t'.join([str(x) for x in [chroms,hap,deletion[0],deletion[1],'-1']])))
                for amp in amps_chrom_hap:
                    if deletion[0]<=amp[0]<deletion[1] or deletion[0]<amp[1]<=deletion[1]:
                        raise TrunkVarError('Those variants below are conflict:\n'+
                            '{}\n'.format('\t'.join([str(x) for x in [chroms,hap,amp[0],amp[1],'+'+str(amp[2])]]))+
                            '{}\n'.format('\t'.join([str(x) for x in [chroms,hap,deletion[0],deletion[1],'-1']])))

class TrunkVarError(Exception):
    pass

