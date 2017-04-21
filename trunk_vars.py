#!/usr/bin/env python

#########################################################################
# Author: Hechuan Yang
# Created Time: 2017-04-20 10:10:59
# File Name: trunk_vars.py
# Description: 
#########################################################################

import logging

def classify_vars(vars_file,ploid,leaves_number):
    '''
    Input vars should be in this format
    chr start end indicator
    chr:   which copy of chrs the var locates in (start from 0)
    start: the start of the var
    end:   the end of the var
    int:   an interger. 0: SNV, -1: deletion, +int: amplification
    '''
    snvs={}
    amps={}
    dels={}
    cnvs={}

#classify vars into different categories
    with open(vars_file) as f:
        for line in f:
            var=line.split()
            chrom=int(var[0])
            start=float(var[1])
            end=float(var[2])
            copy=int(var[3])
            if copy==0:
                if not chrom in snvs:
                    snvs[chrom]=[]
                snvs[chrom].append(start)
            else:
                if not chrom in cnvs:
                    cnvs[chrom]=[]
                cnvs[chrom].append({'seg': [0,1],
                                     'start': start,
                                     'end': end,
                                     'copy': copy,
                                     'leaves_count': leaves_number,
                                     'pre_snvs': [],
                                     'new_copies': [],
                                    })

                if copy==-1:
                    if not chrom in dels:
                        dels[chrom]=[]
                    dels[chrom].append([start,end])
                elif isinstance(copy,int) and copy>0:
                    if not chrom in amps:
                        amps[chrom]=[]
                    amps[chrom].append([start,end])
                else:
                    raise CnvCopyError

    check_vars(snvs,amps,dels,ploid)
    
    logging.debug('trunk SNVs:%s',snvs)
    logging.debug('trunk AMPs:%s',amps)
    logging.debug('trunk DELs:%s',dels)
    logging.debug('trunk CNVs:%s',cnvs)
    return snvs,dels,cnvs

def check_vars(snvs,amps,dels,ploid):
    '''
    Check: 1. whether all var's ploid is valid
           2. whether any snv/amp overlap with deletion.
    '''
    invalid=set(snvs).union(set(amps)).union(set(dels))-set(range(ploid))
    if invalid:
        raise PloidyError
    for chrom in dels:
        for deletion in dels[chrom]:
            if chrom in snvs:
                for snv in snvs[chrom]:
                    if deletion[0]<=snv<deletion[1]:
                        raise VarInDelError
            if chrom in amps:
                for amp in amps[chrom]:
                    if deletion[0]<=amp[0]<deletion[1] or deletion[0]<amp[1]<=deletion[1]:
                        raise VarInDelError

class VarInDelError(Exception):
    pass

class PloidyError(Exception):
    pass

class CnvCopyError(Exception):
    pass
