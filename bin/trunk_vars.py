#!/usr/bin/env python

#########################################################################
# Author: Hechuan Yang
# Created Time: 2017-04-20 10:10:59
# File Name: trunk_vars.py
# Description: 
#########################################################################

import logging
import copy as cp

def classify_vars(vars_file,ploid,seq_length,leaves_number,tree):
    '''
    Input vars should be in this format
    chr start end indicator
    chr:         which copy of chrs the var locates in (0 based)
    start:       the start of the var
    end:         the end of the var
    indicator:   an interger. 0: SNV, -1: deletion, +int: amplification
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
            var=line.split()
            chrom=int(var[0])
            start=int(var[1])
            end=int(var[2])
            copy=int(var[3])
            if not 0<=start<seq_length or not 0<=end<seq_length: 
                raise VarCoordinateError

            if copy==0:
                if not chrom in snvs:
                    snvs[chrom]=[]
                snvs[chrom].append(start)
            else:
                if not chrom in cnvs:
                    cnvs[chrom]=[]
#construct cnv
                cnvs[chrom].append({'seg': [0,seq_length],
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

#right now, copy must be -1 or a positive integer
                for i in range(copy):
                    segment=cp.deepcopy(tree)
                    cnvs[chrom][-1]['new_copies'].append(segment)

    check_vars(snvs,amps,dels,ploid)
    
    logging.debug('trunk SNVs:%s',snvs)
    logging.debug('trunk AMPs:%s',amps)
    logging.debug('trunk DELs:%s',dels)
    logging.debug('trunk CNVs:%s',cnvs)
    return snvs,dels,cnvs

def check_vars(snvs,amps,dels,ploid):
    '''
    Check: 1. whether all variants' ploid is valid.
           2. whether any snv/amp overlap with deletion.
    '''
    invalid=set(snvs).union(set(amps)).union(set(dels))-set(range(ploid))
    if invalid:
#there are variants on chr x, which is not in [0,ploid). 
#P.S. the index of chr is 0 based.
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

class VarCoordinateError(Exception):
    pass

class VarInDelError(Exception):
    pass

class PloidyError(Exception):
    pass

class CnvCopyError(Exception):
    pass
