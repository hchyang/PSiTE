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
    There should be at least 5 columns for each varians in the input file,
    The 6th column is optional.
    chr:      the chromosome of the variants
    hap:      which haplotype of the chromosome the variant locates in (0 based).
              It depends on the parental of the chromosome. If the parental is '011',
              then there are 3 haplotypes: 0,1,2
    start:    the start of the variant
    end:      the end of the variant
    var:      an integer. 0/1/2: SNV, -1: deletion, +int: amplification
    focal_cp: optional. an integer or some integers separeted by comma (only for SNV). 
              0: the SNV is on the original copy
              N: the SNV is on the copy N of this segment
              Without this column: the SNV is on the original copy and all duplicated copy
    P.S. start and end are 0 based. And the region of each var is like in bed: [start,end).
    '''
    snvs={}
    cnvs={}

#classify vars into different categories
    with open(vars_file) as input:
        header=next(input)
        if not header.startswith('#'):
            raise TrunkVarError("The first line should be a header line that starts with '#'!")
        header=header.lstrip('#')
        header=header.rstrip()
        header=header.split()
        if header!=['chr','hap','start','end','var','focal_cp'] or \
            header!=['chr','hap','start','end','var']:
            raise TrunkVarError('The format of your trunk variants file is not right!')
        for line in input:
            cols=line.split()
            if len(cols)==5:
                chroms,hap,start,end,var=cols
                focal_cp=None
            elif len(cols)==6:
                chroms,hap,start,end,var,focal_cp=cols
                focal_cp=[int(x) for x in focal_cp.split(',')]
                if var not in ('0','1','2'):
                    raise TrunkVarError('Only SNVs can have the focal_cp (6th) column.'+
                        'Check the record below:\n{}'.format(line))
            else:
                raise TrunkVarError('There should be 5 or 6 columns in your --trunk_vars file.\n'+
                    'Check the record below:\n{}'.format(line))

            hap=int(hap)
            start=int(start)
            end=int(end)
            if chroms not in chroms_cfg['order']:
                raise TrunkVarError('The chr of the variant below is not in the genome:\n{}'.format(line))
            if not 0<=hap<len(chroms_cfg[chroms]['parental']):
                raise TrunkVarError('The haplotype of the variant below is out of range:\n{}'.format(line))
            if not (0<=start<chroms_cfg[chroms]['length'] and 0<=end<chroms_cfg[chroms]['length']): 
                raise TrunkVarError('The coordinant of the variant below is out of range:\n{}'.format(line))
            if not start<end: 
                raise TrunkVarError('The start of the variant should be less than its end:\n{}'.format(line))

            if var.startswith(('+','-')):
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
                                         'pre_snvs': {},
                                         'new_copies': [],
                                        })

                if copy==-1:
                    cnvs[chroms][hap][-1]['type']='DEL'
                    cnvs[chroms][hap][-1]['pre_snvs']={0:[]}
                elif copy>0:
                    cnvs[chroms][hap][-1]['type']='AMP'
                    for i in range(copy): 
                        segment=cp.deepcopy(tree)
                        cnvs[chroms][hap][-1]['new_copies'].append(segment)
                        cnvs[chroms][hap][-1]['pre_snvs'][i+1]=[]
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
                                         'focal_cp':focal_cp
                                        })

    snvs,cnvs=check_overlap(snvs,cnvs)
    logging.debug('trunk SNVs:%s',snvs)
    logging.debug('trunk CNVs:%s',cnvs)
    return snvs,cnvs

def check_overlap(snvs,cnvs):
    '''
    Check: 1) whether any CNV overlaps with another CNV (Error)
           2) SNV overlap a deletion (Error)
           3) SNV overlap an amplification (add it to the list of pre_snvs of that amplification)
    '''
#check the CNVs and the SNVs in CNVs
    for chroms in sorted(cnvs.keys()):
        snvs[chroms]=snvs.get(chroms,{})
        for hap in sorted(cnvs[chroms].keys()):
            snvs[chroms][hap]=snvs[chroms].get(hap,[])
            for i in range(len(cnvs[chroms][hap])):
                cnv1=cnvs[chroms][hap][i]
#compare a CNV with other CNVs
                if i<len(cnvs[chroms][hap])-1:
                    for cnv2 in cnvs[chroms][hap][i+1:]:
                        if cnv1['start']<=cnv2['start']<cnv1['end'] or cnv1['start']<cnv2['end']<=cnv1['end']:
                            raise TrunkVarError('These variants below are in conflict with each other:\n'+
                                '{}\n'.format('\t'.join([str(x) for x in [chroms,hap,cnv1['start'],cnv1['end'],str(cnv1['copy'])]]))+
                                '{}\n'.format('\t'.join([str(x) for x in [chroms,hap,cnv2['start'],cnv2['end'],str(cnv2['copy'])]])))
#compare a CNV with SNVs
                for snv in snvs[chroms][hap]:
                    if cnv1['start']<=snv['start']<cnv1['end']:
                        if cnv1['copy']==-1: #deletion
                            raise TrunkVarError('These variants below are in conflict with each other:\n'+
                                '{}\n'.format('\t'.join([str(x) for x in [chroms,hap,snv['start'],snv['end'],snv['mutation']]]))+
                                '{}\n'.format('\t'.join([str(x) for x in [chroms,hap,cnv1['start'],cnv1['end'],'-1']])))
                        else: #amplification
                            if snv['focal_cp']!=None:
                                if 0 in snv['focal_cp']:
                                    if snv['focal_cp']==[0]:
                                        del(snv['focal_cp'])
                                        continue
                                    else:
                                        for j in [x for x in snv['focal_cp'] if x!=0]:
                                            if j not in cnv1['pre_snvs']:
                                                cnv1['pre_snvs'][j]=[]
                                            cnv1['pre_snvs'][j].append(snv)
                                else:
                                    for j in snv['focal_cp']:
                                        if j not in cnv1['pre_snvs']:
                                            cnv1['pre_snvs'][j]=[]
                                        cnv1['pre_snvs'][j].append(snv)
                                    snv['not_on_original']=True
                            else:
                                for j in range(1,cnv1['copy']+1):
                                    if j not in cnv1['pre_snvs']:
                                        cnv1['pre_snvs'][j]=[]
                                    cnv1['pre_snvs'][j].append(snv)
                            del(snv['focal_cp'])
                snvs[chroms][hap]=[snv for snv in snvs[chroms][hap] if not snv.get('not_on_original',False)]
#check other SNVs
#check whether there are any snv with focal_cp information
#but without overlapping with any cnv
    for chroms in sorted(snvs.keys()):
        for hap in sorted(snvs[chroms].keys()):
            for snv in snvs[chroms][hap]:
                if 'focal_cp' in snv:
                    if snv['focal_cp']==None or snv['focal_cp']==[0]:
                        del(snv['focal_cp'])
                    else:
                        raise TrunkVarError('The SNV below is not covered by any CNV:\n'+
                            '{}\n'.format('\t'.join([str(x) for x in [chroms,hap,snv['start'],snv['end'],snv['mutation'],snv['focal_cp']]])))
    return snvs,cnvs

class TrunkVarError(Exception):
    pass

