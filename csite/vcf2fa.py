#!/usr/bin/env python3

#########################################################################
# Author: Hechuan Yang
# Created Time: 2017-09-27 10:50:14
# File Name: vcf2fa.py
# Description: 
#########################################################################

import os
import sys
import argparse
import logging
import pyfaidx

#handle the error below
#python | head == IOError: [Errno 32] Broken pipe 
from signal import signal, SIGPIPE, SIG_DFL 
signal(SIGPIPE,SIG_DFL) 

def main(progname=None):
    parse=argparse.ArgumentParser(
        description='Generate perturbed genome by integrating SNPs from a VCF.',
        prog=progname if progname else sys.argv[0])
    parse.add_argument('-r','--reference',required=True,
        help='a fasta file contains the reference genome')
    default=None
    parse.add_argument('-v','--vcf',type=str,default=default,
        help='a VCF file contains germline SNPs [{}]'.format(default))
    default='normal'
    parse.add_argument('-o','--output',type=str,default=default,
        help='output directory [{}]'.format(default))
    default=2
    parse.add_argument('-H','--haplotype',type=int,default=default,choices=[1,2],
        help='number of haplotypes to generate [{}]'.format(default))
    args=parse.parse_args()

    reference=pyfaidx.Fasta(args.reference)
    genome_profile=fai_info(fai=args.reference+'.fai')
    add_vcf_vars(profile=genome_profile,vcf=args.vcf,haplotype=args.haplotype)

    os.mkdir(args.output) 
    for i in range(args.haplotype):
        with open('{}/normal_hap{}.fa'.format(args.output,i),'w') as output:
            for chroms in genome_profile['order']:
                start=0
                segments=[]
                for snp in genome_profile[chroms]['hap_vars'][i]:
                    segments.append(reference[chroms][start:(snp[0]-1)].seq)
                    segments.append(snp[1])
                    start=snp[0]
                if start<=genome_profile[chroms]['length']-1:
                    segments.append(reference[chroms][start:].seq)
                output.write('>{}\n'.format(chroms,i))
                for outputline in pyfaidx.wrap_sequence(genome_profile[chroms]['linebases'],''.join(segments)):
                    output.write(outputline)

def fai_info(fai=None):
    '''
    Extract fasta information from genome.fa.fai file.
    Will return a list with the structure:
    {'order':[chroms1,chrom2,...],
     chroms1:{length,linebases,hap_vars:[[],[]]},
     chroms2:{length,linebases,hap_vars:[[],[]]},
     ...
    }
    '''
    profile={'order':[]}
    with open(fai,'r') as fai_file:
        for line in fai_file:
            line=line.rstrip()
            chroms,length,linebases=[line.split('\t')[x] for x in [0,1,3]]
            length=int(length)
            linebases=int(linebases)
            profile['order'].append(chroms)
            profile[chroms]={'length':length,
                             'linebases':linebases,
                             'hap_vars':[[],[]]}
    return profile

def add_vcf_vars(profile=None,vcf=None,haplotype=None):
    '''
    Extract variants on each copy of each chromosome in vcf file.
    And fill in the list hap_vars in profile.
    '''
    with open(vcf,'r') as vcf_file:
        for line in vcf_file:
            line=line.strip()
            if line.startswith('#CHROM'):
                if len(line.split())!=10:
                    raise VcfInputError('There should be ONE sample in the VCF file')
            elif not line.startswith('##'):
                chroms,pos,ref,alt,tags,values=[line.split('\t')[x] for x in [0,1,3,4,8,9]]
                pos=int(pos)
                alleles=[ref]
                alleles.extend(alt.split(','))
                tags_list=tags.split(':')
                values_list=values.split(':')
                indiv_info={}
                for i in range(len(tags_list)):
                    indiv_info[tags_list[i]]=values_list[i]
                if 'GT' not in indiv_info:
                    raise VcfInputError('Can not find GT information for {}:{}'.format(chroms,pos))
                if '|' not in indiv_info['GT']:
                    raise VcfInputError('Not phased on position {}:{}'.format(chroms,pos))
                gt=indiv_info['GT'].split('|')
                gt=[int(x) for x in gt]
                for i in range(haplotype):
                    if gt[i]!=0:
                        profile[chroms]['hap_vars'][i].append([pos,alleles[gt[i]]])

if __name__=='__main__':
    main()

