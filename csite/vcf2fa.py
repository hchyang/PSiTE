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
import re
import gzip
import pyfaidx

#handle the error below
#python | head == IOError: [Errno 32] Broken pipe 
from signal import signal, SIGPIPE, SIG_DFL 
signal(SIGPIPE,SIG_DFL) 

nucleotide_re=re.compile('^[atcgnATCGN]$')

def main(progname=None):
    parse=argparse.ArgumentParser(
        description='Build normal genome by integrating germline SNPs from a VCF file.',
        prog=progname if progname else sys.argv[0])
    parse.add_argument('-v','--vcf',type=str,required=True,
        help='a VCF file contains germline SNPs')
    parse.add_argument('-r','--reference',type=str,required=True,
        help='a fasta file of reference genome')
    default='normal_fa'
    parse.add_argument('-o','--output',type=str,default=default,
        help='output directory [{}]'.format(default))
    default=2
    parse.add_argument('-H','--haplotype',type=int,default=default,choices=[1,2],
        help='number of haplotypes to generate [{}]'.format(default))
    args=parse.parse_args()

    reference=pyfaidx.Fasta(args.reference)
#build the data structure: genome_profile
    genome_profile=fai_info(fai=args.reference+'.fai',haplotype=args.haplotype)
#fill in the list hap_vars in genome_profile
    add_vcf_vars(profile=genome_profile,vcf=args.vcf,haplotype=args.haplotype)

    try:
        os.mkdir(args.output) 
    except FileExistsError:
        exit('Folder {} exists. Delete it or try another folder.'.format(args.output))
    except FileNotFoundError:
        exit("Can't create folder {}. Please creat parent directories first.".format(args.output))

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

def fai_info(fai=None,haplotype=None):
    '''
    Extract fasta information from genome.fa.fai file.
    Will return a list with the structure:
    {'order':[chroms1,chrom2,...],
     chroms1:{length,linebases,hap_vars:[[],...]},
     chroms2:{length,linebases,hap_vars:[[],...]},
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
                             'hap_vars':[]}
            for i in range(haplotype):
                profile[chroms]['hap_vars'].append([])

    return profile

def add_vcf_vars(profile=None,vcf=None,haplotype=None):
    '''
    Extract variants on each copy of each chromosome in vcf file.
    And fill in the list hap_vars in profile.
    '''
    if vcf.endswith('.gz'):
        vcf_file=gzip.open(vcf,'rb')
    else:
        vcf_file=open(vcf,'r')
    for line in vcf_file:
        if isinstance(line,bytes):
            line=line.decode('utf-8')
        line=line.strip()
        if line.startswith('#'):
            if line.startswith('#CHROM'):
                if len(line.split())!=10 and haplotype==2:
                    raise VcfInputError('When --haploype is 2, only ONE sample is acceptable.')
        else:
            record=line.split('\t')
            chroms=record[0]
            pos=int(record[1])
            ref=record[3]
            alt=record[4]
            alleles=[ref]
            alleles.extend(alt.split(','))
            for n in alleles:
                if not nucleotide_re.match(n):
                    print(line)
                    raise VcfInputError('Only SNPs are acceptable!')
            if haplotype==1:
                if len(alt)==1:
                    profile[chroms]['hap_vars'][0].append([pos,alt])
                elif len(alt)>1:
                    print(line)
                    raise VcfInputError('When --haplotype is 1, only one alternative allele is acceptable.')
                else:
                    raise ShouldNotBeHereError
            elif haplotype==2:
                tags=record[8]
                values=record[9]
                tags_list=tags.split(':')
                values_list=values.split(':')
                indiv_info={}
                for i in range(len(tags_list)):
                    indiv_info[tags_list[i]]=values_list[i]
                if 'GT' not in indiv_info:
                    print(line)
                    raise VcfInputError('Can not find GT information for {}:{}'.format(chroms,pos))
                if '|' not in indiv_info['GT']:
                    print(line)
                    raise VcfInputError('Not phased on position {}:{}'.format(chroms,pos))
                gt=indiv_info['GT'].split('|')
                gt=[int(x) for x in gt]
                for i in range(haplotype):
                    if gt[i]!=0:
                        profile[chroms]['hap_vars'][i].append([pos,alleles[gt[i]]])
            else:
                raise ShouldNotBeHereError
    vcf_file.close()

class VcfInputError(Exception):
    pass

class ShouldNotBeHereError(Exception):
    pass

class FolderExistsError(Exception):
    pass

class ParentNotFoundError(Exception):
    pass

if __name__=='__main__':
    main()

