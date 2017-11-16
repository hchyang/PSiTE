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
#Let's fix at one diploid sample in VCF, which means haplotype==2
haplotype=2

def check_sex(chrs=None):
    if len(chrs)==0:
        sex_chr=[]
    else:
        sex_chr=chrs.split(',')
        if len(sex_chr)!=haplotype:
            raise argparse.ArgumentTypeError('{} is an invalid value for --sex_chr.\n'.format(chrs)+
                'Please specify two sex chromosomes. If there are two copies of the same chromosome, \n'+
                'just write it twice and seprate them by a comma! e.g. --sex_chr X,X \n')
    return sex_chr

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
    default=None
    parse.add_argument('--sex_chr',type=check_sex,default=default,
        help='sex chromosomes of the genome (seperated by comma) [{}]'.format(default))
    args=parse.parse_args()
    if args.sex_chr==None:
        args.sex_chr=[]

#build the data structure: genome_profile
    reference=pyfaidx.Fasta(args.reference)
    genome_profile=fai_info(fai=args.reference+'.fai',sex_chr=args.sex_chr)
#fill in the list hap_vars in genome_profile
    add_vcf_vars(profile=genome_profile,vcf=args.vcf,sex_chr=args.sex_chr)

    try:
        os.mkdir(args.output) 
    except FileExistsError:
        exit('Folder {} exists. Delete it or try another folder.'.format(args.output))
    except FileNotFoundError:
        exit("Can't create folder {}. Please create its parent directories first.".format(args.output))

    for i in range(haplotype):
        with open('{}/normal_hap{}.fa'.format(args.output,i),'w') as output:
            for chroms in genome_profile['order']:
                if i<len(genome_profile[chroms]['hap_vars']):
                    start=0
                    segments=[]
                    for snp in genome_profile[chroms]['hap_vars'][i]:
                        try:
                            segments.append(reference[chroms][start:(snp[0]-1)].seq)
                        except ValueError:
                            if snp[0]-start==1:
#This snp and the previous one is adjacent snps. e.g. chr1 45 (previous) and chr1 46 (current)
#If you retrive by reference['chr1'][45:(46-1)].seq, the return is not '', an error will pop actually.
                                pass
                            else:
                                raise
                        segments.append(snp[1])
                        start=snp[0]
                    if start<genome_profile[chroms]['length']:
                        segments.append(reference[chroms][start:].seq)
                    output.write('>{}\n'.format(chroms,i))
                    for outputline in pyfaidx.wrap_sequence(genome_profile[chroms]['linebases'],''.join(segments)):
                        output.write(outputline)

def fai_info(fai=None,sex_chr=None):
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
        if sex_chr:
            for chroms in sex_chr:
                assert chroms in profile, 'Can not find {} in your fasta file!'.format(chroms)
#There are two different sex chromosomes. So each of them should only have one 
#haplotype
            if sex_chr[0]!=sex_chr[1]:
                profile[sex_chr[0]]['hap_vars']=[[]]
                profile[sex_chr[1]]['hap_vars']=[[]]
    return profile

def add_vcf_vars(profile=None,vcf=None,sex_chr=None):
    '''
    Extract variants on each copy of each chromosome in vcf file.
    And fill in the list hap_vars in profile.
    '''
    if vcf.endswith('vcf.gz'):
        vcf_file=gzip.open(vcf,'rb')
    elif vcf.endswith('vcf'):
        vcf_file=open(vcf,'r')
    else:
        exit('For --vcf, only vcf/vcf.gz file are acceptable!')
    for line in vcf_file:
        if isinstance(line,bytes):
            line=line.decode('utf-8')
        line=line.strip()
        if line.startswith('#'):
            if line.startswith('#CHROM') and len(line.split())!=10:
                raise VcfInputError('Only ONE sample in VCF is acceptable.')
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
                    raise VcfInputError('Only SNPs are acceptable! Check the record below:\n{}\n'.format(line))
            if chroms in sex_chr and sex_chr[0]!=sex_chr[1]:
                if len(alt)==1:
                    profile[chroms]['hap_vars'][0].append([pos,alt])
                else:
                    raise VcfInputError('There is only one copy of chromosome: {},But multiple alternative alleles found in the record below:\n{}\n'.format(chroms,line))
            else:
                tags=record[8]
                values=record[9]
                tags_list=tags.split(':')
                values_list=values.split(':')
                indiv_info={}
                for i in range(len(tags_list)):
                    indiv_info[tags_list[i]]=values_list[i]
                if 'GT' not in indiv_info:
                    raise VcfInputError('Can not find GT information in the record below:\n{}\n'.format(line))
                if '|' not in indiv_info['GT']:
                    raise VcfInputError('Not phased genotype in record below:\n{}\n'.format(line))
                gt=indiv_info['GT'].split('|')
                gt=[int(x) for x in gt]
                for i in range(haplotype):
                    if gt[i]!=0:
                        profile[chroms]['hap_vars'][i].append([pos,alleles[gt[i]]])
    vcf_file.close()

class VcfInputError(Exception):
    pass

class FolderExistsError(Exception):
    pass

class ParentNotFoundError(Exception):
    pass

if __name__=='__main__':
    main()

