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
import re
import gzip
import pyfaidx
import time

#handle the error below
#python | head == IOError: [Errno 32] Broken pipe
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL)

nucleotide_re=re.compile('^[atcgnATCGN]$')


def check_vcf(vcf=None):
    if not vcf.endswith(('vcf','VCF','vcf.gz','VCF.gz','VCF.GZ','vcf.GZ')):
        raise argparse.ArgumentTypeError('For --vcf, only vcf/vcf.gz file are acceptable!')
    return vcf

def check_autosomes(autosomes_str=None):
    tmp=autosomes_str.split(',')
    for i in tmp:
        n=i.count('..')
        if n==0:
            pass
        elif n==1:
            m=re.search('^(.*?)([0-9]+)\.\.(\\1)?([0-9]+)$',i)
            if m:
                start=int(m.group(2))
                end=int(m.group(4))
                if start>=end:
                    raise argparse.ArgumentTypeError("The string '{}' is not valid in --autosomes".format(i))
            else:
                raise argparse.ArgumentTypeError("The string '{}' is not valid in --autosomes".format(i))
        else:
            raise argparse.ArgumentTypeError("The string '{}' is not valid in --autosomes".format(i))
    return autosomes_str

def check_sex(chrs=None):
    sex_chr=chrs.split(',')
    if len(sex_chr)!=2:
        raise argparse.ArgumentTypeError("'{}' is an invalid value for --sex_chr.\n".format(chrs)+
            'Please specify two sex chromosomes. If there are two copies of the same chromosome, \n'+
            'just write it twice and seprate them by a comma! e.g. --sex_chr X,X \n')
    return chrs

def check_output_folder(directory=None):
    good_charactors=re.compile('^[0-9a-zA-Z/_\-.]+$')
    if not good_charactors.match(directory):
        raise argparse.ArgumentTypeError("'{}' is an invalid string for --output. ".format(directory)+
            "Please only the combination of numbers, alphabets and ._/- as the directory name.")
    if os.path.exists(directory):
        raise argparse.ArgumentTypeError("'{}' exists already. Delete it or use another name instead.".format(directory))
    return directory


def main(progname=None):
    t0 = time.time()
    prog=progname if progname else sys.argv[0]
    parser=argparse.ArgumentParser(
        description='Build normal genome by integrating germline SNPs from a VCF file.',
        prog=prog)
    parser.add_argument('-v','--vcf',type=check_vcf,required=True,metavar='FILE',
        help='a VCF file containing germline SNPs')
    parser.add_argument('-r','--reference',type=str,required=True,metavar='FILE',
        help='a fasta file of reference genome')
    default='normal_fa'
    parser.add_argument('-o','--output',type=check_output_folder,default=default,metavar='DIR',
        help='output directory [{}]'.format(default))
    parser.add_argument('-a','--autosomes',type=check_autosomes,required=True,metavar='STR',
        help='autosomes of the genome (e.g. 1,2,3,4,5 or 1..4,5)')
    default=None
    parser.add_argument('-s','--sex_chr',type=check_sex,default=default,metavar='STR',
        help='sex chromosomes of the genome (separated by comma) [{}]'.format(default))
    args=parser.parse_args()
    if args.sex_chr==None:
        args.sex_chr=[]
    else:
        args.sex_chr=args.sex_chr.split(',')
    autosomes=parse_autosomes(args.autosomes)

#build the data structure: genome_profile
    reference=pyfaidx.Fasta(args.reference)
    genome_profile=fai_info(fai=args.reference+'.fai',autosomes=autosomes,sex_chr=args.sex_chr)
#fill in the list hap_vars in genome_profile
    add_vcf_vars(profile=genome_profile,vcf=args.vcf)

    os.mkdir(args.output,mode=0o755)

    for i in range(2):
        with open(os.path.join(args.output,'normal.parental_{}.fa'.format(i)),'w') as output:
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
                    output.write('>{}\n'.format(chroms))
                    for outputline in pyfaidx.wrap_sequence(genome_profile[chroms]['linebases'],''.join(segments)):
                        output.write(outputline)
    t1 = time.time()
    print ("Total time running {}: {} seconds".format
       (prog, str(t1-t0)))


def parse_autosomes(autosomes_str=None):
    '''
    Parse the --autosomes string into a set of individual chromsomes.
    e.g. chr1..3,chr8,chr11..chr13,14..16
         => {chr1,chr2,chr3,chr8,chr11,chr12,chr13,14,15,16}
    '''
    tmp=autosomes_str.split(',')
    autosomes=[]
    for i in tmp:
        n=i.count('..')
        if n==0:
            autosomes.append(i)
        elif n==1:
            m=re.search('^(.*?)([0-9]+)\.\.(\\1)?([0-9]+)$',i)
            if m:
                prefix=m.group(1)
                start=int(m.group(2))
                end=int(m.group(4))
                autosomes.extend([prefix+str(x) for x in range(start,end+1)])
            else:
                raise argparse.ArgumentTypeError("The string '{}' is not valid in --autosomes".format(i))
        else:
            raise argparse.ArgumentTypeError("The string '{}' is not valid in --autosomes".format(i))
    return set(autosomes)

def fai_info(fai=None,autosomes=None,sex_chr=None):
    '''
    Extract fasta information from the file genome.fa.fai.
    Will return a list with the structure:
    {'order':[chroms1,chrom2,...],
     chroms1:{length,linebases,hap_vars:[[],[]]},
     chroms2:{length,linebases,hap_vars:[[],[]]},
     ...
    }
    '''
    profile={'order':[]}
    want=autosomes.union(sex_chr)
    with open(fai,'r') as fai_file:
        for line in fai_file:
            field=line.rstrip().split('\t')
            chroms=field[0]
            length=int(field[1])
            linebases=int(field[3])
            if chroms in want:
                profile['order'].append(chroms)
                profile[chroms]={'length':length,
                                 'linebases':linebases,
                                 'hap_vars':[[],[]]}
                if chroms in sex_chr and sex_chr[0]!=sex_chr[1]:
                    profile[chroms]['hap_vars']=[[]]
    not_found=want-set(profile['order'])
    if not_found:
        raise ChrNotFoundError("Couldn't find chromosome '{}' in the reference file!".format(not_found))
    return profile

def add_vcf_vars(profile=None,vcf=None):
    '''
    Extract variants on each copy of each chromosome in vcf file.
    And fill in the list hap_vars in profile.
    '''
    gz=False
    if vcf.endswith(('gz','GZ')):
        vcf_file=gzip.open(vcf,'rb')
        gz=True
    else:
        vcf_file=open(vcf,'r')
    for line in vcf_file:
        if gz:
            line=line.decode('utf-8')
        line=line.strip()
        if line.startswith('#'):
            if line.startswith('#CHROM') and len(line.split())!=10:
                raise VcfInputError('Only VCF containing ONE sample is acceptable.')
        else:
            field=line.split('\t')
            chroms=field[0]
            if chroms in profile['order']:
                pos=int(field[1])
                ref=field[3]
                alt=field[4]
                alleles=[ref]+alt.split(',')
                for n in alleles:
                    if not nucleotide_re.match(n):
                        raise VcfInputError('Only SNPs are acceptable! Check the record below:\n{}\n'.format(line))
                if len(profile[chroms]['hap_vars'])==1:
                    if len(alt)==1:
                        profile[chroms]['hap_vars'][0].append([pos,alt])
                    else:
                        raise VcfInputError('There is only one copy of chromosome: {}, '.format(chroms)+
                            'but multiple alternative alleles found in the record below:\n{}\n'.format(line))
                else:
                    tags=field[8]
                    values=field[9]
                    tags_list=tags.split(':')
                    values_list=values.split(':')
                    indiv_info={}
                    for i in range(len(tags_list)):
                        indiv_info[tags_list[i]]=values_list[i]
                    if 'GT' not in indiv_info:
                        raise VcfInputError("Couldn't find GT information in the record below:\n{}\n".format(line))
                    if '|' not in indiv_info['GT']:
                        raise VcfInputError('Not phased genotype in record below:\n{}\n'.format(line))
                    gt=indiv_info['GT'].split('|')
                    gt=[int(x) for x in gt]
                    for i in range(2):
                        if gt[i]!=0:
                            profile[chroms]['hap_vars'][i].append([pos,alleles[gt[i]]])
    vcf_file.close()

class ChrNotFoundError(Exception):
    pass

class VcfInputError(Exception):
    pass

class FolderExistsError(Exception):
    pass

class ParentNotFoundError(Exception):
    pass

if __name__=='__main__':
    main()
