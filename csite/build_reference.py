#!/usr/bin/env python3

#########################################################################
# Author: Hechuan
# Created Time: 2017-09-12 16:28:34
# File Name: build_reference.py
# Description: 
#########################################################################

import os
import sys
import re
import argparse
import pyfaidx
import glob
import numpy

#handle the error below
#python | head == IOError: [Errno 32] Broken pipe 
from signal import signal, SIGPIPE, SIG_DFL 
signal(SIGPIPE,SIG_DFL) 

def check_perturbed_folder(directory=None):
    if not os.path.isdir(directory):
        raise argparse.ArgumentTypeError("{} is not exist or not a folder.".format(directory))
    return directory
    
def check_reference_file(reference=None):
    if not os.path.isfile(reference):
        raise argparse.ArgumentTypeError("{} is not exist or not a file.".format(reference))
    return reference
    
def check_sequence_folder(directory=None):
    good_charactors=re.compile('^[0-9a-zA-Z_-]+$') 
    if not good_charactors.match(directory):
        raise argparse.ArgumentTypeError("{} is an invalid string for --sequence. ".format(directory)+
            "Please only use number, alphabet, - and _ in the directory name.")
    if os.path.exists(directory):
        raise argparse.ArgumentTypeError("{} is already exist. Delete it or use another name instead.".format(directory))
    return directory
    
def main():
    parse=argparse.ArgumentParser(description='Build reference genomes for ART to simulate short reads')
    parse.add_argument('-p','--perturbed',required=True,type=check_perturbed_folder,
        help='the folder contain the configure file of perturbed genomes')
    parse.add_argument('-r','--reference',required=True,type=check_reference_file,
        help='reference file')
    default='perturbed_seqs'
    parse.add_argument('-s','--sequence',default=default,type=check_sequence_folder,
        help='the folder to save the built genome sequences [{}]'.format(default))
    default=60
    parse.add_argument('-w','--width',default=default,type=int,
        help='the line width of output fasta [{}]'.format(default))
    args=parse.parse_args()

    reference=pyfaidx.Fasta(args.reference)
    os.mkdir(args.sequence,mode=0o755)
    for node_cfg in glob.glob(args.perturbed+'/node*.cfg'):
        with open(args.sequence+'/'+os.path.basename(node_cfg)+'.fa','w') as outputf:
            with open(node_cfg) as inputf:
                seq_name=''
                seq=[]
                for line in inputf:
                    line=line.rstrip()
                    if line.startswith('>'):
                        if seq:
                            outputf.write('>{}\n'.format(seq_name))
                            for outputline in pyfaidx.wrap_sequence(args.length,''.join(seq)):
                                outputf.write(outputline)
                        seq_name=line[1:]
                        seq=[]
                    else:
                        record=line.split()
                        chroms=record[0]
                        start=int(record[1])
                        end=int(record[2])
                        seq_type=record[3]
                        segment=''
                        if seq_type=='ref':
                            segment=reference[chroms][start:end].seq
                        elif seq_type=='SNV':
                            ref=reference[chroms][start:end].seq
                            m=Mutation(ref=ref,form=record[4])
                            segment=m.alternative
                        elif seq_type=='DEL':
                            pass
                        else:
                            raise ShouldNotBeHereError
                        seq.append(segment)
                if seq:
                    outputf.write('>{}\n'.format(seq_name))
                    for outputline in pyfaidx.wrap_sequence(args.length,''.join(seq)):
                        outputf.write(outputline)
                            
class Mutation:
    '''
    Mutation form are fixed in configure file. We just need to retrieve the alternative
    allele according the reference nuleotide.
    If we do not fix the mutation form first, the same mutation (occured in the common 
    ancestor) in different individuals will have different alternative alleles.
    '''
    _mutation_matrix={'N':['N','N','N'],
                      'A':['G','C','T'],
                      'G':['A','C','T'],
                      'C':['T','A','G'],
                      'T':['C','A','G']}

    def __init__(self,ref=None,form=None):
        self.ref=ref.upper()
        self.form=int(form)
        try:
            self.alternative=Mutation._mutation_matrix[self.ref][self.form]
        except KeyError as e:
            raise Exception('{} is not a nucleotide, which is found in your reference fasta file.'.format(ref)) from e
        except IndexError as e:
            raise Exception('{} is not a valid mutation form, which is found in your genome configure file.'.format(form)+
                '\nOnly 0 (transition), 1 and 2 (transversion) are acceptable.') from e

class ShouldNotBeHereError(Exception):
    pass
         
if __name__ == '__main__':
    main()
