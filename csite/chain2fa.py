#!/usr/bin/env python3

#########################################################################
# Author: Hechuan Yang
# Created Time: 2017-09-12 16:28:34
# File Name: chain2fa.py
# Description: 
#########################################################################

import os
import sys
import re
import argparse
import pyfaidx
import glob
import numpy
from csite.vcf2fa import check_output_folder

#handle the error below
#python | head == IOError: [Errno 32] Broken pipe 
from signal import signal, SIGPIPE, SIG_DFL 
signal(SIGPIPE,SIG_DFL) 

def check_chain_folder(directory=None):
    if not os.path.isdir(directory):
        raise argparse.ArgumentTypeError("'{}' doesn't exist or isn't a folder.".format(directory))
    return directory
    
def check_normal_fastas(fastas=None):
    for fa in fastas.split(','):
        if not os.path.isfile(fa):
            raise argparse.ArgumentTypeError("'{}' doesn't exist or isn't a file.".format(fa))
    return fastas
    
def main(progname=None):
    parse=argparse.ArgumentParser(
        description='Build tumor genomes from somatic variants (encoded in the chain file)',
        prog=progname if progname else sys.argv[0])
    parse.add_argument('-c','--chain',required=True,type=check_chain_folder,metavar='DIR',
        help='the folder containing the chain files of tumor genomes')
    parse.add_argument('-n','--normal',required=True,type=check_normal_fastas,metavar='FILES',
        help='two fasta files (seperated by comma) of normal genome')
    default='tumor_fa'
    parse.add_argument('-o','--output',default=default,type=check_output_folder,metavar='DIR',
        help='output directory [{}]'.format(default))
    default=50
    parse.add_argument('-w','--width',default=default,type=int,metavar='INT',
        help='the line width of output fasta files [{}]'.format(default))
    args=parse.parse_args()

    refs=[]
    for fa in args.normal.split(','):
        refs.append(pyfaidx.Fasta(fa))
    os.mkdir(args.output,mode=0o755)
    parentalre=re.compile('^parental:[01]$')
    for node_chain in glob.glob(args.chain+'/node*.chain'):
        node=os.path.basename(node_chain)
        node=node.split('.')[0]
        outputf=[]
        for parental in 0,1:
            outputf.append(open('{}/{}.parental_{}.fa'.format(args.output,node,parental),'w'))
        reference=None
        with open(node_chain) as inputf:
            seq_name=None
            parental=None
            seq=[]
            for line in inputf:
                line=line.rstrip()
                if line.startswith('>'):
                    if seq:
                        outputf[parental].write('>{}\n'.format(seq_name))
                        for outputline in pyfaidx.wrap_sequence(args.width,''.join(seq)):
                            outputf[parental].write(outputline)
                    seq_name,parental=line[1:].split()
                    if parentalre.match(parental):
                        parental=int(parental.split(':')[1])
                        try:
                            reference=refs[parental]
                        except IndexError:
                            raise FastaMissingError('There is no parental {} avalible,\n'.format(parental)+
                                'which is required in the record ({}):\n{}\n'.format(node_chain,line))
                    else:
                        raise ChainFileError('The format of this line below from the chain file '+
                            '({}) is not correct:\n{}\n'.format(node_chain,line))
                    seq=[]
                else:
                    column=line.split()
                    chroms=column[0]
                    start=int(column[1])
                    end=int(column[2])
                    seq_type=column[3]
                    segment=''
                    if seq_type=='REF':
                        segment=reference[chroms][start:end].seq
                    elif seq_type=='SNV':
                        ref=reference[chroms][start:end].seq
                        m=Mutation(ref=ref,form=column[4])
                        segment=m.alternative
                    elif seq_type=='DEL':
                        pass
                    else:
                        raise ChainFileError('Can not recognize the sequence type ({}) '.format(seq_type)+
                            'of the record below from the chain file ({}):\n{}\n'.format(node_chain,line))
                    seq.append(segment)
            if seq:
                outputf[parental].write('>{}\n'.format(seq_name))
                for outputline in pyfaidx.wrap_sequence(args.width,''.join(seq)):
                    outputf[parental].write(outputline)
        for parental in 0,1:
            outputf[parental].close()
                        
class Mutation:
    '''
    Mutation form are fixed in chain files. We just need to retrieve the alternative
    allele according the reference nuleotide.
    If we do not fix the mutation form first, the same mutation (occured in the common 
    ancestor) in different individuals may have different alternative alleles.
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
            raise FastaFileError("'{}' is not a nucleotide, but it's found in your normal fasta file.".format(ref)) from e
        except IndexError as e:
            raise ChainFileError("'{}' is not a valid mutation form, but it's found in your chain file.".format(form)+
                '\nOnly 0 (transition), 1 and 2 (transversion) are acceptable.') from e

class ChainFileError(Exception):
    pass
         
class FastaMissingError(Exception):
    pass
         
class FastaFileError(Exception):
    pass
         
if __name__ == '__main__':
    main()
