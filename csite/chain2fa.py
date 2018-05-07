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
import multiprocessing
from csite.vcf2fa import check_output_folder

#handle the error below
#python | head == IOError: [Errno 32] Broken pipe 
from signal import signal, SIGPIPE, SIG_DFL 
signal(SIGPIPE,SIG_DFL) 

def check_folder(directory=None):
    if not os.path.isdir(directory):
        raise argparse.ArgumentTypeError("'{}' doesn't exist or isn't a folder.".format(directory))
    return directory
    
def check_normal_fastas(fastas=None):
    for fa in fastas.split(','):
        if not os.path.isfile(fa):
            raise argparse.ArgumentTypeError("'{}' doesn't exist or isn't a file.".format(fa))
    return fastas
    
def main(progname=None):
    parser=argparse.ArgumentParser(
        description='Build tumor genomes from somatic variants (encoded in the chain file)',
        prog=progname if progname else sys.argv[0])
    parser.add_argument('-c','--chain',required=True,type=check_folder,metavar='DIR',
        help='the folder containing the chain files of tumor genomes')
    parser.add_argument('-n','--normal',required=True,type=check_normal_fastas,metavar='FILES',
        help='two fasta files (separated by comma) of normal genome')
    default='tumor_fa'
    parser.add_argument('-o','--output',default=default,type=check_output_folder,metavar='DIR',
        help='output directory [{}]'.format(default))
    default=50
    parser.add_argument('-w','--width',default=default,type=int,metavar='INT',
        help='the line width of output fasta files [{}]'.format(default))
    default=1
    parser.add_argument('--cores',type=int,default=default,metavar='INT',
        help='number of cores used to run the program [{}]'.format(default))

    args=parser.parse_args()

    os.mkdir(args.output,mode=0o755)
    normal_fa=args.normal.split(',')
    for fa in normal_fa:
        pyfaidx.Faidx(fa)
    pool=multiprocessing.Pool(processes=args.cores)
    results=[]
    for node_chain in glob.glob(os.path.join(args.chain,'node*.chain')):
        results.append(pool.apply_async(build_fasta,args=(args.output,node_chain,normal_fa,args.width)))
    pool.close()
    pool.join()
#handle exceptions if any
    for result in results:
        result.get()

def build_fasta(output=None,chain=None,normal_fa=None,width=None):
    refs=[]
    for fa in normal_fa:
        refs.append(pyfaidx.Fasta(fa))
    parentalre=re.compile('^parental:[01]$')
    node=os.path.basename(chain)
    node=node.split('.')[0]
    outputf=[]
    for parental in 0,1:
        outputf.append(open(os.path.join(output,'{}.parental_{}.fa'.format(node,parental)),'w'))
    reference=None
    with open(chain) as inputf:
        seq_name=None
        parental=None
        seq=[]
        for line in inputf:
            line=line.rstrip()
            if line.startswith('>'):
                if seq:
                    outputf[parental].write('>{}\n'.format(seq_name))
                    for outputline in pyfaidx.wrap_sequence(width,''.join(seq)):
                        outputf[parental].write(outputline)
                seq_name,parental=line[1:].split()
                if parentalre.match(parental):
                    parental=int(parental.split(':')[1])
                    try:
                        reference=refs[parental]
                    except IndexError as e:
                        raise FastaMissingError('There is no parental {} avalible,\n'.format(parental)+
                            'which is required in the record ({}):\n{}\n'.format(chain,line)) from e
                else:
                    raise ChainFileError('The format of this line below from the chain file '+
                        '({}) is not correct:\n{}'.format(chain,line))
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
                    try:
                        form=column[4]
                    except IndexError:
                        raise ChainFileError('Can not found mutation form in the record below ({}):\n{}'.format(chain,line))
                    ref=reference[chroms][start:end].seq
                    m=Mutation(ref=ref,form=form)
                    segment=m.alternative
                    if segment==KeyError:
                        raise FastaFileError("'{}' is not a nucleotide, ".format(ref)+
                            "but it's found in your normal fasta file ({}[{}:{}]).".format(normal_fa[parental],chroms,end))
                    elif segment==IndexError:
                        raise ChainFileError("'{}' is not a valid mutation form of SNV,\n".format(form)+
                            "but it's found in your chain file ({}):\n{}".format(chain,line))
                elif seq_type=='DEL':
                    pass
                else:
                    raise ChainFileError('Can not recognize the sequence type ({}) '.format(seq_type)+
                        'of the record below from the chain file ({}):\n{}\n'.format(chain,line))
                seq.append(segment)
        if seq:
            outputf[parental].write('>{}\n'.format(seq_name))
            for outputline in pyfaidx.wrap_sequence(width,''.join(seq)):
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
        except KeyError:
            self.alternative=KeyError
        except IndexError:
            self.alternative=IndexError

class ChainFileError(Exception):
    pass
         
class FastaMissingError(Exception):
    pass
         
class FastaFileError(Exception):
    pass
         
if __name__ == '__main__':
    main()
