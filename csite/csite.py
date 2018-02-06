#!/usr/bin/env python3

#########################################################################
# Author: Hechuan Yang
# Created Time: 2017-09-27 10:50:14
# File Name: csite.py
# Description: 
#########################################################################

import os
import sys
import logging
import pyfaidx

#handle the error below
#python | head == IOError: [Errno 32] Broken pipe 
from signal import signal, SIGPIPE, SIG_DFL 
signal(SIGPIPE,SIG_DFL) 

#If you simulate large genome for a tumor sample with many cells,
#the disk space occupied by tumor genome (fasta) will be huge.
#Because ART does not support compressed fasta, we can not output
#compressed genome fasta to save disk space.

__version__='0.9.0'

def usage():
    print("")
    print("Program: csite.py (a Coalescent Simulator for Tumor Evolution)")
    print("Version: {}".format(__version__))
    print("")
    print("Usage:   csite.py <command> [options]")
    print("")
    print("Command: vcf2fa     build normal genome from input germline vcf file")
    print("         phylovar   simulate somatic variants along a phylogeny")
    print("         chain2fa   build tumor genomes from somatic variants (encoded in the chain file)")
    print("         fa2wgs     simulate WGS reads from normal and tumor genomes (in fasta format)")
    print("         allinone   a wrapper for NGS reads simulation combining all previous steps")
    print("")

def main():
    if len(sys.argv)==1 or sys.argv[1]=='-h' or sys.argv[1]=='--help':
        usage()
    else:
        command=sys.argv[1]
        progname='csite.py '+command
        del sys.argv[1]
        if command=='vcf2fa':
            import csite.vcf2fa 
            csite.vcf2fa.main(progname=progname)
        elif command=='phylovar':
            import csite.phylovar
            csite.phylovar.main(progname=progname)
        elif command=='chain2fa':
            import csite.chain2fa
            csite.chain2fa.main(progname=progname)
        elif command=='fa2wgs':
            import csite.fa2wgs
            csite.fa2wgs.main(progname=progname)
        elif command=='allinone':
            import csite.allinone
            csite.allinone.main(progname=progname)
        else:
            print("[csite.py] Unrecognized command: '{}'".format(command))
            sys.exit()
        
if __name__ == '__main__':
    main()

