#!/usr/bin/env python3

#########################################################################
# Author: Hechuan Yang
# Created Time: 2017-09-27 10:50:14
# File Name: psite.py
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

__version__='0.9.1'

def usage():
    print("")
    print("Program: psite.py (a Phylogeny guided Simulator for Tumor Evolution)")
    print("Version: {}".format(__version__))
    print("")
    print("Usage:   psite.py <command> [options]")
    print("")
    print("Command: vcf2fa     build normal genome from input germline vcf file")
    print("         phylovar   simulate somatic variants along a phylogeny")
    print("         chain2fa   build tumor genomes from somatic variants (encoded in chain files)")
    print("         fa2wgs     simulate WGS reads from normal and tumor genomes (in fasta format)")
    print("         fa2wes     simulate WES reads from normal and tumor genomes (in fasta format)")
    print("         allinone   a wrapper for NGS reads simulation by combining all individual steps")
    print("")

def main():
    if len(sys.argv)==1 or sys.argv[1]=='-h' or sys.argv[1]=='--help':
        usage()
    else:
        command=sys.argv[1]
        progname='psite.py '+command
        del sys.argv[1]
        if command=='vcf2fa':
            import psite.vcf2fa
            psite.vcf2fa.main(progname=progname)
        elif command=='phylovar':
            import psite.phylovar
            psite.phylovar.main(progname=progname)
        elif command=='chain2fa':
            import psite.chain2fa
            psite.chain2fa.main(progname=progname)
        elif command=='fa2wgs':
            import psite.fa2wgs
            psite.fa2wgs.main(progname=progname)
        elif command=='fa2wes':
            import psite.fa2wes
            psite.fa2wes.main(progname=progname)
        elif command=='allinone':
            import psite.allinone
            psite.allinone.main(progname=progname)
        else:
            print("[psite.py] Unrecognized command: '{}'".format(command))
            sys.exit()


if __name__ == '__main__':
    main()
