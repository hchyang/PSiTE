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

__version__='0.9.0'

def usage():
    print("")
    print("Program: csite.py (a Coalescent Simulator for Tumor Evolution)")
    print("Version: {}".format(__version__))
    print("")
    print("Usage:   csite.py <command> [options]")
    print("")
    print("Command: vcf2fa      integrate germline snp")
    print("         phylovar    simulate somatic vars")
    print("         draft2fa    build reference")
    print("         quaternity  a wrapper for the whole pipeline")
    print("")

def main():
    if len(sys.argv)==1 or sys.argv[1]=='-h':
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
        elif command=='draft2fa':
            import csite.draft2fa
            csite.draft2fa.main(progname=progname)
        elif command=='quaternity':
            import csite.quaternity
            csite.quaternity.main(progname=progname)
        else:
            print('Do not have this command in csite: {}'.format(command))
            exit()
        
if __name__ == '__main__':
    main()

