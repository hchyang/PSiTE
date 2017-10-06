#!/usr/bin/env python3

#########################################################################
# Author: Hechuan Yang
# Created Time: 2017-09-27 10:50:14
# File Name: csite.py
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

__version__='0.9.0'

def usage():
    print("")
    print("Program: csite.py (a Coalescent Simulator for Tumor Evolution)")
    print("Version: {}".format(__version__))
    print("")
    print("Usage:   csite.py <command> [options]")
    print("")
    print("Command: germlinev   integrate germline snp")
    print("         somaticv    simulate somatic vars")
    print("         buildref    build reference")
    print("")

def main():
    if len(sys.argv)==1:
        usage()
    else:
        progname=' '.join(sys.argv[0:2])
        if sys.argv[1]=='germlinev':
            import csite.integrate_germline_snp 
            csite.integrate_germline_snp.main(progname=progname)
        elif sys.argv[1]=='somaticv':
            import csite.somatic_sim
            csite.somatic_sim.main(progname=progname)
        elif sys.argv[1]=='buildref':
            import csite.build_reference
            csite.build_reference.main(progname=progname)
        else:
            print('Do not have this command in csite: {}'.format(sys.argv[1]))
            exit()
        
if __name__ == '__main__':
    main()

