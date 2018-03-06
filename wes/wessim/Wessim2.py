#!/usr/bin/env python

import sys
import bisect
import pysam
import gzip
import numpy
from time import time, localtime, strftime
import argparse
from multiprocessing import Process
import os
import math
import random

inds={'A':0,'T':1,'G':2,'C':3,'N':4,'a':0,'t':1,'g':2,'c':3,'n':4}

MAX_INT = 2**16
def random_int():
    return random.randint(0, MAX_INT)

def subprogram(command, name):
    os.system(command)
    print("exiting subprocess " + str(name))

def main(argv):
    t0 = time()
    arguline = " ".join(argv)
    parser = argparse.ArgumentParser(description='Wessim2: Whole Exome Sequencing SIMulator 2 (Probe-based version)', prog='Wessim2', formatter_class=argparse.RawTextHelpFormatter)
    group1 = parser.add_argument_group('Input options')
    group1.add_argument('-R', metavar = 'FILE', dest='reference', required=True, help='faidx-indexed (R)eference genome FASTA file or meta description file (.meta)')
    group1.add_argument('-P', metavar = 'FILE', dest='probe', required=True, help='(P)robe sequence FASTA file')
    group1.add_argument('-B', metavar = 'FILE', dest='probeblat', required=True, help='(B)lat matched probe regions .PSL file')
    default = 0
    group1.add_argument('-s', metavar='INT', type=int, dest='random_seed', help='The seed for random number generator [{}]'.format(default))

    group2 = parser.add_argument_group('Parameters for exome capture')
    group2.add_argument('-f', metavar = 'INT', type=int, dest='fragsize', required=False, help='mean (f)ragment size. this corresponds to insert size when sequencing in paired-end mode. [200]', default=200)
    group2.add_argument('-d', metavar = 'INT', type=int, dest='fragsd', required=False, help='standard (d)eviation of fragment size [50]', default=50)
    group2.add_argument('-m', metavar = 'INT', type=int, dest='fragmin', required=False, help='(m)inimum fragment length [read_length + 20]')
    group2.add_argument('-y', metavar = 'PERCENT',type=int, dest='bind', required=False, help='minimum required fraction of probe match to be h(y)bridized [50]', default=50)
    group2.add_argument('-w', metavar = 'INT', type=int, dest='weight', required=False, help='penalty (w)eight for indel in the hybridization [2]', default=2)

    group3 = parser.add_argument_group('Parameters for sequencing')
    group3.add_argument('-p', action='store_true', help='generate paired-end reads [single]')
    group3.add_argument('-n', metavar = 'INT', type=int, dest='readnumber', required=True, help='total (n)umber of reads')
    group3.add_argument('-l', metavar = 'INT', type=int, dest='readlength', required=True, help='read (l)ength (bp)')
    group3.add_argument('-M', metavar = 'FILE', dest='model', required=True, help='GemSim (M)odel file (.gzip)')
    group3.add_argument('-t', metavar = 'INT', type=int, dest='threadnumber', required=False, help='number of (t)hreaded subprocesses [1]', default=1)

    group4 = parser.add_argument_group('Output options')
    group4.add_argument('-o', metavar = 'FILE', dest='outfile', help='(o)utput file header. ".fastq.gz" or ".fastq" will be attached automatically. Output will be splitted into two files in paired-end mode', required=True)
    group4.add_argument('-z', action='store_true', help='compress output with g(z)ip [false]')
    group4.add_argument('-q', metavar = 'INT', type=int, dest='qualbase', required=False, help='(q)uality score offset [33]', default=33)
    group4.add_argument('-v', action='store_true', help='(v)erbose; print out intermediate messages.')

    args = parser.parse_args()
    reffile = args.reference
    probefile = args.probe
    alignfile = args.probeblat

    if args.random_seed == None:
        seed = random_int()
    else:
        seed = args.random_seed
    print('Random seed for Wessim2: {}'.format(seed))
    random.seed(seed)

    isize = args.fragsize
    isd = args.fragsd
    imin = args.fragmin
    bind = args.bind

    paired = args.p
    readlength = args.readlength
    readnumber = args.readnumber
    threadnumber = args.threadnumber
    if imin==None:
        if paired:
            imin = readlength + 20
        else:
            imin = readlength + 20
    if isize < imin:
        print("too small mean fragment size (" + str(isize) + ") compared to minimum length (" + str(imin) + "). Increase it and try again.")
        sys.exit(0)
    model = args.model

    outfile = args.outfile
    compress = args.z
    qualbase = args.qualbase
    verbose = args.v

    print("-------------------------------------------")
    print("Reference:", reffile)
    print("Probeset:", probefile)
    print("Probematch:", alignfile)
    print("Fragment:",isize, "+-", isd, ">", imin)
    print("Paired-end mode?", paired)
    print("Sequencing model:", model)
    print("Read length:", readlength, "Read number:", readnumber)
    print("Output File:", outfile)
    print("Gzip compress?", compress)
    print("Quality base:", qualbase)
    print("Thread number:", threadnumber)
    print("Job started at:", strftime("%Y-%m-%d %H:%M:%S", localtime()))
    print("-------------------------------------------")


    processes = []
    for t in range(0, threadnumber):
        readstart = int(float(readnumber) / float(threadnumber) * t) + 1
        readend = int(float(readnumber) / float(threadnumber) * (t+1))
        script_dir = os.path.dirname(os.path.realpath(__file__))
        script = os.path.join(script_dir , "__sub_wessim2.py")
        rgid = t + 1
        # command = "python " + script + " " + arguline + " -1 " + str(readstart) + " -2 " + str(readend) + " -i " + str(t+1)
        seed = random_int()
        command = "python {} {} -1 {} -2 {} -i {} -s {}".format(script, arguline, readstart, readend, rgid, seed)
        print(command)
        p = Process(target=subprogram, args=(command, t+1))
        p.start()
        processes.append(p)
    for p in processes:
        p.join()
    t1 = time()
    print("Done generating " + str(readnumber) + " reads in %f secs" % (t1 - t0))
    print("Merging subresults...")
    wread = None
    wread2 = None
    if paired and compress:
        f1 = outfile + "_1.fastq.gz"
        wread = gzip.open(f1, 'wb')
        f2 = outfile + "_2.fastq.gz"
        wread2 = gzip.open(f2, 'wb')
    elif paired and not compress:
        wread = open(outfile + "_1.fastq", 'w')
        wread2 = open(outfile + "_2.fastq", 'w')
    elif not paired and compress:
        wread = gzip.open(outfile + ".fastq.gz", 'wb')
    else:
        wread = open(outfile + ".fastq", 'w')
    if not paired:
        for t in range(0, threadnumber):
            suboutfile = outfile + "-" + str(t+1)
            fread = None
            if compress:
                suboutfile += ".fastq.gz"
                fread = gzip.open(suboutfile, 'rb')
            else:
                suboutfile += ".fastq"
                fread = open(suboutfile, 'r')
            line = fread.readline()
            while line:
                wread.write(line)
                line = fread.readline()
            fread.close()
            os.remove(suboutfile)
        wread.close()
    else:
        for t in range(0, threadnumber):
            suboutfile1 = outfile + "-" + str(t+1) + "_1"
            suboutfile2 = outfile + "-" + str(t+1) + "_2"
            fread1 = None
            fread2 = None
            if compress:
                suboutfile1 += ".fastq.gz"
                suboutfile2 += ".fastq.gz"
                fread1 = gzip.open(suboutfile1, "rb")
                fread2 = gzip.open(suboutfile2, "rb")
            else:
                suboutfile1 += ".fastq"
                suboutfile2 += ".fastq"
                fread1 = open(suboutfile1, "r")
                fread2 = open(suboutfile2, "r")
            line1 = fread1.readline()
            line2 = fread2.readline()
            while line1 and line2:
                wread.write(line1)
                wread2.write(line2)
                line1 = fread1.readline()
                line2 = fread2.readline()
            fread1.close()
            fread2.close()
            os.remove(suboutfile1)
            os.remove(suboutfile2)
        wread.close()
        wread2.close()
    sys.exit(0)


if __name__=="__main__":
    main(sys.argv[1:])
