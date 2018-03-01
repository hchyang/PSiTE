#!/usr/bin/env python

import argparse


def main():
    parser = argparse.ArgumentParser(description='probe2fa: a program to convert probe sequences to FASTA format', prog='probe2fa', formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i', metavar='FILE', dest='probe_file', required=True,
                        help='The file containing sequences of probes for targeted capture sequencing.')
    parser.add_argument('-o', metavar='FILE', dest='fasta_file', required=True,
                        help='The probe file in FASTA format.')
    args = parser.parse_args()

    with open(args.probe_file,'r') as fin, open(args.fasta_file, 'w') as fout:
        fin.readline()
        fin.readline()
        for line in fin:
            values = line.strip().split("\t")
            if len(values) < 3:
                line = fin.readline()
                continue
            target = values[0]
            seqid = values[1]
            seq = values[2]
            fout.write(">" + seqid + "-" + target + "\n")
            fout.write(seq + "\n")

if __name__=="__main__":
    main()
