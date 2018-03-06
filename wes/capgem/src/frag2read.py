#!/usr/bin/env python

# Generate short reads from a fragment selected from a genome
import sys
import random
import bisect
import gzip
try:
   import cPickle as pickle
except:
   import pickle
import numpy
from time import time
import argparse
import math

inds = {'A': 0, 'T': 1, 'G': 2, 'C': 3, 'N': 4,
        'a': 0, 't': 1, 'g': 2, 'c': 3, 'n': 4}

MAX_INT = 2**16
def random_int():
    return random.randint(0, MAX_INT)


def main():
    t0 = time()
    parser = argparse.ArgumentParser(description='frag2read: a program to simulate short reads from fragments', prog='frag2read', formatter_class=argparse.RawTextHelpFormatter)

    group1 = parser.add_argument_group('Input options')
    group1.add_argument('-f', metavar='FILE', dest='fragment_file', required=True,
                        help='The file containing sequences of selected (f)ragment.')
    group1.add_argument('-S', metavar='INT', type=int,
                        dest='fragsize', default=200, help='mean (f)ragment size. this corresponds to insert size when sequencing in paired-end mode. 0 for empirical.')
    # group1.add_argument('-t', metavar='FILE', dest='haplotype_file',
    #                     help='The haplo(t)ype_file file specifying location and frequency of snps')
    default = 0
    group1.add_argument('-s', metavar='INT', type=int, dest='random_seed', help='The seed for random number generator [{}]'.format(default))

    group2 = parser.add_argument_group('Parameters for sequencing')
    group2.add_argument('-p', action='store_true',
                        help='generate (p)aired-end reads [single].')
    group2.add_argument('-l', metavar='INT', type=int,
                        dest='readlength', required=True, help='read (l)ength (bp). 0 for empirical distribution.')
    group2.add_argument('-M', metavar='FILE', dest='model',
                        required=True, help='GemSim (M)odel file (.gzip).')

    group3 = parser.add_argument_group('Output options')
    group3.add_argument('-o', metavar='FILE', dest='outfile',
                        help='(o)utput file header. ".fastq.gz" or ".fastq" will be attached automatically. Output will be splitted into two files in paired-end mode.', required=True)
    group3.add_argument('-z', action='store_true',
                        help='compress output with g(z)ip [false].')
    group3.add_argument('-q', metavar='INT', type=int, dest='qualbase',
                        required=False, help='(q)uality score offset [33].', default=33)

    args = parser.parse_args()

    fragment_file = args.fragment_file
    # haplotype_file = args.haplotype_file
    fragsize = args.fragsize

    if args.random_seed == None:
        seed = random_int()
    else:
        seed = args.random_seed
    print('Random seed: {}'.format(seed))
    random.seed(seed)

    paired = args.p
    readlength = args.readlength
    model = args.model

    outfile = args.outfile
    compress = args.z
    qualbase = args.qualbase

    # Parse known error models
    if paired:
        mx1, mx2, insD1, insD2, delD1, delD2, intervals, gQualL, bQualL, iQualL, mates, rds, rdLenD = parseModel(
            model, paired, readlength)
        m0 = float(mates[0])
        m1 = float(mates[1])
        rd0 = float(rds[0])
        rd1 = float(rds[1])
        unAlign0 = (m0 * rd1 - m1 * m0) / (rd0 * rd1 - m1 * m0)
        unAlign1 = 1.0 - (unAlign0 / (m0 / rd0))
        keys = intervals.keys()
        # keys.sort()
        keys = sorted(keys)
        if fragsize == 0:
            inters = []
            for k in keys:
                inters.append((k, intervals[k]))
            interval = bisect_choiceTUP(inters)
        # inserts1and2
        insDict1 = mkInserts(mx1, insD1)
        insDict2 = mkInserts(mx2, insD2)
        # deletions1and2
        delDict1 = mkDels(mx1, delD1)
        delDict2 = mkDels(mx2, delD2)
    else:
        mx1, insD1, delD1, gQualL, bQualL, iQualL, readCount, rdLenD = parseModel(
            model, paired, readlength)
        insDict = mkInserts(mx1, insD1)
        # deletions
        delDict = mkDels(mx1, delD1)

    gens = genRef('')

    # choose good quality bases
    gQList = []
    for i in (gQualL):
        gL = []
        keys = i.keys()
        # keys.sort()
        keys = sorted(keys)
        for k in keys:
            gL.append((chr(k + qualbase), i[k]))
        gQList.append(bisect_choiceTUP(gL))

    # choose bad quality bases
    bQList = []
    for i in (bQualL):
        bL = []
        keys = i.keys()
        # keys.sort()
        keys = sorted(keys)
        for k in keys:
            bL.append((chr(k + qualbase), i[k]))
        bQList.append(bisect_choiceTUP(bL))

    # choose qualities for inserts
    iQList = []
    for i in (iQualL):
        iL = []
        keys = i.keys()
        # keys.sort()
        keys = sorted(keys)
        for k in keys:
            iL.append((chr(k + qualbase), i[k]))
        iQList.append(bisect_choiceTUP(iL))

    # Generate!
    wread = None
    wread2 = None
    if paired and compress:
        wread = gzip.open(outfile + "_1.fastq.gz", 'wb')
        wread2 = gzip.open(outfile + "_2.fastq.gz", 'wb')
    elif paired and not compress:
        wread = open(outfile + "_1.fastq", 'w')
        wread2 = open(outfile + "_2.fastq", 'w')
    elif not paired and compress:
        wread = gzip.open(outfile + ".fastq.gz", 'wb')
    else:
        wread = open(outfile + ".fastq", 'w')
    dirtag = ('', '+', '-')
    count = 0
    seqgenome = "g1"
    with open(fragment_file, 'r') as fin:
        i = 0
        for seq in fin:
            seq = seq.strip()
            i += 1
            fragment_chrom, fragment_start, ref = seq.split('\t')
            refLen = len(ref)
            if refLen < readlength:
                continue
            if not paired:
                readLen = readlength
                read1, pos, dir, quals1 = readGen1(ref, refLen, readLen, gens(
                ), readLen, mx1, insDict, delDict, gQList, bQList, iQList, qualbase)
                if read1 is None:
                    continue
                head1 = '@' + 'r' + str(i) + '_from_' + seqgenome + ";" + fragment_chrom + \
                    "_" + str(int(fragment_start) + pos + 1) + "_" + dirtag[dir]
            else:
                val = random.random()
                ln1 = readlength
                ln2 = readlength
                inter = fragsize
                read1, pos1, dir1, quals1, read2, pos2, dir2, quals2 = readGenp(
                    ref, refLen, ln1, ln2, gens(), mx1, insDict1, delDict1, gQList, bQList, iQList, qualbase)
                if read1 is None or read2 is None:
                    continue
                p1 = fragment_chrom + "_" + \
                    str(int(fragment_start) + pos1 + 1) + "_" + dirtag[dir1]
                p2 = fragment_chrom + "_" + \
                    str(int(fragment_start) + pos2 + 1) + "_" + dirtag[dir2]
                if val > unAlign0 + unAlign1:
                    pass
                elif val > unAlign1:
                    read2 = 'N' * ln2
                    quals2 = chr(0 + qualbase) * ln2
                    p2 = '*'
                else:
                    read1 = 'N' * ln1
                    quals1 = chr(0 + qualbase) * ln1
                    p1 = '*'
                head1 = '@' + 'r' + str(i) + '_from_' + \
                    seqgenome + ";" + p1 + ":" + p2 + "/1"
                head2 = '@' + 'r' + str(i) + '_from_' + \
                    seqgenome + ";" + p1 + ":" + p2 + "/2"
            if compress:
                wread.write(head1.encode() + b'\n')
                wread.write(read1.upper().encode() + b'\n')
                wread.write(b'+\n')
                wread.write(quals1.encode() + b'\n')
                if paired:
                    wread2.write(head2.encode() + b'\n')
                    wread2.write(read2.upper().encode() + b'\n')
                    wread2.write(b'+\n')
                    wread2.write(quals2.encode() + b'\n')
            else:
                wread.write(head1 + '\n')
                wread.write(read1.upper() + '\n')
                wread.write('+\n')
                wread.write(quals1 + '\n')
                if paired:
                    wread2.write(head2 + '\n')
                    wread2.write(read2.upper() + '\n')
                    wread2.write('+\n')
                    wread2.write(quals2 + '\n')
            count += 1
            if count % 1000000 == 0 and count != 0:
                t = time() - t0
                print("{} reads have been generated... in {} secs".format(count, t))

    wread.close()
    if paired:
        wread2.close()


def comp(sequence):
    """ complements a sequence, preserving case. Function imported from GemSim"""
    d = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'a': 't',
         't': 'a', 'c': 'g', 'g': 'c', 'N': 'N', 'n': 'n'}
    cSeq = ''
    for s in sequence:
        if s in d.keys():
            cSeq += d[s]
        else:
            cSeq += 'N'
    return cSeq


def parseModel(gzipFile, paired, readlen):
    """prepares error models for input to mkErrors."""
    file = gzip.open(gzipFile, 'rb')
    if paired:
        modReadLen = pickle.load(file)
        if readlen != 'd' and readlen > modReadLen:
            print("Inappropriate read length chosen for model. Maximum for this model: {}".format(modReadLen))
            file.close()
            sys.exit()
        mx1 = pickle.load(file)
        mx2 = pickle.load(file)
        insD1 = pickle.load(file)
        insD2 = pickle.load(file)
        delD1 = pickle.load(file)
        delD2 = pickle.load(file)
        intD = pickle.load(file)
        gQualL = pickle.load(file)
        bQualL = pickle.load(file)
        iQualL = pickle.load(file)
        mates = pickle.load(file)
        rds = pickle.load(file)
        rdLenD = pickle.load(file)
        file.close()
        return mx1, mx2, insD1, insD2, delD1, delD2, intD, gQualL, bQualL, iQualL, mates, rds, rdLenD
    else:
        modReadLen = pickle.load(file)
        if readlen != 'd' and readlen > modReadLen:
            print("Inappropriate read length chosen for model. Maximum for this model: {}".format(modReadLen))
            file.close()
            sys.exit()
        mx = pickle.load(file)
        insD = pickle.load(file)
        delD = pickle.load(file)
        gQualL = pickle.load(file)
        bQualL = pickle.load(file)
        iQualL = pickle.load(file)
        readCount = pickle.load(file)
        rdLenD = pickle.load(file)
        file.close()
        return mx, insD, delD, gQualL, bQualL, iQualL, readCount, rdLenD


def mkInserts(mx, insD):
    """Returns a dictionary consisting of compiled functions to make inserts."""
    insertDict = {}
    posKeys = insD.keys()
    # posKeys.sort()
    posKeys = sorted(posKeys)
    for p in posKeys:
        indicies = p.split('.')
        tot = mx[int(indicies[0])][int(indicies[1])][int(indicies[2])][int(
            indicies[3])][int(indicies[4])][int(indicies[5])][5]
        insertKeys = insD[p].keys()
        # insertKeys.sort()
        insertKeys = sorted(insertKeys)
        insertList = []
        iSum = 0
        for i in insertKeys:
            insertList.append((i, insD[p][i]))
            iSum += 0
        insertList.append(('', tot - iSum))
        insert = bisect_choiceTUP(insertList)
        insertDict[p] = insert
    return insertDict


def mkDels(mx, delD):
    """Returns a dictionary consisting of compiled functions to make deletiosn."""
    deletionDict = {}
    posKeys = delD.keys()
    # posKeys.sort()
    posKeys = sorted(posKeys)
    for p in posKeys:
        indicies = p.split('.')
        tot = mx[int(indicies[0])][int(indicies[1])][int(indicies[2])][int(
            indicies[3])][int(indicies[4])][int(indicies[5])][5]
        items = delD[p]
        items.reverse()
        items.append(tot - sum(items))
        items.reverse()
        delete = bisect_choice(items)
        deletionDict[p] = delete
    return deletionDict


def bisect_choice(items):
    """Returns a function that makes a weighted random choice from items."""
    added_weights = []
    last_sum = 0
    for weight in items:
        last_sum += weight
        added_weights.append(last_sum)

    def choice(rnd=random.random, bis=bisect.bisect):
        return bis(added_weights, rnd() * last_sum)
    return choice


def bisect_choiceTUP(items):
    """Returns a function that makes a weighted random choice from a list of tuples."""
    added_weights = []
    last_sum = 0.0
    for item, weight in items:
        weight = float(weight)
        last_sum += weight
        added_weights.append(last_sum)

    def choice(rnd=random.random, bis=bisect.bisect):
        return items[bis(added_weights, rnd() * last_sum)][0]
    return choice


def ln(length):
    """Returns static length as a funtion."""
    def val():
        return length
    return val


def readGen1(ref, refLen, readLen, genos, inter, mx1, insD1, delD1, gQ, bQ, iQ, qual):
    """Generates a random read of desired length from a reference."""
    assert refLen >= readLen
    ind = 0
    dir = random.randint(1, 2)
    end = ind + readLen
    read = ref[ind:end]

    if dir == 2:
        cRef = comp(ref)[::-1]
        read = cRef[refLen - end:refLen - ind]
    # if genos != '':
    #     read = mutate(read, ind, genos, refLen, 1, readPlus, hd)
    read, quals = mkErrors(read, readLen, mx1, insD1, delD1, gQ, bQ, iQ, qual)

    return read, ind, dir, quals


def readGenp(ref, refLen, readLen1, readLen2, genos, mx1, insD1, delD1, gQ, bQ, iQ, qual):
    """Generates a pair of reads from given DNA fragment."""
    assert refLen >= readLen1 and refLen >= readLen2
    cRef = comp(ref)[::-1]
    ind1 = 0
    ind2 = refLen - readLen2
    end1 = ind1 + readLen1
    dir1 = 1
    dir2 = 2
    read1 = ref[ind1:end1]
    read2 = cRef[ind1:end1]
    read1, quals1 = mkErrors(read1, readLen1, mx1,
                             insD1, delD1, gQ, bQ, iQ, qual)
    read2, quals2 = mkErrors(read2, readLen2, mx1,
                             insD1, delD1, gQ, bQ, iQ, qual)
    pairorder = random.randint(1, 2)
    if pairorder == 1:
        # r1=[read1, ind1, dir1, quals1, read2, ind2, dir2, quals2]
        # print('r1: {}\n'.format(str(r1)))
        return read1, ind1, dir1, quals1, read2, ind2, dir2, quals2
    else:
        # r2=[read2, ind2, dir2, quals2, read1, ind1, dir1, quals1]
        # print('r2: {}\n'.format(str(r2)))
        return read2, ind2, dir2, quals2, read1, ind1, dir1, quals1


def mutate(read, ind, gens, refLen, dir, readLn, hd):
    """Adds predetermined mutations to reads based on known haplotype."""
    d = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'a': 't',
         't': 'a', 'c': 'g', 'g': 'c', 'N': 'N', 'n': 'n'}
    if gens == {}:
        return read
    else:
        chroms = gens.keys()
        if hd not in chroms:
            return read
        else:
            posi = gens[hd].keys()
            if dir == 1:
                for p in posi:
                    if p > ind and p <= (ind + readLn):
                        read1 = read[:p - (ind + 1)] + gens[hd][p]
                        read1 = read1 + read[p - ind:]
                        read = read1
                    elif p <= ind + readLn - refLen:
                        read1 = read[:refLen - ind + p - 1] + gens[hd][p]
                        read1 += read[refLen - ind + p:]
                        read = read1
                return read
            elif dir == 2:
                for p in posi:
                    if p > ind and p <= (ind + readLn):
                        read1 = read[:p - (ind + 1)] + d[gens[hd][p]]
                        read1 = read1 + read[p - ind:]
                        read = read1
                    elif p <= ind + readLn - refLen:
                        read1 = read[:refLen - ind + p - 1] + d[gens[hd][p]]
                        read1 += read[refLen - ind + p:]
                        read = read1
                return read


def genRef(ref):
    """Returns input as function"""
    def r():
        return ref
    return r


def mkErrors(read, readLen, mx, insD, delD, gQ, bQ, iQ, qual):
    """Adds random errors to read."""
    pos = 0
    quals = ''
    qualslist = []
    index = '0.4.4.4.4.' + str(inds[read[0]])
    if index in insD:
        insert = insD[index]()
        read = 'NNNN' + insert + read
        for i in insert:
            #			quals+=iQ[0]()
            qualslist.append(iQ[0]())
            pos += 1
    else:
        read = 'NNNN' + read
    prev = read[pos:pos + 4]
    after = read[pos + 4]
    d0 = pos
    d1 = inds[prev[3]]
    d2 = inds[prev[2]]
    d3 = inds[prev[1]]
    d4 = inds[prev[0]]
    d5 = inds[after]
    pos += 1
    while pos <= readLen and pos < len(read) - 4:
        d0 = pos
        d4 = d3
        d3 = d2
        d2 = d1
        d1 = d5
        d5 = inds[read[pos + 4]]
        index = '.'.join([str(d0), str(d1), str(d2),
                          str(d3), str(d4), str(d5)])
        Mprobs = mx[d0][d1][d2][d3][d4][d5]
        tot = float(Mprobs[5])
        if not tot == 0:
            Mprobs = Mprobs / tot
        val = random.random()
        a = Mprobs[0]
        t = Mprobs[1] + a
        g = Mprobs[2] + t
        c = Mprobs[3] + g
        n = Mprobs[4] + c
        success = False
        if val > n or tot == 0:
            gPos = pos - 1
            while gPos >= 0:
                try:
                    qualslist.append(gQ[gPos]())
                    success = True
                    break
                except:
                    gPos -= 1
            if success == False:
                qualslist.append(chr(30 + qual))
        elif val > c:
            read = read[:pos + 3] + 'N' + read[pos + 4:]
            bPos = pos - 1
            while bPos >= 0:
                try:
                    qualslist.append(bQ[bPos]())
                    success = True
                    break
                except:
                    bPos - 1
                if success == False:
                    qualslist.append(chr(2 + qual))
        elif val > g:
            read = read[:pos + 3] + 'C' + read[pos + 4:]
            bPos = pos - 1
            while bPos >= 0:
                try:
                    qualslist.append(bQ[bPos]())
                    success = True
                    break
                except:
                    bPos - 1
                if success == False:
                    qualslist.append(chr(2 + qual))
        elif val > t:
            read = read[:pos + 3] + 'G' + read[pos + 4:]
            bPos = pos - 1
            while bPos >= 0:
                try:
                    qualslist.append(bQ[bPos]())
                    success = True
                    break
                except:
                    bPos - 1
                if success == False:
                    qualslist.append(chr(2 + qual))
        elif val > a:
            read = read[:pos + 3] + 'T' + read[pos + 4:]
            bPos = pos - 1
            while bPos >= 0:
                try:
                    qualslist.append(bQ[bPos]())
                    success = True
                    break
                except:
                    bPos - 1
                if success == False:
                    qualslist.append(chr(2 + qual))
        else:
            read = read[:pos + 3] + 'A' + read[pos + 4:]
            bPos = pos - 1
            while bPos >= 0:
                try:
                    qualslist.append(bQ[bPos]())
                    success = True
                    break
                except:
                    bPos - 1
                if success == False:
                    qualslist.append(chr(2 + qual))
        if index in delD:
            delete = delD[index]()
            read = read[:pos + 4] + read[pos + delete + 4:]
        if index in insD:
            insert = insD[index]()
            read = read[:pos + 4] + insert + read[pos + 4:]
            for i in insert:
                iPos = pos - 1
                while iPos >= 0:
                    try:
                        qualslist.append(iQ[iPos]())
                        success = True
                        break
                    except:
                        iPos -= 1
                    if success == False:
                        qualslist.append(chr(2 + qual))
            pos += len(insert)
        pos += 1
    qualslist.append(qualslist[-1])
    readback = read
    read = read[4:readLen + 4]
    quals = ''.join(qualslist)[:readLen]
    if len(quals) != len(read):
        print("unexpected stop")
        print('read {}, quality {}'.format(read, quals))
        return None, None
    return read, quals


if __name__ == "__main__":
    main()
