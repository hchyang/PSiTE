#!/usr/bin/env python3

#########################################################################
# Author: Bingxin Lu
# Created Time: 2017-12-18
# File Name: fa2wes.py
# Description: Simulate WES reads from whole genome sequences
#########################################################################

import sys
import os
import argparse
import numpy
import logging
import pyfaidx
import subprocess
import shutil
import glob
import multiprocessing
import pip
from csite.phylovar import check_purity, check_seed, random_int
from csite.fa2wgs import check_folder, check_file, check_depth, merge_fq, OutputExistsError, tipnode_leaves_counting, genomesize

# handle the error below
# python | head == IOError: [Errno 32] Broken pipe
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE, SIG_DFL)


# MAX_READFRAC specifies the maximum fraction of simulated total short reads for a single run of simulation. If a large number of reads are simulated, this allows simulation in several batches, with each batch generating a smaller number of reads. Different batches can run at the same time. In the end, the output files from different batches are merged.
MAX_READFRAC = 0.02
MAX_CHROM = 512000000
RATIO_WESSIM = 0.85
RATIO_CAPGEM = 0.39


def check_normal_fa(normal_dir):
    '''
    There must be one fasta file for each haplotype in the normal dir
    '''
    for parental in 0, 1:
        if not os.path.isfile('{}/normal.parental_{}.fa'.format(normal_dir, parental)):
            raise argparse.ArgumentTypeError('Cannot find normal.parental_{}.fa under directory: {}'.format(
                parental, normal_dir))


def check_tumor_fa(tumor_dir, tip_node_leaves, simulator):
    '''
    Ensure the size of a chromsome is not too large for 'samtools index'
    See https://github.com/samtools/htsjdk/issues/447 for the issues for large chromosomes
    '''
    for tip_node, leaves in tip_node_leaves.items():
        for parental in 0, 1:
            fasta = '{}/{}.parental_{}.fa'.format(tumor_dir, tip_node, parental)
            if not os.path.isfile(fasta):
                raise argparse.ArgumentTypeError('Cannot find {}.parental_{}.fa under directory: {}'.format(
                    tip_node, parental, tumor_dir))
            if (simulator == 'capgem'):
                fa = pyfaidx.Faidx(fasta)
                for chroms in fa.index.keys():
                    chr_len = fa.index[chroms].rlen
                    if(chr_len > MAX_CHROM):
                        raise argparse.ArgumentTypeError('The size of chromsome {} ({}) for {} is larger than 512 M!'.format(
                            chroms, chr_len, fasta))


class TargetAction(argparse.Action):
    # adapted from documentation
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, values)
        defaultsim = getattr(namespace, 'simulator')
        if defaultsim == 'wessim' :
            defaultval = RATIO_WESSIM
        elif defaultsim == 'capgem':
            defaultval = RATIO_CAPGEM
        else:
            defaultval = 0.5
        setattr(namespace, 'ontarget_ratio', defaultval)


def check_program(value):
    if value == "capgem":
        progs=['bowtie2-build', 'bowtie2', 'samtools', 'capsim']
        for prog in progs:
            if shutil.which(prog) is None:
                raise argparse.ArgumentTypeError("Cannot find program '{}'. Please ensure that you have installed it!".format(prog))
    elif value == "wessim":
        progs=['samtools', 'faToTwoBit', 'blat']
        for prog in progs:
            if shutil.which(prog) is None:
                raise argparse.ArgumentTypeError("Cannot find program '{}'. Please ensure that you have installed it!".format(prog))
        package = "pysam"
        try:
            import package
        except:
            pip.main(['install', package])
    else:
        pass
    return value


def check_snakemake(value):
    # Use double quotes around option --cluster to distinguish with the single quotes around --snakemake
    value = value.replace("'",'"')
    return value


def compute_target_size(ftarget):
    '''
    Comute target size from the provided file
    '''
    size = 0
    with open(ftarget,'r') as fin:
        for line in fin:
            if line.startswith("#"):
                pass
            else:
                fields = line.split("\t")
                start = int(fields[1])
                end = int(fields[2])
                size += end - start
    return size


def compute_normal_gsize(normal_dir):
    normal_gsize = 0
    gsize_file = os.path.join(normal_dir, 'normal_gsize.txt')
    if os.path.exists(gsize_file):
        with open(gsize_file,'r') as fin:
            normal_gsize = int(fin.read().strip())
    else:
        for parental in 0, 1:
            normal_gsize += genomesize(
                fasta='{}/normal.parental_{}.fa'.format(normal_dir, parental))
        with open(gsize_file,'w') as fout:
            line = '{}\n'.format(normal_gsize)
            fout.write(line)

    return normal_gsize


def compute_tumor_dna(tumor_dir, tip_node_leaves):
    tumor_dna = 0
    tip_node_gsize = {}
    gsize_file = os.path.join(tumor_dir, 'tumor_gsize.txt')
    if os.path.exists(gsize_file):
        with open(gsize_file,'r') as fin:
            tumor_dna = int(fin.readline().strip())
            for line in fin:
                fields = line.strip().split('\t')
                tip_node_gsize[fields[0]] = [int(fields[1]), int(fields[2]), int(fields[3])]
    else:
        for tip_node, leaves in tip_node_leaves.items():
            # The value of tip_node_gsize[tip_node] is a list of three elements:
            # 0) genomesize of parental 0
            # 1) genomesize of parental 1
            # 2) the sum of parental 0 and 1
            tip_node_gsize[tip_node] = []

            for parental in 0, 1:
                assert os.path.isfile('{}/{}.parental_{}.fa'.format(tumor_dir, tip_node, parental)),\
                'Cannot find {}.parental_{}.fa under the tumor directory: {}'.format(
                    tip_node, parental, tumor_dir)
                tip_node_gsize[tip_node].append(genomesize(
                    fasta='{}/{}.parental_{}.fa'.format(tumor_dir, tip_node, parental)))

            tip_node_gsize[tip_node].append(
                tip_node_gsize[tip_node][0] + tip_node_gsize[tip_node][1])
            tumor_dna += tip_node_gsize[tip_node][2] * tip_node_leaves[tip_node]

        with open(gsize_file,'w') as fout:
            line = '{}\n'.format(tumor_dna)
            fout.write(line)
            for key, val in tip_node_gsize.items():
                line = '{}\t{}\t{}\t{}\n'.format(key, val[0], val[1], val[2])
                fout.write(line)

    return tip_node_gsize, tumor_dna


def merge_normal_sample(args, outdir):
    suffixes = ['fastq.gz', '1.fastq.gz', '2.fastq.gz']
    sample_fq_files = []
    for suffix in suffixes:
        prefix = '{}/{}_reads/normal.parental_[01]/normal_normal.parental_[01]*_'.format(outdir, args.simulator)
        source = glob.glob(prefix + suffix)
        if len(source):
            target_dir = os.path.join(outdir, 'merged')
            if not os.path.exists(target_dir):
                os.makedirs(target_dir)
            target = '{}/{}_normal_{}'.format(target_dir, args.simulator, suffix)
            source.sort()
            sample_fq_files.append([target, source])

    pool = multiprocessing.Pool(processes=args.cores)
    results = []
    for x in sample_fq_files:
        results.append(pool.apply_async(merge_fq, args=x))
    pool.close()
    pool.join()
    for result in results:
        result.get()


def merge_tumor_sample(args, tip_node_leaves, outdir):
    suffixes = ['fastq.gz', '1.fastq.gz', '2.fastq.gz']

    sample_fq_files = []
    for suffix in suffixes:
        if args.separate:
            for tip_node in sorted(tip_node_leaves.keys()) + ['normal']:
                prefix = '{}/{}_reads/{}.parental_[01]/tumor_{}.parental_[01]*_'.format(
                    outdir, args.simulator, tip_node, tip_node)
                source = glob.glob(prefix + suffix)
                if len(source):
                    target_dir = os.path.join(outdir, 'separate')
                    if not os.path.exists(target_dir):
                        os.makedirs(target_dir)
                    target = '{}/{}_{}'.format(target_dir, tip_node, suffix)
                    source.sort()
                    sample_fq_files.append([target, source])
        elif args.single:
            for tip_node in sorted(tip_node_leaves.keys()):
                prefix = '{}/{}_reads/{}.parental_[01]/tumor_{}.parental_[01]*_'.format(
                    outdir, args.simulator, tip_node, tip_node)
                source = glob.glob(prefix + suffix)
                if len(source):
                    target_dir = os.path.join(outdir, 'separate')
                    if not os.path.exists(target_dir):
                        os.makedirs(target_dir)
                    target = '{}/{}_{}'.format(target_dir, tip_node, suffix)
                    source.sort()
                    sample_fq_files.append([target, source])
        else:
            prefix = '{}/{}_reads/*.parental_[01]/tumor_*.parental_[01]*_'.format(outdir, args.simulator)
            source = glob.glob(prefix + suffix)
            if len(source):
                target_dir = os.path.join(outdir, 'merged')
                if not os.path.exists(target_dir):
                    os.makedirs(target_dir)
                target = '{}/{}_tumor_{}'.format(target_dir, args.simulator, suffix)
                source.sort()
                sample_fq_files.append([target, source])

    pool = multiprocessing.Pool(processes=args.cores)
    results = []
    for x in sample_fq_files:
        results.append(pool.apply_async(merge_fq, args=x))
    pool.close()
    pool.join()
    for result in results:
        result.get()



def clean_output(level, outdir):
    '''
    Remove intermediate output of WES simulators according to the specified levels.
    Level 0: keep all the files.
    Level 1: keep files that are necessary for rerunning simulation ('config', 'genome_index', 'mapping', 'merged', and 'separate').
    Level 2: keep only final results ('merged' and 'separate').
    '''
    if level == 0:
        return
    elif level == 1:
        # Used to rerun based on previous mapping results
        dirs_keep = ['config', 'genome_index', 'mapping', 'merged', 'separate']
        for entry in os.scandir(outdir):
            if entry.is_dir():
                if entry.name not in dirs_keep:
                    shutil.rmtree(entry.path)
    elif level == 2:
        # Only keep the final reads
        dirs_keep = ['merged', 'separate']
        for entry in os.scandir(outdir):
            if entry.is_dir():
                if entry.name not in dirs_keep:
                    shutil.rmtree(entry.path)


def prepare_sample_normal(sample_file, args, normal_gsize, target_size):
    '''
    Create a configuration file for simulating normal samples with snakemake
    '''
    if args.single_end:
        rlen = args.rlen
    else:
        rlen = args.rlen * 2
    if args.normal_rdepth > 0:
        total_rnum = int((args.normal_rdepth * target_size) / (rlen * args.ontarget_ratio))
    else:
        total_rnum = args.normal_rnum
    logging.info(' Total number of reads to simulate for normal sample: %d', total_rnum)
    MAX_READNUM = int(total_rnum * MAX_READFRAC)

    with open(sample_file, 'w') as fout:
        fout.write('probe: {}\n'.format(os.path.abspath(args.probe)))
        fout.write('error_model: {}\n'.format(os.path.abspath(args.error_model)))
        fout.write('directory: normal\n')

        fout.write('genomes:\n')
        # two normal cell haplotypes
        for parental in 0, 1:
            ref = '{}/normal.parental_{}.fa'.format(args.normal, parental)
            fullname = os.path.abspath(ref)
            fout.write('  normal.parental_{}: {}\n'.format(parental, fullname))

        fout.write('samples:\n')
        total_num_splits = 0
        # two normal cell haplotypes
        for parental in 0, 1:
            ref = '{}/normal.parental_{}.fa'.format(args.normal, parental)
            proportion = genomesize(fasta=ref) / normal_gsize
            readnum = int(proportion * total_rnum)
            if readnum > MAX_READNUM:
                num_splits = int(numpy.ceil(readnum / MAX_READNUM))
                total_num_splits += num_splits
                for split in range(1, num_splits+1):
                    fout.write('  normal_normal.parental_{}_{}:\n'.format(parental, str(split)))
                    fout.write('    gid: normal.parental_{}\n'.format(parental))
                    fout.write('    proportion: {}\n'.format(str(proportion/num_splits)))
                    fout.write('    split: {}\n'.format(str(split)))
                    split_readnum = int(numpy.ceil(readnum/num_splits))
                    fout.write('    readnum: {}\n'.format(str(split_readnum)))
                    seed = random_int()
                    fout.write('    seed: {}\n'.format(str(seed)))
            else:
                total_num_splits += 1
                fout.write('  normal_normal.parental_{}:\n'.format(parental))
                fout.write('    gid: normal.parental_{}\n'.format(parental))
                fout.write('    proportion: {}\n'.format(str(proportion)))
                fout.write('    readnum: {}\n'.format(str(readnum)))
                seed = random_int()
                fout.write('    seed: {}\n'.format(str(seed)))

    return total_num_splits


def prepare_sample_tumor(sample_file, args, total_cells, normal_cells, normal_gsize, tip_node_leaves, tip_node_gsize, target_size):
    '''
    Create a configuration file for simulating tumor samples with snakemake
    '''
    if args.single_end:
        rlen = args.rlen
    else:
        rlen = args.rlen * 2
    if args.tumor_rdepth > 0:
        total_rnum = int((args.tumor_rdepth * target_size) / (rlen * args.ontarget_ratio))
    else:
        total_rnum = args.tumor_rnum
    logging.info(' Total number of reads to simulate for tumor sample: %d', total_rnum)
    MAX_READNUM = int(total_rnum * MAX_READFRAC)

    with open(sample_file, 'w') as fout:
        fout.write('probe: {}\n'.format(os.path.abspath(args.probe)))
        fout.write('error_model: {}\n'.format(os.path.abspath(args.error_model)))
        fout.write('directory: tumor\n')
        fout.write('genomes:\n')
        # two normal cell haplotypes
        if not args.single:
            for parental in 0, 1:
                ref = '{}/normal.parental_{}.fa'.format(args.normal, parental)
                fullname = os.path.abspath(ref)
                fout.write('  normal.parental_{}: {}\n'.format(parental, fullname))
        # tumor cells haplotypes
        for tip_node in sorted(tip_node_leaves.keys()):
            for parental in 0, 1:
                ref = '{}/{}.parental_{}.fa'.format(
                    args.tumor, tip_node, parental)
                fullname = os.path.abspath(ref)
                fout.write('  {}.parental_{}: {}\n'.format(tip_node, parental, fullname))

        fout.write('samples:\n')
        total_num_splits = 0
        # two normal cell haplotypes, only generated under non-single mode
        if not args.single:
            for parental in 0, 1:
                ref = '{}/normal.parental_{}.fa'.format(args.normal, parental)
                fullname = os.path.abspath(ref)
                cell_proportion = normal_cells / total_cells
                proportion = cell_proportion * genomesize(fasta=ref) / normal_gsize
                readnum = int(proportion * total_rnum)
                if readnum > MAX_READNUM:
                    num_splits = int(numpy.ceil(readnum / MAX_READNUM))
                    total_num_splits += num_splits
                    for split in range(1, num_splits+1):
                        fout.write('  tumor_normal.parental_{}_{}:\n'.format(parental, str(split)))
                        fout.write('    gid: normal.parental_{}\n'.format(parental))
                        fout.write('    cell_proportion: {}\n'.format(str(cell_proportion)))
                        fout.write('    proportion: {}\n'.format(str(proportion/num_splits)))
                        fout.write('    split: {}\n'.format(str(split)))
                        split_readnum = int(numpy.ceil(readnum/num_splits))
                        fout.write('    readnum: {}\n'.format(str(split_readnum)))
                        seed = random_int()
                        fout.write('    seed: {}\n'.format(str(seed)))
                else:
                    total_num_splits += 1
                    fout.write('  tumor_normal.parental_{}:\n'.format(parental))
                    fout.write('    gid: normal.parental_{}\n'.format(parental))
                    fout.write('    cell_proportion: {}\n'.format(str(cell_proportion)))
                    fout.write('    proportion: {}\n'.format(str(proportion)))
                    fout.write('    readnum: {}\n'.format(str(readnum)))
                    seed = random_int()
                    fout.write('    seed: {}\n'.format(str(seed)))

        # tumor cells haplotypes
        for tip_node in sorted(tip_node_leaves.keys()):
            for parental in 0, 1:
                ref = '{}/{}.parental_{}.fa'.format(
                    args.tumor, tip_node, parental)
                fullname = os.path.abspath(ref)
                if args.single:
                    cell_proportion = 1
                else:
                    cell_proportion = tip_node_leaves[tip_node] / total_cells
                proportion = cell_proportion * tip_node_gsize[tip_node][parental] / tip_node_gsize[tip_node][2]
                readnum = int(proportion * total_rnum)
                if readnum > MAX_READNUM:
                    num_splits = int(numpy.ceil(readnum / MAX_READNUM))
                    total_num_splits += num_splits
                    for split in range(1, num_splits+1):
                        fout.write('  tumor_{}.parental_{}_{}:\n'.format(tip_node, parental, str(split)))
                        fout.write('    gid: {}.parental_{}\n'.format(tip_node, parental))
                        fout.write('    proportion: {}\n'.format(str(proportion/num_splits)))
                        fout.write('    split: {}\n'.format(str(split)))
                        split_readnum = int(numpy.ceil(readnum/num_splits))
                        fout.write('    readnum: {}\n'.format(str(split_readnum)))
                        seed = random_int()
                        fout.write('    seed: {}\n'.format(str(seed)))
                else:
                    total_num_splits += 1
                    fout.write('  tumor_{}.parental_{}:\n'.format(tip_node, parental))
                    fout.write('    gid: {}.parental_{}\n'.format(tip_node, parental))
                    fout.write('    cell_proportion: {}\n'.format(str(cell_proportion)))
                    fout.write('    proportion: {}\n'.format(str(proportion)))
                    fout.write('    readnum: {}\n'.format(str(readnum)))
                    seed = random_int()
                    fout.write('    seed: {}\n'.format(str(seed)))
    return total_num_splits


def prepare_sample_all(sample_file, args, total_cells, normal_cells, normal_gsize, tip_node_leaves, tip_node_gsize, target_size):
    '''
    Create a configuration file for simulating both tumor and normal samples with snakemake
    '''
    if args.single_end:
        rlen = args.rlen
    else:
        rlen = args.rlen * 2
    if args.normal_rdepth > 0:
        total_rnum_normal = int((args.normal_rdepth * target_size) / (rlen * args.ontarget_ratio))
    else:
        total_rnum_normal = args.normal_rnum
    logging.info(' Total number of reads to simulate for normal sample: %d', total_rnum_normal)
    MAX_READNUM_NORMAL = int(total_rnum_normal * MAX_READFRAC)
    if args.tumor_rdepth > 0:
        total_rnum = int((args.tumor_rdepth * target_size) / (rlen * args.ontarget_ratio))
    else:
        total_rnum = args.tumor_rnum
    logging.info(' Total number of reads to simulate for tumor sample: %d', total_rnum)
    MAX_READNUM = int(total_rnum * MAX_READFRAC)

    with open(sample_file, 'w') as fout:
        fout.write('probe: {}\n'.format(os.path.abspath(args.probe)))
        fout.write('error_model: {}\n'.format(os.path.abspath(args.error_model)))
        fout.write('genomes:\n')
        # two normal cell haplotypes
        if not args.single:
            for parental in 0, 1:
                ref = '{}/normal.parental_{}.fa'.format(args.normal, parental)
                fullname = os.path.abspath(ref)
                fout.write('  normal.parental_{}: {}\n'.format(parental, fullname))
        # tumor cells haplotypes
        for tip_node in sorted(tip_node_leaves.keys()):
            for parental in 0, 1:
                ref = '{}/{}.parental_{}.fa'.format(
                    args.tumor, tip_node, parental)
                fullname = os.path.abspath(ref)
                fout.write('  {}.parental_{}: {}\n'.format(tip_node, parental, fullname))

        fout.write('samples:\n')
        # Use gid to distinguish normal genomes in normal sample and tumor sample
        total_num_splits = 0

        # two normal cell haplotypes for normal samples
        for parental in 0, 1:
            ref = '{}/normal.parental_{}.fa'.format(args.normal, parental)
            proportion = genomesize(fasta=ref) / normal_gsize
            readnum = int(proportion * total_rnum_normal)
            if readnum > MAX_READNUM_NORMAL:
                num_splits = int(numpy.ceil(readnum / MAX_READNUM_NORMAL))
                total_num_splits += num_splits
                for split in range(1, num_splits+1):
                    fout.write('  normal_normal.parental_{}_{}:\n'.format(parental, str(split)))
                    fout.write('    gid: normal.parental_{}\n'.format(parental))
                    fout.write('    proportion: {}\n'.format(str(proportion/num_splits)))
                    fout.write('    split: {}\n'.format(str(split)))
                    split_readnum = int(numpy.ceil(readnum/num_splits))
                    fout.write('    readnum: {}\n'.format(str(split_readnum)))
                    seed = random_int()
                    fout.write('    seed: {}\n'.format(str(seed)))
            else:
                total_num_splits += 1
                fout.write('  normal_normal.parental_{}:\n'.format(parental))
                fout.write('    gid: normal.parental_{}\n'.format(parental))
                fout.write('    proportion: {}\n'.format(str(proportion)))
                fout.write('    readnum: {}\n'.format(str(readnum)))
                seed = random_int()
                fout.write('    seed: {}\n'.format(str(seed)))

        # two normal cell haplotypes for tumor sample, only generated under non-single mode
        if not args.single:
            for parental in 0, 1:
                ref = '{}/normal.parental_{}.fa'.format(args.normal, parental)
                fullname = os.path.abspath(ref)
                cell_proportion = normal_cells / total_cells
                proportion = cell_proportion * genomesize(fasta=ref) / normal_gsize
                readnum = int(proportion * total_rnum)
                if readnum > MAX_READNUM:
                    num_splits = int(numpy.ceil(readnum / MAX_READNUM))
                    total_num_splits += num_splits
                    for split in range(1, num_splits+1):
                        fout.write('  tumor_normal.parental_{}_{}:\n'.format(parental, str(split)))
                        fout.write('    gid: normal.parental_{}\n'.format(parental))
                        fout.write('    cell_proportion: {}\n'.format(str(cell_proportion)))
                        fout.write('    proportion: {}\n'.format(str(proportion/num_splits)))
                        fout.write('    split: {}\n'.format(str(split)))
                        split_readnum = int(numpy.ceil(readnum/num_splits))
                        fout.write('    readnum: {}\n'.format(str(split_readnum)))
                        seed = random_int()
                        fout.write('    seed: {}\n'.format(str(seed)))
                else:
                    total_num_splits += 1
                    fout.write('  tumor_normal.parental_{}:\n'.format(parental))
                    fout.write('    gid: normal.parental_{}\n'.format(parental))
                    fout.write('    cell_proportion: {}\n'.format(str(cell_proportion)))
                    fout.write('    proportion: {}\n'.format(str(proportion)))
                    fout.write('    readnum: {}\n'.format(str(readnum)))
                    seed = random_int()
                    fout.write('    seed: {}\n'.format(str(seed)))

        # tumor cells haplotypes
        for tip_node in sorted(tip_node_leaves.keys()):
            for parental in 0, 1:
                ref = '{}/{}.parental_{}.fa'.format(
                    args.tumor, tip_node, parental)
                fullname = os.path.abspath(ref)
                if args.single:
                    cell_proportion = 1
                else:
                    cell_proportion = tip_node_leaves[tip_node] / total_cells
                proportion = cell_proportion * tip_node_gsize[tip_node][parental] / tip_node_gsize[tip_node][2]
                readnum = int(proportion * total_rnum)
                if readnum > MAX_READNUM:
                    num_splits = int(numpy.ceil(readnum / MAX_READNUM))
                    total_num_splits += num_splits
                    for split in range(1, num_splits+1):
                        fout.write('  tumor_{}.parental_{}_{}:\n'.format(tip_node, parental, str(split)))
                        fout.write('    gid: {}.parental_{}\n'.format(tip_node, parental))
                        fout.write('    proportion: {}\n'.format(str(proportion/num_splits)))
                        fout.write('    split: {}\n'.format(str(split)))
                        split_readnum = int(numpy.ceil(readnum/num_splits))
                        fout.write('    readnum: {}\n'.format(str(split_readnum)))
                        seed = random_int()
                        fout.write('    seed: {}\n'.format(str(seed)))
                else:
                    total_num_splits += 1
                    fout.write('  tumor_{}.parental_{}:\n'.format(tip_node, parental))
                    fout.write('    gid: {}.parental_{}\n'.format(tip_node, parental))
                    fout.write('    cell_proportion: {}\n'.format(str(cell_proportion)))
                    fout.write('    proportion: {}\n'.format(str(proportion)))
                    fout.write('    readnum: {}\n'.format(str(readnum)))
                    seed = random_int()
                    fout.write('    seed: {}\n'.format(str(seed)))
    return total_num_splits



def run_snakemake(outdir, args, sample_file, snake_file):
    # Copy sample Snakefile to the output directory
    if args.simulator == 'capsim':
        snake_file_copy = os.path.join(outdir, 'config/Snakefile_capsim')
        shutil.copyfile(snake_file, snake_file_copy)
    elif args.simulator == 'wessim':
        snake_file_copy = os.path.join(outdir, 'config/Snakefile_wessim')
        shutil.copyfile(snake_file, snake_file_copy)
        tmp = os.path.join(outdir, 'config/tmp')
        if args.single_end:
            # remove option -p
            with open(snake_file_copy,'r') as fin, open(tmp,'w') as fout:
                for line in fin:
                    if 'Wessim2' in line:
                        nline=line.replace("-p","")
                        fout.write(nline)
                    else:
                        fout.write(line)
            os.remove(snake_file_copy)
            shutil.move(tmp, snake_file_copy)
    else:
        snake_file_copy = os.path.join(outdir, 'config/Snakefile_capgem')
        shutil.copyfile(snake_file, snake_file_copy)
        tmp = os.path.join(outdir, 'config/tmp')
        if args.single_end:
            # remove option -p
            with open(snake_file_copy,'r') as fin, open(tmp,'w') as fout:
                for line in fin:
                    if 'frag2read' in line:
                        nline=line.replace("-p","")
                        fout.write(nline)
                    else:
                        fout.write(line)
            os.remove(snake_file_copy)
            shutil.move(tmp, snake_file_copy)

    # Remove the surrounding quotes around snakemake command.
    orig_params = args.snakemake.split()
    if not ('--cores' in args.snakemake or '--jobs' in args.snakemake or '-j' in args.snakemake):
        # Use the number of cores specified here
        orig_params += ['-j', str(args.cores)]
    if '--cluster' in args.snakemake and '--cluster-config' not in args.snakemake:
        cluster_file = os.path.join(os.path.dirname(os.path.realpath(sys.argv[0])), 'wes/config/cluster.yaml')
        assert os.path.isfile(cluster_file), 'Cannot find cluster.yaml below under the program directory:{}'.format(cluster_file)
        cluster_file_copy = os.path.join(outdir, 'config/cluster.yaml')
        shutil.copyfile(cluster_file, cluster_file_copy)
        orig_params += ['--cluster-config', cluster_file_copy]

    config = ' rlen=' + str(args.rlen)
    final_cmd_params =  orig_params + ['-s', os.path.abspath(snake_file_copy), '-d', os.path.abspath(outdir), '--configfile', os.path.abspath(sample_file), '--config', config]
    logging.info(' Command: %s', ' '.join(final_cmd_params))

    os.system(' '.join(final_cmd_params))


def main(progname=None):
    parser = argparse.ArgumentParser(
        description='a wrapper of simulating targeted capture sequencing from reference genome files',
        prog=progname if progname else sys.argv[0])

    group1 = parser.add_argument_group('Input arguments')
    group1.add_argument('-n', '--normal', metavar='DIR', type=check_folder, required=True,
                       help='The directory of the fasta files of normal genomes')
    group1.add_argument('-t', '--tumor', metavar='DIR', type=check_folder, required=True,
                       help='The directory of the fasta files of tumor genomes')
    group1.add_argument('-m','--map', metavar='FILE', type=check_file, required=True,
                       help='The map file containing the relationship between tip nodes and samples')
    group1.add_argument( '--probe', metavar='FILE', type=check_file, required=True,
                       help='The Probe file containing the probe sequences (FASTA format)')
    group1.add_argument('--target', metavar='FILE', type=check_file, required=True,
                       help='The Target file containing the target regions (BED format)')

    group2 = parser.add_argument_group('Arguments for simulation')
    default = 0.6
    group2.add_argument('-p', '--purity', metavar='FLOAT', type=check_purity, default=default,
                       help='The proportion of tumor cells in simulated sample [{}]'.format(default))
    default = 150
    group2.add_argument('--rlen', metavar='INT', type=int, default=default,
                       help='Illumina: read length [{}]'.format(default))
    group2.add_argument('--single_end', action='store_true',
        help='Simulating single-end reads')
    group = group2.add_mutually_exclusive_group()
    default = 0
    group.add_argument('-d', '--tumor_rdepth', metavar='FLOAT', type=check_depth, default=default,
                       help='The mean depth of tumor sample for simulating short reads [{}]'.format(default))
    default = 0
    group.add_argument('-r', '--tumor_rnum', metavar='INT', type=int, default=default,
                       help='The number of short reads to simulate for tumor sample [{}]'.format(default))
    group = group2.add_mutually_exclusive_group()
    default = 0
    group.add_argument('-D', '--normal_rdepth', metavar='FLOAT', type=check_depth, default=default,
                       help='The mean depth of normal sample for simulating short reads [{}]'.format(default))
    default = 0
    group.add_argument('-R', '--normal_rnum', metavar='INT', type=int, default=default,
                       help='The number of short reads to simulate for normal sample [{}]'.format(default))
    default = None
    group2.add_argument('--random_seed', metavar='INT', type=check_seed,
                       help='The seed for random number generator [{}]'.format(default))
    default = 'wessim'
    group2.add_argument('--simulator', default=default, choices=['wessim','capgem'], action=TargetAction, type=check_program,
                       help='The whole-exome sequencing simulator used for simulating short reads [{}]'.format(default))
    default = RATIO_WESSIM
    group2.add_argument('--ontarget_ratio', metavar='FLOAT', type=float, default=default,
                       help='The percentage that simulated reads are expected to be from the target regions. It is dependent on the simulator. The default value is {} for wessim and {} for capgem [{}]'.format(RATIO_WESSIM, RATIO_CAPGEM, default))
    default = None
    group2.add_argument('--error_model', metavar='FILE', type=check_file,
                       help='The file containing the empirical error model for NGS reads generated by GemErr (It must be provided when capgem or wessim is used for simulation) [{}]'.format(default))
    group2.add_argument('--single', action='store_true',
        help='single cell mode. After this setting, the value of --rnum is for each tumor cell (not the whole tumor sample anymore)')
    default = "snakemake --rerun-incomplete -k --latency-wait 120"
    group2.add_argument('--snakemake', metavar='STR', type=check_snakemake, default=default,
                       help="The snakemake command used for calling a whole-exome sequencing simulator. The Snakefile for a simulator is under the directory 'wes/config' of the source code. Additional parameters for a simulator can be adjusted in the Snakefile ['{}']".format(default))
    default = 1
    group2.add_argument('--cores', type=int, default=default, metavar='INT',
                        help="The number of cores used to run the program (including snakemake). If '--cores' or '--jobs' or '-j' is specified in the options of snakemake, the value specified by '--cores' here will be ignored when snakemake is called [{}]".format(default))

    group3 = parser.add_argument_group('Output arguments')
    default = 'wes_reads'
    group3.add_argument('-o', '--output', metavar='DIR', type=str, default=default,
                       help='The output directory [{}]'.format(default))
    default = 'fa2wes.log'
    group3.add_argument('-g', '--log', metavar='FILE', type=str, default=default,
                       help='The log file to save the settings of each command [{}]'.format(default))
    default = 2
    group3.add_argument('--out_level', type=int, choices=[0, 1, 2], default=default,
                       help="The level used to indicate how many intermediate output files are kept. \
                       Level 0: keep all the files.\
                       Level 1: keep files that are necessary for rerunning simulation ('config', 'genome_index', 'mapping', 'merged', and 'separate'). \
                       Level 2: keep only final results ('merged' and 'separate') [{}]".format(default))
    group3.add_argument('--separate', action='store_true',
                        help='Output the reads of each genome separately')

    args = parser.parse_args()
    check_normal_fa(args.normal)

    # logging and random seed setting
    logging.basicConfig(filename=args.log,
                        filemode='w', format='[%(asctime)s] %(levelname)s: %(message)s',
                        datefmt='%m-%d %H:%M:%S', level='INFO')
    argv_copy = sys.argv[:]
    try:
        snakemake_index = argv_copy.index('--snakemake')
        # Single quotes are required for the snakemake command
        snakemake_str = argv_copy[snakemake_index + 1]
        if "'" in snakemake_str:
            snakemake_str = snakemake_str.replace("'",'"')
        argv_copy[snakemake_index + 1] = "'{}'".format(snakemake_str)
    except ValueError:
        pass
    argv_copy.insert(1, 'fa2wes')
    logging.info(' Command: %s', ' '.join(argv_copy))

    if args.random_seed == None:
        seed = random_int()
    else:
        seed = args.random_seed
    logging.info(' Ontarget ratio: %s', str(args.ontarget_ratio))
    logging.info(' Random seed: %d', seed)
    numpy.random.seed(seed)

    #create output folders
    if os.path.exists(args.output):
        if os.path.isdir(args.output):
            pass
        else:
            raise OutputExistsError("A file in the name of '{}' exists.\nDelete it or try another name as output folder.".format(args.output))
    else:
        os.makedirs(args.output,mode=0o755)

    wes_dir = os.path.join(os.path.dirname(os.path.realpath(sys.argv[0])), 'wes')
    # Add path variables
    if args.simulator == 'capsim':  # Not exposed to user for simplificity
        snake_file = os.path.join(wes_dir, 'config/Snakefile_capsim')
    elif args.simulator == 'wessim':
        snake_file = os.path.join(wes_dir, 'config/Snakefile_wessim')
        wessim_dir = os.path.join(wes_dir, 'wessim')
        os.environ['PATH'] += os.pathsep + wessim_dir
    else: # capgem
        snake_file = os.path.join(wes_dir, 'config/Snakefile_capgem')
        capgem_dir = os.path.join(wes_dir, 'capgem')
        if os.path.exists(os.path.join(capgem_dir, 'bin')):
            os.environ['PATH'] += os.pathsep + os.path.join(capgem_dir, 'bin')
        os.environ['PATH'] += os.pathsep + os.path.join(capgem_dir, 'src')
    assert os.path.isfile(snake_file), 'Cannot find Snakefile {} under the program directory:\n'.format(snake_file)

    normal_gsize = compute_normal_gsize(args.normal)
    target_size = compute_target_size(args.target)
    logging.info(' Size of target region: %s', str(target_size))

    # Simulate normal and tumor sample at the same time
    if (args.tumor_rdepth > 0 or args.tumor_rnum > 0) and (args.normal_rdepth > 0 or args.normal_rnum > 0):
        tip_node_leaves = tip_node_leaves_counting(f=args.map)
        check_tumor_fa(args.tumor, tip_node_leaves, args.simulator)

        outdir = args.output
        configdir = os.path.join(outdir, 'config')
        if not os.path.exists(configdir):
            os.makedirs(configdir)

        if args.single:
            for tip_node in tip_node_leaves:
                assert tip_node_leaves[tip_node]==1,\
                    'In single mode, each tip node should represent only one cell.\n'+\
                    'But {} leaves are found underneath tip node {} in your map file!'.format(tip_node_leaves[tip_node], tip_node)

        tumor_cells = sum(tip_node_leaves.values())
        total_cells = tumor_cells / args.purity
        logging.info(' Number of total cells: %d', total_cells)
        normal_cells = total_cells - tumor_cells
        logging.info(' Number of normal cells: %d', normal_cells)
        normal_dna = normal_gsize * normal_cells
        tip_node_gsize, tumor_dna= compute_tumor_dna(args.tumor, tip_node_leaves)
        total_dna = (normal_dna + tumor_dna)

        sample_file = os.path.join(outdir, 'config/sample.yaml')
        total_num_splits = prepare_sample_all(sample_file, args, total_cells, normal_cells, normal_gsize, tip_node_leaves, tip_node_gsize, target_size)
        logging.info(' Number of splits in simulation: %d', total_num_splits)

        run_snakemake(outdir, args, sample_file, snake_file)
        merge_normal_sample(args, outdir)
        merge_tumor_sample(args, tip_node_leaves, outdir)
        clean_output(args.out_level, outdir)

    # Separate the simulation of tumor and normal samples
    elif args.tumor_rdepth > 0 or args.tumor_rnum > 0:
        tip_node_leaves = tip_node_leaves_counting(f=args.map)
        check_tumor_fa(args.tumor, tip_node_leaves, args.simulator)

        outdir = os.path.join(args.output, 'tumor')
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        configdir = os.path.join(outdir, 'config')
        if not os.path.exists(configdir):
            os.makedirs(configdir)

        if args.single:
            for tip_node in tip_node_leaves:
                assert tip_node_leaves[tip_node]==1,\
                    'In single mode, each tip node should represent only one cell.\n'+\
                    'But {} leaves are found underneath tip node {} in your map file!'.format(tip_node_leaves[tip_node], tip_node)

        tumor_cells = sum(tip_node_leaves.values())
        total_cells = tumor_cells / args.purity
        logging.info(' Number of total cells: %d', total_cells)
        normal_cells = total_cells - tumor_cells
        logging.info(' Number of normal cells: %d', normal_cells)
        normal_dna = normal_gsize * normal_cells
        tip_node_gsize, tumor_dna= compute_tumor_dna(args.tumor, tip_node_leaves)
        total_dna = (normal_dna + tumor_dna)

        sample_file = os.path.join(outdir, 'config/sample.yaml')
        total_num_splits = prepare_sample_tumor(sample_file, args, total_cells, normal_cells, normal_gsize, tip_node_leaves, tip_node_gsize, target_size)
        logging.info(' Number of splits in simulation: %d', total_num_splits)

        run_snakemake(outdir, args, sample_file, snake_file)
        merge_tumor_sample(args, tip_node_leaves, outdir)
        clean_output(args.out_level, outdir)

    elif args.normal_rdepth > 0 or args.normal_rnum > 0:
        outdir = os.path.join(args.output, 'normal')
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        configdir = os.path.join(outdir, 'config')
        if not os.path.exists(configdir):
            os.makedirs(configdir)

        sample_file = os.path.join(outdir, 'config/sample.yaml')
        total_num_splits = prepare_sample_normal(sample_file, args, normal_gsize, target_size)
        logging.info(' Number of splits in simulation: %d', total_num_splits)

        run_snakemake(outdir, args, sample_file, snake_file)
        merge_normal_sample(args, outdir)
        clean_output(args.out_level, outdir)
    else:
        logging.info('Please specify sequening depth!')
