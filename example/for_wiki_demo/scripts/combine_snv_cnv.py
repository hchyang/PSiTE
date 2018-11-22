#!/usr/bin/env python

#####################################################
# Author: Bingxin Lu
# Description: This script is used to combine the output of Mutect and Sequenza, which serves as the input to PyClone
# Usage: Run `python combine_snv_cnv.py -h` for details
#####################################################


from __future__ import division
import argparse
import gzip

def read_cnv(infile, is_dpratio):
    '''
    Read output of Sequenza
    '''
    cnvs=[]
    with open(infile,'r') as fin:
        fin.readline() # Skip header
        for line in fin:
            fields = line.split('\t')
            chrom = fields[0]
            chrom = chrom.replace('"', '')
            start = int(fields[1])
            end = int(fields[2])
            if is_dpratio:
                dpratio = float(fields[6])
                minor_cn = 0
                major_cn = int(round(dpratio, 0))
                # print('{}\t{}'.format(dpratio, major_cn))
            else:
                minor_cn = fields[11]
                major_cn = fields[10]
            # ignore mutations with copy number 0
            if int(minor_cn) > 0 or int(major_cn) > 0:
                record = (chrom, start, end, minor_cn, major_cn)
                cnvs.append(record)
    # print(cnvs)
    return cnvs

def read_type(infile, columns):
    '''
    Read types of samples.
    Format of input file: tab delimited with header -- Sample,Normal.vcf,Tumor,Tumor.vcf,Tumor.vcf.column
    columns: colomn ID of normal and tumor samples
    '''
    normal_samples=set()
    tumor_samples=set()
    with open(infile,'r') as fin:
        fin.readline() # Skip header
        for line in fin:
            fields = line.strip().split('\t')
            idx = int(columns[0]) - 1
            normal = fields[idx]
            normal_samples.add(normal)
            idx = int(columns[1]) - 1
            tumor = fields[idx]
            tumor_samples.add(tumor)
    return (normal_samples, tumor_samples)

def combine_snv_cnv(snv_file, cnvs, sampleID, min_CellFreq, max_PM, min_depth, outfile, is_gzipped, normal_samples, tumor_samples):
    tumor_index = 10    # By default tumor data is at the last column
    tumors = range(1, 10)
    tumor_types = []
    for t in tumors:
        tumor_types.append('0'+str(t))
    # print(tumor_types)

    if is_gzipped:
        fin = gzip.open(snv_file,'rt')
    else:
        fin = open(snv_file,'r')
    fout = open(outfile,'w')

    header = '\t'.join(["mutation_id","ref_counts","var_counts","normal_cn","minor_cn","major_cn"])
    fout.write(header)
    fout.write('\n')

    for line in fin:
        if line.startswith('##'):
            continue
        if line.startswith('#'):
            headers = line.strip().split('\t')
            type1 = headers[-2]
            type2 = headers[-1]
            print('Type at the second last column {}, Type at the last column {}'.format(type1, type2))
            if tumor_samples is not None and type1 in tumor_samples:
                tumor_index = 9
                if normal_samples is not None:
                    assert (type2 in normal_samples)
            if 'tumor' in type1 or 'tumour' in type1 or 'normal' in type2:
                tumor_index = 9
            # https://wiki.nci.nih.gov/display/TCGA/TCGA+barcode
            if 'TCGA' in type1:
                codes1 = type1.split('-')
                sample = codes1[3]
                # print(sample)
                if sample[0:-1] in tumor_types:
                    tumor_index = 9
            # print(tumor_index)
            continue
        fields = line.strip().split('\t')
        chrom = fields[0]
        if chrom == 'X' or chrom == 'Y' or 'M' in chrom:
            continue
        pos = int(fields[1])
        # print(pos)
        # Find whether this position overlaps with a CNV. If yes, output this mutation.
        for (cnv_chrom, start, end, minor_cn, major_cn) in cnvs:
            if chrom==cnv_chrom and pos>=start and pos<=end:
                # print(pos, start, end)
                # Extract "ref_counts","var_counts" from this line
                # mutation_id = 's{}_c{}_p{}'.format(sampleID, int(chrom), pos)
                # https://groups.google.com/forum/#!searchin/pyclone-user-group/binomial|sort:date/pyclone-user-group/gyH5nxJfyVw/nbdibizcBAAJ
                mutation_id = 'chr{}:{}'.format(int(chrom), pos)
                # Select the genotye for tumor
                genotye = fields[tumor_index]
                # print(genotye)
                allele_depth = genotye.split(':')[1].split(',')
                ref_counts = int(allele_depth[0])
                var_counts = int(allele_depth[1])
                cov = ref_counts + var_counts
                if cov < min_depth:
                    print(" SNV(s) excluded due to the number of reads covering the loci ({}) is below {}.".format(cov, min_depth));
                    continue
                freq = var_counts * 100 / (var_counts + ref_counts)
                total_cn = int(minor_cn) + int(major_cn)
                if total_cn > max_PM:
                    print(" SNV(s) excluded due to high-level amplifications (>{} copies) within that region. Consider increasing value of parameter max_PM to include these SNVs, provided high coverage data (> 150 fold) is available".format(max_PM));
                    continue
                if freq * total_cn < min_CellFreq*100:
                    print(" SNV(s) excluded due to AF*CN below {} (SNV can't be explained by an SP present in {}% or more of the sample).".format(min_CellFreq, min_CellFreq*100));
                    continue
                # "mutation_id","ref_counts","var_counts","normal_cn","minor_cn","major_cn"
                normal_cn = 2
                oline = '{}\t{}\t{}\t{}\t{}\t{}\n'.format(mutation_id, ref_counts, var_counts, normal_cn, minor_cn, major_cn)
                fout.write(oline)
                break
    if tumor_index == 10:
        print("Tumor column is the last column")
    else:
        print("Tumor column is the second last column")
    fin.close()
    fout.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", dest="snv_file", default="", required=True, help="The main output of Mutect (e.g. mutect.PASS.vcf)")
    parser.add_argument("-c", dest="cnv_file", default="", required=True, help="The main output of Sequenza (e.g. output_segments.txt)")
    parser.add_argument("-s", dest="sampleID", default='tumor', required=True, help="The ID of sample")
    parser.add_argument("-o", dest="output", default='', required=True, help="The output file")
    parser.add_argument("-t", dest="type_file", default="", help="The file specifying the type of each sample (normal or tumor)")
    parser.add_argument("-I", dest="columns", default="2,4", help="The columns corresponding to normal samples and tumor samples in the file specified by -t")
    parser.add_argument("-g", dest="is_gzipped", action="store_true", default=False, help="The SNV file is gzipped or not")
    parser.add_argument("-r", dest="is_dpratio", action="store_true", default=False, help="Using the estimates of depth ratio from Sequenza")
    parser.add_argument("-f", dest="min_CellFreq", type=float, default=0.1, help="The threshold to ignore low-frequency variants")
    parser.add_argument("-p", dest="max_PM", type=int, default=6, help="The threshold to ignore variants with high copy number")
    parser.add_argument("-d", dest="min_depth", type=int, default=0, help="The threshold to ignore variants with low read coverage")
    
    args = parser.parse_args()

    cnvs = read_cnv(args.cnv_file, args.is_dpratio)

    if args.type_file!="":
        columns=args.columns.split(",")
        normal_samples, tumor_samples = read_type(args.type_file, columns)
    else:
        normal_samples, tumor_samples = set(), set()
    # print('normal samples: {}'.format(normal_samples))
    # print('tumor samples: {}'.format(tumor_samples))
    combine_snv_cnv(args.snv_file, cnvs, args.sampleID, args.min_CellFreq, args.max_PM, args.min_depth, args.output, args.is_gzipped, normal_samples, tumor_samples)
