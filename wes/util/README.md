This folder contains several scripts for common use.

#### probe2fa.py
This script is used to convert probe sequences (TXT format) to FASTA format.

Sample command:
`probe2fa.py -i a.txt -o a.fa`

#### GemErr.py
This script is used to create an error model for a particular sequencing run.

Steps to create an error model:
1. Prepare the FASTA file of the reference genome.
```
wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
gunzip hs37d5.fa.gz
```

2. Prepare a BAM file.  
  - Downloading alignment files for Illumina WGS 2x150bp 300X data:
```
wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/NIST_NA12878_HG001_HiSeq_300x/RMNISTHS_30xdownsample.bam
wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/NIST_NA12878_HG001_HiSeq_300x/RMNISTHS_30xdownsample.bam.bai
```
  - Downloading alignment files for Illumina HiSeq Exome data:
```
wget ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/project.NIST_NIST7035_H7AP8ADXX_TAAGGCGA_1_NA12878.bwa.markDuplicates.bam
wget ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/project.NIST_NIST7035_H7AP8ADXX_TAAGGCGA_1_NA12878.bwa.markDuplicates.bai
```

3. Create SAM file from BAM file.
  - For Illumina WGS 2x150bp 300X alignment file:  
  `samtools view -h RMNISTHS_30xdownsample.bam 22 > RMNISTHS_30xdownsample.chr22.sam`

  - For Illumina HiSeq Exome alignment file:  
  `samtools view -h project.NIST_NIST7035_H7AP8ADXX_TAAGGCGA_1_NA12878.bwa.markDuplicates.bam 1 > project.NIST_NIST7035_H7AP8ADXX_TAAGGCGA_1_NA12878.bwa.markDuplicates.chr1.sam`

4. Prepare a VCF file which contains germline variants (optional).
  - Downloading reference germline variants:  
`wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/latest/GRCh37/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz`

  - Extracting variants on Chr1 and Chr22:
``` 
less HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz | grep -v '^#' | cut -f1,2 | grep '^22' > GV.chr22.pos`
less HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz | grep -v '^#' | cut -f1,2 | grep '^1' > GV.chr1.pos
```

5. Run GemErr.py.  

```
python GemErr.py -r 150 -f hs37d5.fa -s RMNISTHS_30xdownsample.chr22.sam -n RMNISTHS_30xdownsample_chr22 -p -e GV.chr22.pos
python GemErr.py -r 100 -f hs37d5.fa -s project.NIST_NIST7035_H7AP8ADXX_TAAGGCGA_1_NA12878.bwa.markDuplicates.chr1.sam -n NIST_NIST7035_H7AP8ADXX_TAAGGCGA_1_NA12878_chr1 -p -e GV.chr1.pos
```

#### GemStats.py
This script is used to generate statistics on error model files produced by GemErr.py.

Sample command:
```
python GemStats.py -m RMNISTHS_30xdownsample_chr22_p.gzip -p -n RMNISTHS_30xdownsample_chr22_p
python GemStats.py -m NIST_NIST7035_H7AP8ADXX_TAAGGCGA_1_NA12878_chr1_p.gzip -p -n NIST_NIST7035_H7AP8ADXX_TAAGGCGA_1_NA12878_chr1_p
```

#### Note
GemErr.py and GemStats.py are taken from GemSIM (https://sourceforge.net/projects/gemsim).
They are slightly revised to remove their dependencies on python 2.
GemErr.py can be used to create run-specific error models in simulation. 

As shown in the paper describing GemSIM, different error models largely affect downstream analysis, such as SNP/SNV calling. **Please check the properties of the error model before simulating WES short reads! Please also ensure that the error model is consistent with the length of short reads to simulate!**
