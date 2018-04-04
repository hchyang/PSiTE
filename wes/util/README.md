This folder contains several scripts for common use.

#### probe2fa.py
This script is used to convert probe sequences (TXT format) to FASTA format.

Sample command:
`probe2fa.py -i a.txt -o a.fa`

#### GemErr.py
This script is used to create an error model for a particular sequencing run.

Steps to create an error model:
1. Prepare the FASTA file of the reference genome:
`wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz`,
`gunzip hs37d5.fa.gz`

2. Prepare a BAM file:
`wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/NIST_NA12878_HG001_HiSeq_300x/RMNISTHS_30xdownsample.bam`
`wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/NIST_NA12878_HG001_HiSeq_300x/RMNISTHS_30xdownsample.bam.bai`

3. Create SAM file from BAM file:
`samtools view -h RMNISTHS_30xdownsample.bam 22 > RMNISTHS_30xdownsample.chr22.sam`

4. Optionally, prepare a VCF file which contains germline variants:
`wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/latest/GRCh37/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz`
Extract variants on Chr22:
` less HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz | grep -v '^#' | cut -f1,2 | grep '^22' > GV.chr22.pos`

5. Run GemErr.py:
`python GemErr.py -r 150 -f hs37d5.fa -s RMNISTHS_30xdownsample.chr22.sam -n RMNISTHS_30xdownsample_chr22 -p -e GV.chr22.pos`

#### GemStats.py
This script is used to generate statistics on error model files produced by GemErr.py.

Sample command:
`python GemStats.py -m RMNISTHS_30xdownsample_chr22_p.gzip -p -n RMNISTHS_30xdownsample_chr22_p`

#### Note
GemErr.py and GemStats.py are taken from GemSIM (https://sourceforge.net/projects/gemsim).
They are slightly revised to remove their dependencies on python 2.
