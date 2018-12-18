This folder contains several files that are required for running fa2wes.

#### File information
* S03723314_Probes.fa.gz:  
Probe sequence file obtained runing probe2fa.py on file S03723314_Probes.txt that was downloaded from https://earray.chem.agilent.com/suredesign/.
**Please decompress this file before using it.**

* S03723314_Covered_c3.bed:  
Target bed file obtained by extracting the first three columns of file S03723314_Covered.bed that was downloaded from https://earray.chem.agilent.com/suredesign/.

* RMNISTHS_30xdownsample_chr22_p.gzip:  
Error model built by GemErr.py for simulating paired end reads of length 150 bp. Refer to wes/util/README.md to see the procedure of generating this file.  

* NIST_NIST7035_H7AP8ADXX_TAAGGCGA_1_NA12878_chr1_p.gzip:  
Error model built by GemErr.py for simulating paired end reads of length 100 bp. Refer to wes/util/README.md to see the procedure of generating this file.  

#### Note
S03723314 corresponds to SureSelect Human All Exon V4.
