This module is a revision to Wessim (http://sak042.github.io/Wessim/).

#### Original Wessim
Wessim is a simulator for whole exome sequencing.
Based on approaches for fragment generation, there are two modes of Wessim.
One is Wessim1 (Wessim1.py), which uses ideal target approach.
The other is Wessim2 (Wessim2.py), which uses probe hybridization approach.
Wessim2 is highly recommended when the probe sequence is available since it is more realistic.

#### Revisions
For ease of use, we revised the source code of WesSim2 and include it in our CSiTE package.  

Major changes:
1. Replace gfServer and gfClient with blat for mapping probes to the genome
2. Remove the dependencies on  Python2                                                                                                                                                      

#### Install
Several additional packages are required to run WesSim, including:
	* pysam (http://code.google.com/p/pysam/)
	* samtools (http://samtools.sourceforge.net/)
	* faToTwoBit (http://hgdownload.cse.ucsc.edu/admin/exe/ )
	* blat (http://hgdownload.cse.ucsc.edu/admin/exe/)

You may use `pip install pysam` to install pysam. 
