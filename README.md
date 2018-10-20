# Phylogeny guided Simulator for Tumor Evolution (PSiTE)

PSiTE is a phylogeny guided simulator for tumor evolution. It first simulates somatic 
variants (Single Nucleotide Variants (SNVs) and Copy Number Variants (CNVs)) 
along the history of tumor evolution and then generates Next Generation 
Sequencing (NGS) data for tumor samples. PSiTE can simulate a wide range of 
tumor data including single sector, multi-sectoring bulk tumor as well single 
cell data. It provides a powerful approach simulating tumor evolution with 
different demographic histories and allows efficient benchmarking of methods 
for clonality analysis.

## Table of Contents

<!-- toc -->

- [**1. Installation**](#1-installation)
- [**2. Usage**](#2-usage)
  + [**2.1. vcf2fa (module 1)**](#21-vcf2fa-module-1)
    * [2.1.1. Input files](#211-input-files)
    * [2.1.2. Output files](#212-output-files)
    * [2.1.3. Options](#213-options)
  + [**2.2. phylovar (module 2)**](#22-phylovar-module-2)
    * [2.2.1. Input files](#221-input-files)
    * [2.2.2. Output files](#222-output-files)
    * [2.2.3. Options](#223-options)
  + [**2.3. chain2fa (module 3)**](#23-chain2fa-module-3)
    * [2.3.1. Input files](#231-input-files)
    * [2.3.2. Output files](#232-output-files)
    * [2.3.3. Options](#233-options)
  + [**2.4. fa2wgs (module 4)**](#24-fa2wgs-module-4)
    * [2.4.1. Input files](#241-input-files)
    * [2.4.2. Output files](#242-output-files)
    * [2.4.3. Options](#243-options)
  + [**2.5. fa2wes (module 5)**](#25-fa2wes-module-5)
    * [2.5.0. Requirements](#250-requirements)
    * [2.5.1. Input files](#251-input-files)
    * [2.5.2. Output files](#252-output-files)
    * [2.5.3. Options](#253-options)
  + [**2.6. allinone (module 6)**](#26-allinone-module-6)
    * [2.6.1. Input files](#261-input-files)
    * [2.6.2. Output files](#262-output-files)
    * [2.6.3. Options](#263-options)
- [**3. Tutorial**](#3-tutorial)

<!-- tocstop -->

## 1. Installation

PSiTE is written in Python3 (>=3.5). It requires three python libraries: numpy, 
pyfaidx and PyYAML. In order to simulate whole genome sequencing (WGS) data, 
ART is also needed. We recommend using the latest version of ART (MountRainier 
or later), since older versions introduce high levels of sequencing errors. If 
users would like to simulate whole exome sequencing (WES) data, please refer to 
[section 2.5.0](#250-requirements) for the additional requirements.

PSiTE can be downloaded from github:

    git clone https://github.com/hchyang/PSiTE.git

## 2. Usage

There are six modules in PSiTE.

```
Program: psite.py (a Phylogeny guided Simulator for Tumor Evolution)
Version: 0.9.0

Usage:   psite.py <command> [options]

Command: vcf2fa     build normal genome from input germline vcf file
         phylovar   simulate somatic variations along a phylogeny
         chain2fa   build tumor genomes from somatic variants (encoded in chain files)
         fa2wgs     simulate WGS reads from normal and tumor genomes (in fasta format)
         fa2wes     simulate WES reads from normal and tumor genomes (in fasta format)
         allinone   a wrapper for NGS reads simulation by combining all individual steps
```

PSiTE starts the simulation by generating the personal (diploid) genome of an 
individual using the input VCF file (Module 1: vcf2fa). Subsequently, by taking 
the evolutionary history of the sample (For example, using coalescent 
simulations from Population Genetic modeling) and a set of user specified rate 
parameters of somatic events (e.g. SNV and CNV rates), PSiTE stimulates somatic 
variants of a sample along its evolutionary history (Module 2: phylovar). At 
the end of the second module, the simulated somatic variants will be stored in 
a special file format called chain file, which records all the somatic events 
from the root of the tree to focal clone. In the third step, PSiTE integrates 
simulated somatic variants (stored in the chain file) into the personal genome 
of the individual and builds the genomes of tumor clones (possibly down to 
individual cells) (Module 3: chain2fa). In the last step, users can use ART to 
simulate WGS reads from cancer genomes of different clones, taking into account 
their relative proportions in the tumor (Module 4: fa2wgs). Alternatively, users 
can use a similar module to simulate WES data (Module 5: fa2wes). In order to 
facilitate the simulation, PSiTE also provides a wrapper module (Module 6: 
allinone), which allows users to execute all previous steps in one command. 

### 2.1 vcf2fa (module 1)

Given a reference genome (in FASTA format) and a list of phased germline SNVs 
(in VCF format), vcf2fa will build the genomes of normal cells (in FASTA 
format). 

#### 2.1.1 Input files

##### Reference file (-r/--reference)

The reference genome (e.g. hg19.fa, in FASTA format) is specified via 
`-r/--reference`. Note that this file should contain all the chromosomes of 
interest. If a user wants to simulate a female individual, all the autosomes 
and chromosome X should be included in this file). 

##### VCF file (-v/--vcf)

A list of phased germline SNVs in the VCF format is specified via `-v/--vcf`. 
The variants in this file will be spiked into the reference genome to build the 
germline genome of an individual. Only SNVs are acceptable. PSiTE expects phased 
genotype data for the cancer individual. To be more precise, the VCF file should 
contain only one sample's genotype information and two alleles at each locus 
should be separated by '|' in the GT field of the VCF file. 

#### 2.1.2 Output files 

##### Output (-o/--output)

The directory that stores the output of vcf2fa is specified via `-o/--output`. 
After running vcf2fa, two FASTA files will be generated. They are named as 
normal.parental_0.fasta and normal.parental_1.fasta respectively (see 
`--parental` option under 2.2.3 for more details). Each file contains a 
haploid genome created by replacing the reference alleles with germline variants 
of the corresponding haplotype. 

#### 2.1.3 Options

##### -a/--autosomes and -s/--sex_chr

These two options allow users to specify the set of chromosomes to simulate. In 
the `--autosomes` option, the chromosomes should be separated by comma 
(e.g. '1,2,3,4,5'). To avoid long names, users can also use two dots to 
represent a continuous block of chromosomes. For example, 'chr1..4,chr7' is 
equivalent to 'chr1,chr2,chr3,chr4,chr7'. In vcf2fa, it is required to specify 
the values for `--autosomes`, following which vcf2fa generates two parental 
copies for all the specified autosomes. The parameter `--sex_chr` is optional. 
`--sex_chr X,Y` and `--sex_chr X,X` will specify a male and a female individual 
respectively. (Note that users should be careful about the chromosome names. 
Chromosome names not found in the reference file will lead to the early 
termination of the program.) 

### 2.2 phylovar (module 2)

phylovar is the core module of PSiTE. It can jointly simulate SNVs and CNVs of 
a tumor sample along the history of a cell lineage tree (i.e. phylogenetic 
tree). Phylovar can be used either in a standalone mode or in an integrated 
mode to simulate NGS data of a tumor sample together with other modules in PSiTE.  

#### 2.2.1 Input files 

##### Tree file (-t/--tree)

To run the module, a cell lineage tree (i.e. a phylogenetic tree in the Newick 
format) which captures the ancestral relationship of the tumor cells is 
required. The [ms](http://home.uchicago.edu/rhudson1/source/mksamples.html) 
program is recommended to generate the tree (with `-T` option in ms), and it has 
the full-assemblage of options needed to generate complex demographic histories 
of the sample from Population Genetics.

A typical tree is as follows:

    ((2:0.083,4:0.083):0.345,(5:0.322,(1:0.030,3:0.030):0.292):0.105);

##### Trunk variants file (--trunk_vars)

This file contains known trunk variants, specified by `--trunk_vars`. The 
truncal variants can also be simulated randomly using the `--trunk_length` 
option (see Section 2.2.3).

The format of this file is shown below:

    #chr hap start end var target
    1 0 364645 364646 0
    1 1 464646 466646 +3 194327,149250611,5687878
    1 1 464650 464651 1
    1 1 464660 464661 1 0
    1 1 464670 464671 1 2,3
    2 0 465636 465637 1 
    2 1 468637 472633 -1

- **chr**:      The chromosome on which the variant locates.
- **hap**:      The haplotype copy of the chromosome on which the somatic 
variant resides on. The specified haplotype should match what is given in the 
option `--parental`  (see Section 2.2.3). For example, if '0011' is specified 
via `--parental` for a given chromosome, the hap column allows users to specify 
the truncal events to be on one of the four haplotypes (0/1/2/3).  
- **start**:    The starting position of the variant (0-based, inclusive).
- **end**:      The end position of the variant (0-based, exclusive).
- **var**:      The type of the variant. 0/1/2: SNV, -1: deletion, 
+N: amplification. 0/1/2 represent different types of base substitutions (see 
section `--tstv` option under section 2.2.3).
- **target**:   This column is optional and is only relevant when the focal 
variant is a SNV or an amplification.  **1)** For an amplification with N new 
copies (i.e. +N in the var column), the `target` column can be a list with N 
integers (seperated by commas) or omitted. The integer list indicates all the 
insert loci of the new copies in this column. Without this information, the 
amplification will be treated as a tandem amplification. **2)** For a SNV, 
`target` indicates the copy of the segment that carries this variant. If a SNV 
is covered by an amplification with N new copies, the value in this column can 
be either an integer in the range of [0,N], or a list of numbers separated by 
commas representing multiple copies carrying this mutation. If there is no 
value in this column, the focal SNV is assumed to occur early (i.e. if there 
are overlapping CNV events, SNVs are assumed to be earlier than the CNV event). 
In this case, all the new copies of the amplification will carry the SNV. In 
addition to this default setting, users can also specify alternative scenarios 
using this column. For example:

```
#chr hap start  end    var target
1    1   464646 466646 +3 194327,149250611,5687878
1    1   464650 464651 1
1    1   464660 464661 1 0
1    1   464670 464671 1 2,3
```

In this example, there is an amplification in the region 1:464646-466646 
covering all the three SNVs. The first copy of this amplification is 
inserted before 194327. The second one is inserted before 149250611. The 
third copy is inserted before 5687878. And the three SNV records correspond 
to the following three scenarios: 

1. The original copy and all the amplified copies will carry the first SNV 
(1:464650-464651). 
2. Only the original copy (before amplification) carries the second SNV 
(1:464660-464661).
3. Only the second and the third new copies carry the third SNV 
(1:464670-464671).

##### Configuration file (--config)

The configuration file, specified by `--config`, enables users to specify 
mutation parameters for the whole genome or individual chromosome. A typical 
configuration file contains two sections: the genome section and chromosome 
section. In the genome section, users can specify the default parameters shared 
across all chromosomes. In the chromosome section, users can customize 
parameters specific to an individual chromosome. The local/chromosome 
parameters will overwrite the global/genome parameters. Separate parameter 
sections provide a flexible way to simulate somatic variants across the genome 
with different mutation rates.

The configuration file should be specified in the [YAML](http://yaml.org/) 
format. Here is an example of the configuration file (for parameters listed in 
the YAML file, please see section 2.3.3 for further details).

    genome:
        snv_rate: 2.0
        cnv_rate: 0.6
        del_prob: 0.5
        tandem_prob: 1.0
        cnv_length_beta: 200
        cnv_length_max: 400
        copy_parameter: 0.5
        copy_max: 5
        parental: '01'
        tstv: 2.0
        length: 109000
    chromosomes:
        - '1':
            length: 100000
            parental: '00'
        - '2':
            snv_rate: 2.0
            length: 9000

In a configuration file, all the parameters under the genome section must be 
specified. The parameters specified under snv_rate, cnv_rate as well as length 
in the genome section must be the sum of corresponding values across all the 
chromosomes. The parental value and the chromosome names should be quoted with 
quotation marks.

##### Affiliation file (--affiliation)

By default, phylovar simulates the genomic profile of a single sample. Using 
this option, phylovar can simulate multi-sector tumor samples. In the 
affiliation file, users should designate the sample affiliation of the tumor 
cells to different sectors. An example affiliation file is shown below:

    #sector purity depth prune_p cells
    sector1 0.6    -     0.05    1,2,3,4..5000
    sector2 0.6    100   0.05    5001..10000

- **sector**: The id of each sector. 
- **purity**: The purity of each sector (the proportion of cells that is tumor 
cells in each sector) .
- **depth**: The depth of each sector. It's different from the 'depth' in 
sector file (see section 2.4.1 and section 2.5.1). The 'depth' here is used
to simulate read count for the simulated SNVs (see 'SNV files' under section 
2.2.2). It dose not affect the coverage of sequencing reads in module 
fa2wgs/fa2wes. User can use '-' here to disable the read count simulation.
- **prune_p**: The pruning proportion for each sector. Subtrees that have 
descendants fewer than the specified proportion will be trimmed away (see 2.2.3 
--prune for details).
- **cells**: The list of cells in each sector. It can be specified as a 
list of cell names separated by commas or a continuous block of cells separated 
by two dots (e.g. 'cell1..4,cell7' is equivalent to 
'cell1,cell2,cell3,cell4,cell7').

##### CNV length distribution file (--cnvl_dist)

By default, phylovar simulates CNVs whose length distribution is an exponential
distribution (check --cnv_length_beta in section 2.2.3). With this option, users 
can specify a file, which contains an empirical distribution of CNV length.
An example input file is shown below:

    #low high prob
    10000 20000 0.2
    20000 80000 0.6
    80000 100000 0.2

- **low**: The lower bound of the bin (inclusive)
- **high**: The upper bound of the bin (exclusive)
- **prob**: The probability of a simulated CNV whose length falls in the range of
[low,high).

Each record is specifying a bin. In each bin the length of CNVs follows an
uniform distribution. Sum of the probability of all bins should be 1. This
CNV length setting will override other settings of CNVs' length (i.e. 
cnv_length_beta and cnv_length_max settings in command line or in the 
configuration file).

#### 2.2.2 Output files

Module phylovar can output multiple files to facilitate benchmarking of methods 
used in cancer genomics. Since the evolutionary history of the sample can be 
very complex, outputting all the files can be quite resource intensive. Phylovar 
allows user to toggle different options for outputting the mutational 
information. Output files that are not foremost will be labelled as optional. 

##### SNV files (-S/--snv)

This option specifies the directory to store the SNV files of all sectors. 
The SNV file contains the frequency information of simulated SNVs in each 
sector. In addition to the SNV file for each sector, there is also a SNV file 
for the whole tumor sample named 'tumor.snv' under this folder. There are at
least five columns in the SNV file.

- **chr**: The chromosome on which the SNV locates.
- **start**: The start position of the SNV (0-based, inclusive).
- **end**: The end position of the SNV (0-based, exclusive).
- **form**: The type of base substitution (0 is transition, 1 and 2 are two 
types of transversions, see the details for option --tstv under section 2.2.3).
- **frequency**: The frequency of the alternative allele in tumor sample 
(taking into account tumor purity).

There are two more columns in SNV files if users simulate read count
for each SNV by specifiying 'depth' in affiliation file or `--depth`.

- **rcount**: The read count of SNVs. The read count of each sector is 
simulated with the depth specified in affiliation file, and the read count of 
the SNVs in whole tumor is simulated with the depth specified by `--depth` in 
section 2.2.3. The vales in this column is in the foramt of 'N:M'. N is the 
read count of alternative allele. M is the read count of the total reads 
covering this loci. They are simulated by taking purity and CNVs into account.
- **rfreq**: The frequency of the alternative allele in sectors or the whole 
tumor sample. It's calculated by N/M. 


##### CNV files (-V/--cnv) 

This option specifies the directory to store the CNV files across all sectors. 
It contains the information of all simulated CNVs. In addition to the CNV file 
for each sector, there is also a CNV file for the whole tumor sample named 
'tumor.cnv' under this folder. There are five columns in the CNV files.

- **chr**: The chromosome on which the CNV locates.
- **start**: The start position of the CNV (0-based, inclusive).
- **end**: The end position of the CNV (0-based, exclusive).
- **copy**: The copy changes of the CNV (e.g. -1 stands for deletion, +N stands
for an amplification with N new copies).    
- **carrier**: The number of tumor cells within the sample carrying the CNV

##### SNV genotype file (--snv_genotype) (optional)

The output file specified by `--snv_genotype` contains the genotype information 
of each tumor cell. Each SNV has one record in this file. The first three 
columns are the coordinates of the SNV. The fourth column is the form of the SNV 
event (transition or transversion, see section 2.2.3). Subsequently, there is 
one column (genotype) per tumor cell. The SNV genotype is in the form of 'M:N', 
in which M denotes the number of alternative alleles and N denotes the number of 
reference alleles. The columns are: 

- **chr**: The chromosome of the SNV.
- **start**: The start position of the SNV (0-based, inclusive).
- **end**: The end position of the SNV (0-based, exclusive).
- **form**: The mutational type of the SNV (0 is transition, 1 and 2 are two 
types of transversions, see section 2.2.3).
- **cell1**: The genotype of cell 1.
- **cell2**: The genotype of cell 2.
- **etc ...**

Users should be careful whether to toggle this option since the output of this 
option can be very large. 

##### CNV genotype file (--ind_cnvs) (optional)

The output file specified by `--ind_cnvs` contains the CNVs on each parental 
copy of each cell in the sample.

    #cell parental chr start     end       copy
    1       0      1   7912422   7930111   +2
    1       1      1   43110140  43341629  +1
    2       0      1   2255734   2299608   -1
    2       0      2   22660687  22788472  -1
    2       1      2   59756841  61142076  +3

There are six columns in this file:

- **cell**: The id of the cell in the sample.
- **parental**: The parental copy in which the variant locates (0 means one of 
the parental copy, and 1 means the other copy).
- **chr**: The chromosome on which the CNV locates.
- **start**: The start position of the CNV (0-based, inclusive).
- **end**: The end position of the CNV (0-based, exclusive).
- **copy**: The copy number of the CNV. -1: deletion; +N amplification.
Users should be careful whether to toggle this option since the output of this 
option can be very large. 

##### Total CNV profile file (--cnv_profile) (optional)

The output file specified by `--cnv_profile` contains the total CNV profile 
across the whole tumor sample. It outputs a step function of copy numbers across 
the genome. There are four columns in this file:

- **chr**: The chromosome of the segment
- **start**: The start position of the region (0-based, inclusive).
- **end**: The end position of the region (0-based, exclusive).
- **local_cp**: The copy number of the local region. For example, if there are 
1000 cells in the tumor sample, the local copy number for most regions is 2000 
as cells are diploid. For a genome region whose parental 0 copy are lost in 200 
cells, the local copy number will be 1800. 

##### Variant tree file (--nhx/--NHX) (optional)

The variant tree file (in NHX format) outputs locations of the somatic variants 
on the phylogenetic tree.  The NHX format contains IDs of all the nodes in the 
phylogenetic tree and somatic variants along all the branches. This allows users 
to reconstruct the entire history of somatic events in tumor evolution. If users 
specify the option with `--nhx`, phylovar will output the pruned version of the 
tree (see `--prune` option under 2.2.3 for details) together with all somatic 
variants. If users specify the option with `--NHX`, the output contains the 
original tree (before pruning) together with all the somatic variants. In the 
file, the somatic variants are represented in the format of 'chr#start#end#var'. 
Here, entry 'var' has the same format as column 'var' in the trunk variant file 
(see section 2.2.1).

##### Node variant file (--nodes_vars) (optional)

The node variant file, specified by `--nodes_vars`, contains the somatic 
variants (SNVs/CNVs) occurring on the branch leading to each node in the 
phylogenetic tree. This is the collapsed version of the variant tree file where 
the phylogenetic information is simply replaced with the node information. 
There are six columns in this file:

- **node**: The ID of the focal node. Each node has the format 'nodeX', in which 
X is an integer starting from 1.
- **chr**: The chromosome on which the variant locates.
- **hap**: The haplotype copy of the chromosome on which the somatic variant 
locates. (see `--trunk_vars` option under section 2.2.1) 
- **start**: The start position of the variant (0-based, inclusive).
- **end**: The end position of the variant (0-based, exclusive).
- **var**: The type of the variant. 0/1/2: SNV, -1: deletion, +int: amplification 
(see `--trunk_vars` option under section 2.2.1).

##### Node CCF file (--nodes_ccf) (optional)

The node CCF file, specified by `--nodes_ccf`, contains the Cancer Cell 
Fraction (CCF) information of each node (all the nodes after pruning) in each 
sector. For example, for a sector containing 1 million cells, if its purity is 
0.6, the tumor cells in the whole sector is 600,000. Then, if an inner node, 
nodeX, in the tree has 150,000 descendants (leaf nodes) belonging to this 
sector. The CCF of nodeX in this sector is 0.25.

There are number_of_sectors+2 columns in this file:

- **node**: The ID of the focal node. Each node has the format 'nodeX', in which 
X is an integer starting from 1.
- **sector1**: The CCF of the focal node in sector 1.
- **sector2**: The CCF of the focal node in sector 2.
- **...**
- **sectorX**: The CCF of the focal node in sector X.
- **tumor**: The CCF of the focal node in the whole tumor sample.

##### Chain file (--chain) (optional)

The chain files store somatic events of the focal clone. If specified with the 
`--chain` option, phylovar will output the chain files to a folder. The chain 
files are required by module chain2fa to convert the normal genome into tumor 
genomes. The chain file contains mutational events from the root of the tree to 
each individual tip node (corresponding to tumor clones or single cells). There 
is one chain file per tip node. Below is an example of the chain file:

    >1_Hap0 parental:0
    1       0       33      REF
    1       55      99      AMP     +2/2
    1       33      35      REF
    1       35      55      DEL     -1
    1       55      99      AMP     +1/2
    1       55      99      REF
    1       55      100000  REF
    >1_Hap1 parental:0
    1       0       33      REF
    1       33      34      SNV     1
    1       34      100000  REF

In the chain file, there is one block for each haplotype of a chromosome. Within 
each block, there is a header line and a body section. Lines which start with 
'>' are header lines. All other lines are components of the body section.

There are two fields in a header line. The first is the haplotype name and the 
second is the parental copy of that haplotype. Each haplotype is in the format 
of '{chrName}\_Hap{index}'. For example, if a user simulates three copies of 
sequence 1 with parental 001, there will be another block with header 
'>1_Hap2 parental:1' in addition to the two blocks in this example.

There are five columns in the body section:

- **chr**:  The chromosome on which the segment locates.
- **start**: The start position of the segment (0-based, inclusive).
- **end**: The end position of the segment (0-based, exclusive).
- **type**: The type of the segment (SNV: single nucleotide variation, AMP: 
amplification, DEL: deletion, REF: reference). The amplified regions are 
represented by multiple overlapping records. For example, line 3,6 and 7 in the 
example chain file indicate that there is an amplification in region 1:55-99.
- **var**: The type of the variant. 0/1/2: SNV, -1: deletion. For 
amplifications, this column is in the format of `+N/M`. N is the copy index and 
M is the total number of new copies of the amplification events. In the example
chain file, the first copy of the amplification (1:55-99) is inserted before 
1:55, and the second one is inserted before 1:33. 

##### Tipnode map file (--map) (optional)

The tipnode map file is stored in the folder specified by `--map`. It contains 
the descendant information of each tip node (corresponding to a tumor clone). 
These files will be used later by module fa2wgs/fa2wes to compute the sequencing 
coverage for each tumor genome. If multi-sector samples are simulated, there 
will be one tipnode map file for each sector. They will be named as 
'{sectorname}.tipnode.map'. By default (single sector), there will be only one 
map file in the folder, named 'tumor.tipnode.map'. 

There are three columns in a node map file:

- **tip_node**: The id of the tip node in the phylogenetic tree. If pruning 
option is turned on, the tip node refers to the leaves after the pruning step.
- **cell_count**: The number of cells (i.e. number of leaves) descending from 
this node.
- **cells**: Names of the cells (i.e. names of leaves) descending from this node.

##### Log file (-g/--log) 

The log file, specified by `-g/--log`, contains the log information including 
commands to call phylovar and the random seed. Users can use the information in 
this file to repeat the simulation.

#### 2.2.3 Options

##### --name

This option specifies the name of the simulated sequence (e.g. chromosome name).

##### --length

This option specifies the length of the simulated sequence. 

##### --snv_rate and --cnv_rate
These two options set two most important parameters in the simulation: the 
mutation rates of SNVs (specified by `--snv_rate`) and CNVs (specified by 
`--cnv_rate`) along the history of tumor evolution. Note that these rates 
represent the total mutation rate of the target segment. The mutational events 
are then simulated according to a Poisson process with user-specified rates 
(see Notes for extra discussions).

##### --tstv

This option allows users to specify the rate of transitions to transversions.  
phylovar only simulates somatic events (before referring to the reference 
genome). It uses 0/1/2 to represent different types of substitutions. 0 stands 
for transition. 1 and 2 stand for two types of transversions. After linking 
somatic changes to the reference genome, these SNVs will later be translated to 
actual nucleotides in the module chain2fa. The rules of translating 0/1/2 to 
nucleotide changes are summarized in the following table. 

reference | 0 | 1 | 2
----------|---|---|---
 N | N | N | N
 A | G | C | T
 G | A | C | T
 C | T | A | G
 T | C | A | G

##### --cnv_length_beta

This option specifies the mean length (parameter beta or mean of the exponential 
distribution) of the simulated CNV events. Note that phylovar assumes the length 
of CNVs follows an exponential distribution.

##### --cnv_length_max

This option sets the upper limit of CNVs' length. Setting an upper bound will 
effectively truncate the exponential distribution at this limit. The probability 
distribution will be renormalized taking into account the truncated probability.

##### --del_prob 

A CNV event can be either a deletion or an amplification. Through this option, 
users can specify the probability that a CNV event is a deletion.

##### --copy_max 

This option sets the upper bound of the copy number of an amplification event. 

##### --copy_parameter

When an amplification event happens, phylovar randomly picks a copy number 
according to a geometric-like distribution with Pr(n+1)=p\*Pr(n). The parameter 
p is specified by '-c/--copy_parameter'. The overall distribution will be 
normalized so that the total probability is  1.

##### --parental 

Since aneuploidy at the chromosomal level is widespread in tumor cells, phylovar 
allows users to specify the local ploidy of a cancer genome via --parental. The 
default value of this option is '01', with which phylovar simulates two copies 
of the chromosome (one copy from each of the normal haploid genomes). '001' 
specifies two copies from the first haploid genome ('0') and one copy from the 
second haploid genome ('1'). Note that aneuploidies are simulated as truncal 
events shared among all tumor cells.

All of the above parameters of phylovar can be specified in the configuration 
file (`--config`). This is quite convenient if you want to simulate somatic 
variants in a genome, which contains multiple chromosomes. The format of the 
config file is described in section 2.2.1 Input files (Configuration file).

##### --trunk_length/--trunk_vars 

The phylogenetic tree only describes the polymorphic part of the sample. Users can 
also simulate the truncal part of tumor evolution (tumorigenesis) through two 
different approaches: a) specifying the trunk length using `--trunk_length`.  
b) specify the events on the trunk in a file through `--trunk_vars` (see section 
2.2.2). Option `--trunk_length` accepts a single float number that represents 
the length of the branch from the root of the phylogenetic tree to the most common 
ancestor of the sample. 

##### --prune

When the number of simulated cells is very large, computational load can be very 
heavy. Given the fact that mutational events on terminal branches are extremely 
rare, a pruning algorithm is implemented to trim the tree. The parameters of the 
algorithm are specified via `--prune`. Given a tree for a population of 
1,000,000 cells, all subtrees with <10,000 tips will be trimmed away after 
setting `--prune 0.01`. This means that there will be no polymorphic variants on 
subtrees with <10,000 tips. Namely, terminal cells belonging to the same subtree 
with <10,000 tips will have the same genome. 

##### -s/--sex_chr

This is the same option as option `--sex_chr` in module vcf2fa. Setting this 
parameter will affect the ploidy of the sex chromosomes.

##### --purity

This option specifies the proportion of tumor cells in the whole sample. With 
this option, phylovar can output true frequencies of each variant in the whole 
sample by taking into account tumor purity (section 2.2.2 output file: SNVs 
file). 

##### --depth

This option specifies the depth of the whole tumor sample for read count 
simulation. With this option, phylovar can output read count of each variant 
in the whole sample by taking into account tumor purity and CNVs effect
(section 2.2.2 output file: SNVs file). 

##### --random_seed

This option sets the seed for the random number generator. This random seed 
should be an integer between 0 and 2\*\*31-1. 

##### --loglevel [DEBUG, INFO]

This option specified the verbosity level of the log file. If the level is 
DEBUG, the steps in the simulation will be written to the log file with higher 
verbosity.

#### Notes

By default, the branch length of trees generated from the ms program is measured 
in 4N generations. If we set `-r` to 100 in phylovar (PSiTE), this is equivalent 
to setting the population rescaled parameter -t to 100 in ms (see ms manual for 
details).

### 2.3 chain2fa (module 3)

After generating genomes of normal cells with vcf2fa and chain files of tumor 
cells with phylovar, chain2fa then builds the genome sequences for each tumor 
cell (or clone). 

#### 2.3.1 Input files

##### Normal fasta (-n/--normal)

The genomes of an individual (with two parental haplotypes) are specified via 
`-n/--normal`. These genomes are generated by module vcf2fa using the germline 
variants. In other words, they serve as  the genetic background of the 
individual. Two haploid genomes of the normal cell should be specified as a list 
separated by comma, such as, `--normal normal.parental_0.fa,normal.parental_1.fa`.  
Here, the order of the fasta files is important. The first fasta will always be 
treated as parental 0, and the second one will be treated as parental 1. In the 
example given above, normal.parental_0.fa will be used to build the tumor 
haplotype corresponding to parental haplotype 0 and normal.parental_1.fa will be 
used to build the tumor haplotype from parental haplotype 1.

##### Chain file (-c/--chain)

The directory which contains all the chain files is specified by option 
`-c/--chain`. All the chain files named as 'node\*.chain' within this folder 
will be used to build genomes of tumor clones or single cells.
 
#### 2.3.2 Output files 

##### Output (-o/--output)

The folder storing simulated tumor genomes is specified via option 
`-o/--output`. By default, there are two FASTA files (corresponding to parental 
0 and parental 1 respectively) per chain file. They are named as 
'node\*.parental\*.fa'. 

#### 2.3.3 Options

##### --width

This option specifies the per line width of the genome sequence. Setting this 
to be 100 will output the genome sequences in chunks of 100 bp. 

##### --cores

This option specifies the number of cores used to run this module.

### 2.4 fa2wgs (module 4)

After running the first three modules (vcf2fa, phylovar and chain2fa), PSiTE 
have generated the genomes of normal/tumor clones (cells). Module fa2wgs then calls 
[ART](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm) 
to simulate whole-genome sequencing data from these genomes. fa2wgs first 
computes the sequencing coverage of each genome according to its proportion in 
the tumor sample and then calls ART to simulate short reads from each genome. 

#### 2.4.1 Input files

##### Normal fasta (-n/--normal) and tumor fasta (-t/--tumor)

The folder containing the FASTA files of normal genomes is specified by 
`-n/--normal`. The corresponding tumor folder is specified by `-t/--tumor`. 
fa2wgs expects two fasta files (normal.parental_0.fa/normal.parental_1.fa) under 
the normal fasta folder and two fasta files 
(node\*.parental_0.fa/node\*.parental_1.fa) for each tip node under the tumor 
fasta folder.

##### Tipnode map files (-m/--map)

The folder containing the tip node map files is specified via `-m/--map`. Each 
file contains the compositional information of a tip node, including the number 
of cells in the sample that have the same tumor genome represented by the node. 
Please see the 'Tipnode map files' under 2.2.2 for more details.

##### Sector file (--sectors)

The sector file is specified by `--sectors`. In this file, users can set 
different purity and sequencing depth for different tumor sectors. WGS data will 
be simulated only for sectors listed in this file. An example of the sectors 
file is shown below:

    #sector purity depth
    sector1 0.6 100
    sector2 0.6 110

There are three columns in a sectors file:

- **sector**: The ID of the focal sector. 
- **purity**: The purity of sectors (the cancer cell proportion in a sector).
- **depth**:  The depth of WGS data for sectors to simulate.

When the sector file is supplied, fa2wgs will ignore options `-d/--tumor_depth` 
and `-p/--purity`. Without this option, fa2wgs will simulate WGS data for all 
the sectors which has a map file in the folder specified by `-m/--map`. In this 
case, the purity and sequencing depth for each sector are taken from 
`-p/--purity` and `-d/--tumor_depth`. In addition to the sectors specified in 
the affiliation file, users can also specify 'tumor' in the sectors file to make 
fa2wgs to simulate the WGS data for the whole tumor sample.

#### 2.4.2 Output files

##### Output folder (-o/--output)

The folder containing the short reads generated by ART is specified by 
`-o/--output`. The FASTQ files of tumor sample are stored in a subfolder called 
'tumor'. The FASTQ files of normal sample are stored in a subfolder called 
'normal'. 

##### Log file (-g/--log)

For each simulation, fa2wgs calls ART multiple times to simulate NGS reads from 
individual tumor clone (i.e. cancer genome). This file, specified by `-g/--log`, 
stores command settings of each command and can be used to replicate the 
simulation process.

#### 2.4.3 Options

##### -d/--tumor_depth

This option specifies the mean depth of the tumor sample. 

##### -D/--normal_depth

For most of the cancer genomic analysis, paired normal sample is also needed. 
This option specifies the mean depth of normal sample.

##### -p/--purity 

This option specifies the purity of tumor samples. Purity refers to the 
percentage of cells that are tumor cells. 

##### --rlen

This option specifies the length of the simulated short reads.

##### --art

This option specifies the command to call ART. fa2wgs provides a very flexible 
way to call ART. Users can pass all parameters to ART by `--art` except 
parameters handled by fa2wgs itself (rlen, fcov, in, out, id, and rndSeed). For 
example, users can use `--art '/path/to/ART/art_illumina --noALN --quiet 
--paired --mflen 500 --sdev 20'` if ART is not installed system-wide. Users can 
also use `--art 'echo art_illumina --noALN --quiet --paired--mflen 500 --sdev 
20'` to just print out the command to call ART.

##### --cores

This option specifies the number of cores used to run fa2wgs. With this option, 
users can use multiple CPUs to reduce simulation time.

##### --separate

By default, fa2wgs merges the simulated NGS data of all tumor genomes into one 
file. This option allows fa2wgs to store the individual fastq files separately. 

##### --single
The option indicates NGS simulation will be in single cell mode. After 
specifying this option, the value of `--tumor_depth` is the depth of each tumor 
cell (not the total depth of tumor sample anymore).

##### -s/--random_seed

This option sets the seed for random number generator. This number is passed to 
ART (i.e. `--rndSeed` option in ART) to initiate the simulation of short reads.

### 2.5 fa2wes (module 5)

In addition to simulating WGS data, PSiTE can also simulate WES data (Illumina) 
by running module fa2wes. Similar to module fa2wgs, fa2wes first computes 
sequencing coverage (or read number) for each genome in the tumor/normal sample 
and then calls a WES simulator to simulate short reads. 

To avoid reinventing the wheel, we integrate two popular WES simulators: 
[WesSim](http://sak042.github.io/Wessim) and 
[CapSim](https://github.com/Devika1/capsim). Both of them simulate shorts reads 
by combining all the three stages of capture sequencing: fragmentation, fragment 
capture and in silico sequencing. WesSim uses an existing NGS simulator 
[GemSim](https://sourceforge.net/projects/gemsim) to generate short reads. 
CapSim simulates reads by copying from captured fragments with random errors 
introduced and fixed quality scores. The coverage distribution of reads 
generated by CapSim has been found to resemble real exome capture data. To get 
realistic quality values for the short reads, we created a hybrid simulator 
CapGem by combining the first two steps of CapSim (fragmentation and fragment 
capture) with the last step of WesSim (sequencing by GemSim). Specifically, the 
probe-based version of WesSim (WesSim2) and CapGem are supported in fa2wes. 

Most options for module fa2wes are similar to those for the module fa2wgs. 

#### 2.5.0 Requirements

Due to the complexity of simulating exome data, this module requires the 
installation of the following software.

##### Snakemake
Unlike WGS simulation, there are often several steps for simulating WES data. 
To create a smooth workflow, we use 
[snakemake](http://snakemake.readthedocs.io/en/latest/) (version >= 3.8), 
a scalable workflow management system.  

There are several ways to install snakemake. Please follow the instructions 
[here](http://snakemake.readthedocs.io/en/latest/getting_started/installation.html).

##### WesSim

For ease of use, we revised the source code of WesSim2 (probe-based version) and 
include it in our PSiTE package.  

Several additional packages are required to run WesSim:

- [pysam](http://code.google.com/p/pysam/)
- [samtools](http://samtools.sourceforge.net/)
- [faToTwoBit](http://hgdownload.cse.ucsc.edu/admin/exe/)
- [blat](http://hgdownload.cse.ucsc.edu/admin/exe/)

##### CapGem

CapGem requires a Java Runtime Environment (Java Runtime Environment >=1.8). 
Several additional packages are required to use CapGem for simulation:

- [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- [samtools](http://samtools.sourceforge.net/) (version >= 1.5)

To install CapGem, one can navigate to the folder containing the source code of 
CapGem (PSiTE/wes/capgem) and then run `make`. 

A sample command is `make install INSTALL_DIR=./ MXMEM=8000m SERVER=true`. If 
`INSTALL_DIR` is specified to be a directory other than './', this directory 
should be added to environment variable PATH. After installation, user can also 
adjust the maximum memory used by Java by adding `export _JAVA_OPTIONS="-Xmx12g"` 
to the file .bashrc under home directory.

Due to the dependence on samtools to build index for tumor genomes, CapGem can 
only be used for genomes with longest chromosome being less than 512 M.

##### Computing resources 

The module fa2wes is computationally intensive, since it includes several steps: 
building index, mapping probes, and simulating reads. The process can be 
accelerated by using multiple cores in a single computer or using multiple 
cluster nodes. 

The computational resources used by different simulators are quite different. 
Given a haploid genome and a computer with 8 GB memory, WesSim (with default 
parameters in the configure file) takes about 10 minutes (1 GB memory, 700 MB 
disk space) to build an index and about 40 minutes (4 GB memory, 70 MB disk 
space) to map probe sequences onto the genome, whereas CapGem (with default 
parameters in the configure file) takes about 5 hours (7.8 GB memory, 4 GB disk 
space) to build an index and about 30 minutes (7 GB memory, 900 MB disk space) 
to map probe sequences. The time used for simulating reads will dependent on the 
number of reads to simulate. In our experiences, WesSim is often much faster 
than CapGem in generating short reads. 

In principle, one can specify the read depth/number for tumor and normal samples 
in a single command. To save time, one can simulate tumor and normal samples in 
parallel by running two fa2wes commands at the same time, with one command 
specifying the read number of tumor (normal) sample to be zero while simulating 
normal (tumor) sample.

#### 2.5.1 Input files 

Similar to module fa2wgs, this module requires the FASTA files of normal and 
tumor genomes as well as tumor chain files. Additionally, a probe file and a 
target file are required.

##### FASTA files of normal genomes (-n/--normal)

These fasta files for the normal genomes are specified via `-n/--normal`. 
There should be two FASTA files (normal.parental_0.fa/normal.parental_1.fa) 
under this folder.

##### FASTA files of tumor genomes (-t/--tumor)

These fasta files for the tumor genomes are specified via `-t/--tumor`. 
There should be two FASTA files (node\*.parental_0.fa/node\*.parental_1.fa) 
for each tumor genome under this folder.

##### Tipnode map files (-m/--map)

These files are specified via `-m/--map`. Each file contains the compositional 
information of a tumor genome (i.e. which cell belongs to this tumor 
clone/genome). Please see the 'Tipnode map file' under 2.2.3 for more details of 
this file.

##### Probe file (--probe)

The probe file (FASTA format) is specified via `--probe`. This file contains the 
probe sequences for target capture. Probe files for different sequencing 
platforms can be obtained from the respective vendor's website. For example, the 
probe sequences for Agilent's SureSelect platforms can be downloaded from 
[here](https://earray.chem.agilent.com/suredesign/) (One has to register and log 
in the website for downloading). The download option can be found in the menu by 
first clicking 'Find Designs' and then 'SureSelect DNA', followed by clicking 
'Agilent Catalog Designs'. To obtain a FASTA file from the downloaded TXT file, 
user can use the script probe2fa.py provided in our PSiTE package (under 
directory wes/util/). We provide an example probe file 
wes/example/S03723314_Probes.fa.gz for users (need extration before using). 
Check the [wes/example/README.md](wes/example/README.md) file for details.

##### Target file (--target)

The target file (BED format) is specified via `--target`. This file contains the 
positions of target regions along the genome. Similar to the probe file, Target 
files for different sequencing platforms can also be obtained from vendor's 
website. Users can find an example target file file wes/example folder of PSiTE
package (wes/example/S03723314_Covered_c3.bed). 
Check the [wes/example/README.md](wes/example/README.md) file for details.

##### Error model file (--error_model)

This option specifies the file containing the empirical error model for NGS 
sequencing generated by GemErr.py (please see wes/util/README.md for details). 
GemErr.py was adapted from the GemSim package and it can tabulate an error model 
from real sequencing data. With the error model learned from the real data, 
users can simulate more realistic short reads. Users can use a default error 
model provided under wes/example folder of PSiTE package 
(wes/example/RMNISTHS_30xdownsample_chr22_p.gzip).
Check the [wes/example/README.md](wes/example/README.md) file for details.

##### Sector file (--sectors)

The sector file is specified via `--sectors`. This file contains purity and 
depth of WES data for different tumor sectors. Please see the 'Sector file 
(--sectors)' under 2.4.1 for more details of this file. 

#### 2.5.2 Output files

##### Output of WES simulator (-o/--output)

A WES simulator typically generates many files corresponding to different tumor 
genomes. These files are under the folder specified via `-o/--output`. There can 
be multiple subfolders generated under this folder. If the simulation is only 
done for the tumor or normal sample, there will be an additional folder 'tumor' 
or 'normal' which contains similar hierarchies of subfolders. 

In each subfolder, there are several subfolders that contain the intermediate 
and final (underlined) results in simulation:

- **config**: This folder contains the configuration files of snakemake for 
running different simulators. Users can change the parameters of tools used in 
simulation directly in the respective configuration file for specific needs.
- **genome_index**: This folder contains the index of reference genomes. 
- **mapping**: This folder contains the results of mapping the probe sequence to 
the reference genomes. 
- **wessim_reads/capgem_reads**: This folder contains simulated short reads for 
different genomes in the sample. 
- **frags**: This folder contains the sequencing fragments captured by CapGem.
- **merged**: This folder contains the simulated short reads for the 
tumor/normal sample. For the tumor sample, multiple tumor genomes (i.e. clones) 
are merged to represent the tumor sample. 
- **separate**: This folder contains simulated short reads for each genome (i.e. 
tumor clone) in the sample. It is generated when using option `--separate` 
or `--single` (see section 2.5.3). When multi-sector data are simulated, there 
will be a subfolder for each sector under this folder.

All the other subfolders contain intermediate results and can be deleted, 
including 'log' which store files indicating the completeness of snakemake 
tasks, 'stdout' which stores the standard output and standard errors of 
snakemake commands, and '.snakemake' which stores the files generated by 
snakemake.

To simulate short reads with different parameters without changing the 
underlying genomes, one can delete all the subfolders except 'config', 
'genome_index' and 'mapping', before rerunning fa2wes.

##### Log file (-g/--log)

The log file is specified via `-g/--log`. It stores the commands for calling the 
WES simulator, which can be used to replicate the simulation process.

#### 2.5.3 Options 
 
##### -p/--purity 

This option specifies the purity of the tumor sample.

##### --ontarget_ratio 

This option specifies the percentage of simulated reads that are expected to be 
from the target regions. Because of off-target effect in exome capture 
sequencing, this option is used to adjust the read number such that the we can 
obtain the desired number of on-target reads. This option is specific to each 
capture probe/target region. Users need to be careful changing this option. 

##### --read_length 

This option specifies the length of simulated short reads. 

##### --single_end

This option allows users to simulate single-end reads (instead of paired-end 
reads).

##### -d/--tumor_rdepth 

This option specifies the mean depth of the tumor sample. The number of 
generated short reads is computed by the following formula: 
(target_size\*rdepth)/(rlen\*ontarget_ratio). Here, 'target_size' is the size of 
target regions specified in the target file; 'rlen' is the simulated read length 
for the short reads (single-end or pair-end sequencing); 'target_ratio' is the 
percentage of on-target reads and specified by `--target_ratio`.

##### -r/--tumor_rnum

Instead of specifying the coverage of the target region, fa2wes also allows user 
to directly specify the number of simulated reads. This option specifies the 
number of short reads to simulate for tumor sample. 

##### -D/--normal_rdepth 

This option specifies the mean depth of the normal sample. In principle, one can 
specify the mean depth for tumor and normal samples in a single command. To save 
time, one can simulate tumor and normal samples in parallel by running two 
fa2wes commands at the same time, with one command specifying the depth of tumor 
(normal) sample to be 0 when simulating normal (tumor) sample.

##### -R/--normal_rnum 

This option allows users to specify the number of short reads to simulate for 
the normal sample. 

##### --simulator 

This option specifies the WES simulator to use. Available choices of the WES 
simulators include: wessim (default) and capgem.

##### --snakemake

This option specifies the snakemake command used to call an external WES 
simulator. Users can pass all parameters to a selected simulator by revising the 
Snakefile of the corresponding simulator (under directory wes/config/), except 
the option specifying read length. To accelerate the simulation process, users 
can submit each job to a cluster. For example, one can use the following option 
to submit the jobs to a Univa Grid Engine queuing system: `--snakemake 
'snakemake --rerun-incomplete -k --latency-wait 120 --cluster-config 
config/cluster.yaml --cluster "qsub -V -l 
mem_free={cluster.mem},h_rt={cluster.time} -pe OpenMP {cluster.n} -o stdout/ 
-e stdout/"'`. 

##### --cores

This option specifies the number of cores used to run the program (including 
snakemake). If `--cores` or `--jobs` or `-j` is specified in the options of 
snakemake command, the value specified by `--cores` here will be ignored when 
snakemake is called. Because WES simulation involves multiple steps, it is 
recommended to use a larger number of cores to save time.

##### --single

With this option, WES simulation will be in single cell mode, which simulates 
single cell exome sequencing. Specifically, the value of --tumor_depth 
(--tumor_rnum) refers to the depth (read number) of each tumor cell. 

##### --separate

With this option, fa2wes will keep the short reads of each genome separately 
(instead of mixing them together to create the tumor sample). 

##### --out_level [0,1,2]

This option specifies the level used to indicate how many intermediate output 
files are kept. 

- **Level 0**: keep all the files.  
- **Level 1**: keep files that are necessary for rerunning simulation ('config', 
'genome_index', 'mapping', 'merged', and 'separate'). 
- **Level 2**: keep only final results ('merged' and 'separate' folders).

### 2.6 allinone (module 6)

allinone is a convenient wrapper for simulating NGS reads for the normal and 
tumor sample in one command. It first calls module vcf2fa to construct the 
germline genome. Subsequently, it calls module phylovar and chain2fa to build 
the genomes of tumor cells. At the end, it calls module fa2wgs or fa2wes to 
simulate short reads for the normal and tumor sample. 

#### 2.6.1 Input files

##### Reference file (-r/--reference) and VCF file (-v/--vcf)

These two files will be passed to module vcf2fa to construct the germline genome 
of an individual.  Check section 2.1 for details.

##### Tree file (-t/--tree), configuration file (-c/--config) and affiliation file (--affiliation)

These three files will be passed to module phylovar to simulate somatic variants 
along the evolutionary history of the sample. Check section 2.2 for details.

##### Sectors file (--sectors)

This file will be passed to module fa2wgs and fa2wes. Check section 2.4 for details. 

##### --probe/--target/--error_model

These files will be passed to module fa2wes. Check section 2.5 for details.

#### 2.6.2 Output files

##### Output (--output)

The output folder is specified via `--output`. All files generated in the 
simulation will be placed in this folder. allinone will create four or five 
subfolders under this folder: normal_fa, tumor_chain, tumor_fa, wgs_reads and/or 
wes_reads. They contain normal FASTA files, tumor chain files, tumor FASTA files 
and short reads for normal/tumor samples respectively. 

##### Log file (-g/--log)

The log file is specified via `-g/--log`, which records random seeds and all the 
sub-commands, along with their parameters. 

#### 2.6.3 Options

##### --autosomes/--sex_chr

These two options will be passed to module vcf2fa. Check section 2.1 for details.

##### --trunk_vars/--trunk_length/--prune

These four options will be passed to module phylovar. Check section 2.2 for 
details.

##### --cores

This option will be passed to module chain2fa, fa2wgs and fa2wes. Check section 
2.3, 2.4 and 2.5 for details.

##### --purity/--rlen/--separate/--single

These options will be passed to module fa2wgs or fa2wes. Check section 2.4 and 
2.5 for details.

##### --art/--tumor_depth/--normal_depth

These options will be passed to module fa2wgs. Check section 2.4 for details.

##### --tumor_rdepth/--tumor_rnum/--normal_rdepth/--normal_rnum/--simulator/--snakemake/--out_level

These options will be passed to module fa2wes. Check section 2.5 for details.

##### --random_seed

This option specifies the seed used for random number generation. The generated 
random number will be used as a random seed by module phylovar, fa2wgs and 
fa2wes.

##### --start

The option `--start step_ID` allows users to skip some of the steps which have 
been finished already. The program will carry out the simulation from the 
step_ID th step. The step_ID (inclusive) can be an integer number from 1 to 4.

- **1**: vcf2fa
- **2**: phylovar
- **3**: chain2fa
- **4**: fa2wgs/fa2wes

## 3. Tutorial

In this section, we will demonstrate how to use PSiTE.

1. At first, several input files have to be prepared with the following steps.

    - Simulate a coalescent tree of 1000 tumor cells, which are sampled from 
    an exponentially growing tumor. (Please check the manual of ms for more 
    information)
  
        ```
        ms 1000 1 -T -G 1 |tail -n1 > ms_tree.txt
        ```
  
    - Download the fasta file of human reference genome from the website of 1000 
    genomes.
  
        ```
        wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz
        gunzip human_g1k_v37.fasta.gz
        ```
  
    - Download variants data of NA12878 (Genome in a bottle consortium) from 
    NCBI.
  
        ```
        wget -O NA12878.raw.vcf.gz ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/latest/GRCh37/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz
        ```
  
    - Filter the raw variants to get the phased SNPs as the germline variants of 
    the sample to simulate.
  
        ```
        zcat NA12878.raw.vcf.gz \
            |awk '/^#/ || ($NF~/^[01]\|[01]/ && length($4)==1 && length($5)==1)' \
            |gzip -c > NA12878.phased_snp.vcf.gz
        ```
    
2. Next, vcf2fa simulates the (normal) germline genomes of a female individual, 
by integrating germline SNPs into the human reference genome.

      ```
      psite.py vcf2fa -r human_g1k_v37.fasta -v NA12878.phased_snp.vcf.gz \
          --autosomes 1..22 --sex_chr X,X -o normal_fa
      ```

3. Subsequently, phylovar simulates somatic variants of the sample. There are 
multiple options in phylovar that allow a flexible simulation, as shown below.

    - Simulate variants based on a configuration file. Check 
    cfg_template_female.yaml in PSiTE package (under folder example_doc)  for 
    detailed settings.
  
        ```
        psite.py phylovar -t ms_tree.txt --config cfg_template_female.yaml \
            --purity 0.8 --sex_chr X,X
        ```
  
    - With the option --trunk_length, we can simulate truncal mutations of the 
    tumor sample.
  
        ```
        psite.py phylovar -t ms_tree.txt --config cfg_template_female.yaml \
            --purity 0.8 --sex_chr X,X --trunk_length 2.0
        ```
  
    - The above two phylovar commands will simulate somatic variants of 1000 
    tumor genomes. The fasta files of these genomes can consume a lot of space, 
    because we need to store the genomes of every individual cell. Since most 
    low frequency variants are not informative, users can employ `--prune` to 
    trim the infrequent lineages.
  
        ```
        psite.py phylovar -t ms_tree.txt --config cfg_template_female.yaml \
            --purity 0.8 --sex_chr X,X --trunk_length 2.0 --prune 0.05
        ```
  
    - When running phylovar, users can choose to generate chain files and map 
    files with option `--chain` and `--map` respectively. The command will 
    generate the chain files for each tip node under the folder tumor_chain. 
    These two files can then be used for generating genomes of tumor cells and 
    the sequencing data from those cells.
  
        ```
        psite.py phylovar -t ms_tree.txt --config cfg_template_female.yaml \
            --purity 0.8 --sex_chr X,X --trunk_length 2.0 --prune 0.05 \
            --chain tumor_chain --map map
        ```
  
4. With the chain files generated by phylovar, chain2fa can then build the 
genomes of tumor cells in the sample. 

    ```
    psite.py chain2fa -c tumor_chain \
        -n normal_fa/normal.parental_0.fa,normal_fa/normal.parental_1.fa \
        -o tumor_fa --cores 16
    ```

5. After generating the tumor genomes, NGS reads can be simulated by calling 
fa2wgs or fa2wes. Below are several examples showing how to generate different 
types of NGS reads.
    - Simulate the NGS reads of a tumor sample with the purity of 0.8 and the 
    coverage of 50X. At the same time, generate the paired normal sample at 
    coverage 30X. Use 16 CPUs to run the simulation.
  
        ```
        psite.py fa2wgs -n normal_fa -t tumor_fa -m map --purity 0.8 \
            --tumor_depth 50 --normal_depth 30  -o wgs_reads --cores 16
        ```
  
    - Simulate the WES reads of a tumor sample with the purity of 0.8 and the 
    coverage of 100X as well as the paired normal sample with 100X coverage. 
    Use wessim and 16 CPUs to run the simulation. The probe file, target file 
    and error model file used in this command can be found under directory 
    wes/example of PSiTE package (**Note:** The probe file 
    PSiTE/wes/example/S03723314_Probes.fa is compressed in '.gz' format when 
    distributing. Users should extract it before using).
  
        ```
        psite.py fa2wes -n normal_fa -t tumor_fa -m map \
            --probe PSiTE/wes/example/S03723314_Probes.fa \
            --target PSiTE/wes/example/S03723314_Covered_c3.bed \
            --error_model PSiTE/wes/example/RMNISTHS_30xdownsample_chr22_p.gzip \
            --purity 0.8 --tumor_rdepth 100 --normal_rdepth 100 \
            --simulator wessim -o wes_reads --cores 16
        ```
  
    - Simulate the WES reads of a tumor sample with the purity of 0.8 and the 
    coverage of 100X. Use capgem and clusters to run the simulation. Please ensure 
    that there are enough memory and disk space to run capgem (see section 2.5.0 
    for the requirements).
  
        ```
        psite.py fa2wes -n normal_fa -t tumor_fa -m map \
            --probe PSiTE/wes/example/S03723314_Probes.fa \
            --target PSiTE/wes/example/S03723314_Covered_c3.bed \
            --error_model PSiTE/wes/example/RMNISTHS_30xdownsample_chr22_p.gzip \
            --purity 0.8 --tumor_rdepth 100 --normal_rdepth 0 \
            --simulator capgem -o wes_reads \
            --snakemake 'snakemake -j 200 --rerun-incomplete --latency-wait 120 \
            -k --cluster "qsub -V -l mem_free={cluster.mem},h_rt={cluster.time} \
            -pe OpenMP {cluster.n} -o wes_reads/stdout/ -e wes_reads/stdout/"' 
        ```
  
6. Finally, users can just use a single allinone command to run the whole 
pipeline to generate the WGS or WES data. 

    - Run the whole pipeline to generate the WGS data.
  
        ```
        psite.py allinone \
            -r human_g1k_v37.fasta -v NA12878.phased_snp.vcf.gz \
            -t ms_tree.txt -c cfg_template_female.yaml -o output \
            --autosomes 1..22 --sex_chr X,X --trunk_length 2.0 \
            -d 50 -D 30 -p 0.8 -x 0.05 --cores 16
        ```
      
    - Run the whole pipeline to generate the WES data.
  
        ```
        psite.py allinone --type WES \
            -r human_g1k_v37.fasta -v NA12878.phased_snp.vcf.gz \
            -t ms_tree.txt -c cfg_template_female.yaml -o output \
            --autosomes 1..22 --sex_chr X,X --trunk_length 2.0 \
            -p 0.8 -x 0.05 \
            --probe PSiTE/wes/example/S03723314_Probes.fa \
            --target PSiTE/wes/example/S03723314_Covered_c3.bed \
            --error_model PSiTE/wes/example/RMNISTHS_30xdownsample_chr22_p.gzip \
            --rlen 150 --tumor_rdepth 100 --normal_rdepth 100 \
            --simulator wessim --cores 16
        ```
  
    - Run the whole pipeline to generate both WGS and WES data.
  
        ```
        psite.py allinone --type BOTH \
            -r human_g1k_v37.fasta -v NA12878.phased_snp.vcf.gz \
            -t ms_tree.txt -c cfg_template_female.yaml -o output \
            --autosomes 1..22 --sex_chr X,X --trunk_length 2.0 \
            -d 50 -D 30 -p 0.8 -x 0.05 \
            --probe PSiTE/wes/example/S03723314_Probes.fa \
            --target PSiTE/wes/example/S03723314_Covered_c3.bed \
            --error_model PSiTE/wes/example/RMNISTHS_30xdownsample_chr22_p.gzip \
            --rlen 150 --tumor_rdepth 100 --normal_rdepth 100 \
            --simulator wessim --cores 16
        ```

## Authors

* [Hechuan Yang](https://github.com/hchyang)
* [Bingxin Lu](https://github.com/icelu)

## License

This project is licensed under the GNU GPLv3 License - see the 
[LICENSE](LICENSE) file for details.

