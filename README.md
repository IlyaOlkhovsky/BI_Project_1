# Bioinformatics Institute Project 1: mutations identification
This tool is developed to ease the analysis of mutations in bacterial genomes.
## Requeirements
The following packages need to be installed for our tool to work correctly:
 - [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
 - [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
 - [BWA](https://bio-bwa.sourceforge.net/)
 - [SAMtools](http://www.htslib.org/)
 - [VarScan](http://dkoboldt.github.io/varscan/)
## Usage
The input data for our tool consists of a single file in `.fna` format with a reference genome and two files in `.fastq` format with sequencing data (reads in forward and reverse direction). It is assumed that the sequenced data is obtained with [ILLUMINA](https://www.illumina.com/). Our tool is used as follows:
```
 python3 main.py [-h] [--trimmomatic TRIMMOMATIC] [--varscan VARSCAN] 
 --ref_genome REF_GENOME --reads_forward READS_FORWARD --reads_reverse READS_REVERSE
```
The result of the latter command is the `.vcf` files with the VarScan output, which can be found in the `varscan_output` folder. These can then be used with [IGV](https://software.broadinstitute.org/software/igv/) to view the sequenced reads alignment on the reference genome and putative mutations.
