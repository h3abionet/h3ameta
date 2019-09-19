# Still under development

## Pipeline for short reads
A Nextflow pipeline script for processing metagenomic short reads (Illumina). The nextflow includes cleaning reads from human DNA, aligning against a viral genomes and taxonomic classification of the reads. At the end of the pipeline, an HTML report will be generated comparing the results of the statistics of the alignment against the candidate viral sequences and the results of the taxonomic classification. The report provides an interactive pie chart for  taxa investigation.


#Software requirements would be:

* [Nextflow](https://www.nextflow.io/)
* [Singularity](https://www.sylabs.io/guides/3.0/user-guide/installation.html)
* [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
* [Samtools](http://www.htslib.org/download/)
* [Kraken2](https://ccb.jhu.edu/software/kraken2/)
* [Pairfq](https://github.com/sestaton/Pairfq)
* [Fastq screen](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/_build/html/index.html)
* [Krona](https://github.com/marbl/Krona/wiki/Installing)
* [Python2.7](https://www.python.org/download/releases/2.7/)

# Python libraries requirements would be: 

pip2 install os-sys
pip2 install regex
pip2 install glob2
pip2 install python-csv
Pip2 install pandas 
Pip2 install numpy 
Pip2 install argparse
Pip2 install shutil

# Databases requirements would be:
[Viral genomes] : ???
[Minikraken2](https://ccb.jhu.edu/software/kraken2/index.shtml?t=downloads)
[Human genome](ftp://ftp.ensembl.org/pub/current/fasta/homo_sapiens/dna/)


# Configuration
1.The pipeline works with 2 FASTQ files  2 input files (for forward and reverse reads). Each paired end reads need to have a common basename. to differentiate between forward and reverse, the two files must contain the consecutive prefix "R1.fastq" and "R2.fastq".

2. The human genome and viral genomes  need  to be formatted using bowtie2-build:
	Example : 
  ```
 bowtie2-build all.fasta viralgenomes
 bowtie2-build GRCh38  humangenome
```
3. Set the path to the human reference database in `fastq_screen.conf`
4. All other parameters can be configured in `nextflow.config` or as parameters when running.

## Running test data
```
git clone https://github.com/h3abionet/h3ameta.git
cd h3ameta/viraldetect/short-reads/
nextflow main_illumina.nf --reads ‘test_data/*_{1,2}.fastq' -c template_nextflow.config
```

## Example reads are here

1) Forward reads are [here](http://web.cbio.uct.ac.za/~gerrit/downloads/viraldetect/illumina/gut/in/ERR1600426_1.100.fastq.gz)
2) Reverse reads are [here](http://web.cbio.uct.ac.za/~gerrit/downloads/viraldetect/illumina/gut/in/ERR1600426_2.100.fastq.gz)
3) Pipeline output is [here](http://web.cbio.uct.ac.za/~gerrit/downloads/viraldetect/illumina/gut/out)

## Running with Docker

Alternatively, Viral detect pipeline can be run using a Docker image. 
At the first step, Docker need to be installed. Installation instructions can be found [here](https://docs.docker.com/install/)


## To note
* Bowtie2 databases need to be build beforehand. Should include this process in the pipeline at some point. Need to give current instructions here.
