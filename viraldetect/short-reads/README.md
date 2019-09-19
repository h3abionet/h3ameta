# Still under development

## Pipeline for short reads
A Nextflow pipeline script for processing metagenomic short reads (Illumina). The nextflow includes cleaning reads from human DNA, aligning against a viral genomes and taxonomic classification of the reads. At the end of the pipeline, an HTML report will be generated comparing the results of the statistics of the alignment against the candidate viral sequences and the results of the taxonomic classification. The report provides an interactive pie chart for  taxa investigation.


## Software requirements would be:

* [Nextflow](https://www.nextflow.io/)
* [Singularity](https://www.sylabs.io/guides/3.0/user-guide/installation.html)

## Configuration
1. Set the path to the human reference database in `fastq_screen.conf`
2. All other parameters can be configured in `nextflow.config` or as parameters when running.

## Running

```
nextflow run main.nf
```
## Example reads are here

1) Forward reads are [here](http://web.cbio.uct.ac.za/~gerrit/downloads/viraldetect/illumina/gut/in/ERR1600426_1.100.fastq.gz)
2) Reverse reads are [here](http://web.cbio.uct.ac.za/~gerrit/downloads/viraldetect/illumina/gut/in/ERR1600426_2.100.fastq.gz)
3) Pipeline output is [here](http://web.cbio.uct.ac.za/~gerrit/downloads/viraldetect/illumina/gut/out)

## To note
* Bowtie2 databases need to be build beforehand. Should include this process in the pipeline at some point. Need to give current instructions here.
