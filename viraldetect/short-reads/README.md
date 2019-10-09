# Still under development

## Pipeline for short reads
Still in development. Working on Mariems's code to get it up and running. Basic structure in place. Need to test Nextflow code and build final container for report.

## Software requirements would be:

* [Nextflow](https://www.nextflow.io/)
* [Singularity](https://www.sylabs.io/guides/3.0/user-guide/installation.html)

## Configuration
1. All parameters can be configured in `nextflow.config` or as parameters when running.

## Running

```
nextflow run main.nf
```
## Example reads are here

1) Forward reads are [here](http://web.cbio.uct.ac.za/~gerrit/downloads/viraldetect/illumina/gut/in/ERR1600426_1.fastq.gz)
2) Reverse reads are [here](http://web.cbio.uct.ac.za/~gerrit/downloads/viraldetect/illumina/gut/in/ERR1600426_2.fastq.gz)
3) Pipeline output is [here](http://web.cbio.uct.ac.za/~gerrit/downloads/viraldetect/illumina/gut/out)

## To note
* Bowtie2 databases need to be build beforehand.
