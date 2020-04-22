

# StrainSifter

This is a NextFlow translation of Fiona Tamburini's [StrainSifter workflow](https://github.com/tamburinif/StrainSifter/)

## Installation





### Dependencies

We recommend running StrainSifter in the provided conda environment. If you wish to run StrainSifter without using the conda environment, the following tools must be installed and in your system PATH:
* [Burrows-Wheeler Aligner (BWA)](http://bio-bwa.sourceforge.net)
* [Samtools](http://www.htslib.org)
* [Bamtools](https://github.com/pezmaster31/bamtools)
* [Bedtools](http://bedtools.readthedocs.io/en/latest/)
* [MUSCLE](https://www.drive5.com/muscle/)
* [FastTree](http://www.microbesonline.org/fasttree/)
* [Python3](https://www.python.org/downloads/)
* [R](https://www.r-project.org)

## Running StrainSifter

Due to the computing demands of the StrainSifter pipeline, we recommend running on a computing cluster if possible.
Instructions to enable Snakemake to schedule cluster jobs with SLURM can be found at https://github.com/bhattlab/slurm

### Input files

* Reference genome assembly in fasta format (can be a draft genome or a finished reference genome)
Acceptable file extensions: ".fasta", ".fa", ".fna"

* Two or more short read datasets in fastq format (metagenomic reads or isolate reads), optionally gzipped
Acceptable file extensions: ".fq", ".fastq", ".fq.gz", ".fastq.gz"

Short read data can be paired- or single-end.

You will need to indicate input files in the config file for each sample you wish to run StrainSifter on. This is described below in the *Config* section:

### Config

You must update the ssifter.config ss file as follows:

**reference:** Path to reference genome (fasta format)

**reads:** Samples and the file path(s) to the input reads.

<br>

Optionally, you can update the following parameters:

**prefix:** (optional) desired filename for output files. If blank, the name of the reference genome will be used.

**mapq:** minimum mapping quality score to evaluate a read aligment

**n_mismatches:** consider reads with this many mismatches or fewer

**min_cvg:** minimum read depth to determine the nucleotide at any given postion

**min_genome_percent:** the minimum fraction of bases that must be covered at min_cvg or greater to process an sample

**base_freq:** minimum frequency of a nucleotide to call a base at any position

<br>

Example config.yaml:

    ##### input files #####

    # reference genome (required)
    reference: /path/to/ref.fna

    # short read data (at least two samples required)
    reads:
    sample1:
    [/path/to/sample1_R1.fq,
    /path/to/sample1_R2.fq]
    sample2:
    [/path/to/sample2_R1.fq,
    /path/to/sample2_R2.fq]
    sample3: /path/to/sample3.fq

    # prefix for output files (optional - can leave blank)
    prefix:


    ##### workflow parameters #####

    # alignment parameters:
    mapq: 60
    n_mismatches: 5

    # variant calling parameters:
    min_cvg: 5
    min_genome_percent: 0.5
    base_freq: 0.8


### Running StrainSifter

To run StrainSifter, the config file must be present in the directory in which you wish to run the workflow.
You should then be able to run StrainSifter as follows:

#### Phylogeny

To generate a phylogenetic tree showing all of the input samples that contain your strain of interest at sufficient coverage to profile:

    snakemake {prefix}.tree.pdf

#### SNV counts

To generate a list of pairwise SNV counts between all input samples:

    snakemake {prefix}.dist.tsv

### FAQ

Q: Can StrainSifter be used for non-bacterial genomes (e.g. yeast)?

A: At present, we recommend StrainSifter for bacteria only. In theory, StrainSifter should work for yeast if a haploid reference genome is provided.
