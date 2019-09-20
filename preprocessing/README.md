# h3ameta - preprocessing

## Purpose

This is a simple workflow that does QC. Input data is given as a set of fastq file and QC is done by

* running fastqc and producing reports
* running multiqc and producing reports
* running Trimmomatic
* re-rerunnning fastqc and multiqc on the trimmed data


## Running the workflow


The basic way to run the workflow is

```
nextflow run preprocessing.nf
```

## Config file 

The following set up needs to be done in the `nextflow.config` file

* The Trimmomatic jar file must be specified by setting the `trimmomatic_jar` variable. The default value
  is the location in the Docker image that we provide. If you use our Docker image using Docker or Singularity 
  then you do not need to change `trimmomatic_jar`. If you do not use our container set the variable -- full, absolute   path inclusing file name

* The queue name you want to run your job on if you are using slurm or pbs -- otherwise you can ignore.

* give the names of the input data -- as a Unix-style glob. Note that
  * if you are doing paired-end QC, the glob must be of this form
     `data/dir1/dir2/*_{1,2}.fq.gz`. The key part is the star (*) and
      the braces { } -- file suffix can be anything sensible and does
       not have to be compressed
  * `paired` : state whether paired-end QC (true) or single end QC should be done (false)
  * `out_dir` : where the directory that the results will go to: Thre will be four separate directories  `fastqc`,  `fastqc_trim`,  `multiqc` and   `multiqc_trim`
  * `num_cpus` : how many CPUs you want for `fastqc`. This is per-file figure, so a modest number is probably adequate.
* Then specify the [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) parameters that you want to use. Sensible default parameters are given but you should really apply your mind.