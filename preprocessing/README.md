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
nextflow run main.nf
```

If you want to run it using Singularity

```
nextflow run main.nf -profile singularity
```

If you want to use it using Docker


```
nextflow run main.nf -profile docker
```

If you want to use it using Slurm

```
nextflow run main.nf -profile slurm
```

If you want to use it using Slurm and Docker

```
nextflow run main.nf -profile docker,slurm
```






## Config file 

Configuration is by default done in the `nextflow.config` file. However, see our recommended way of running later

* The Trimmomatic jar file must be specified by setting the `trimmomatic_jar` variable. The default value
  is the location in the Docker image from bioicontairs. If you use this Docker image using Docker or Singularity 
  then you do not need to change `trimmomatic_jar`. If you do not use our container set the variable -- full, absolute   path including file name.

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


## Recommended configuration

Nextflow supports _multiple_ config file. We recommend that the `nextflow.config` file that we provide should be minimally changed so as to support running at your site. Typically this would mean changing the `trimmomatic_jar` and `job_queue` variables and that't all. You could even delete the `input` variable. This could then be used by multiple users and/or over multiple runs.

Then for each run, create a separate  config file, which is very simple e.g. `run17.config`.

```
params {
       paired = false
       input    = "quality_control/*_{1,2}.fq.gz"
    
}   
```

You can then run your workflow


`nextflow run -c run17.config -c nextflow.config  main.nf`


What this does is take your `nextflow.config` file as a base and then add or _change_ any entries that are found in `run17.config`. 

