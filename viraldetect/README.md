
# Viral Detection Pipeline

Thanks to the increased cost-effectiveness of high-throughput technologies, the number of studies focusing on microorganisms (bacteria, archaea, microbial eukaryotes, fungi, and viruses) and on their connections with human health and diseases has surged, and, consequently, a plethora of approaches and software has been made available for their study, making it difficult to select the best methods and tools. 

Here we present a viral detection pipeline that, starting from the raw sequencing data and having a strong focus on quality control, allows, within hours, the data processing up to the functional annotation.

It is constructed on [Nextflow](https://github.com/nextflow-io/nextflow), a framework based on the dataflow programming model, which allows writing workflows that are highly parallel, easily portable (including on distributed systems), and very flexible and customisable. 

Viral dection pipeline is accompanied by a [Docker container](https://www.docker.com/), that saves the users from the hassle of installing the required software, increasing, at the same time, the reproducibility of the pipeline results (see [Using Docker or Singularity](#using-docker-or-singularity)). 



## Table of contents

- [Citation](#Citation)
- [Dependencies](#dependencies)
- [Installation](#installation)
- [Other requirements](#other-requirements)
- [Usage](#usage)
- [Using Docker or Singularity](#using-docker-or-singularity)
- [Troubleshooting](#Troubleshooting)
- [Changelog](#changelog)
- [License](#license)
- [Acknowledgements](#acknowledgements)


## Citation

Please cite viral Detection as:

> citation goes here

## Dependencies

There is need install Nextflow  (version 0.26.x or higher), as explained [here](https://www.nextflow.io/docs/latest/getstarted.html). Please note that Nextflow requires BASH and [Java 7](http://www.oracle.com/technetwork/java/javase/downloads/index.html) or higher to be installed. Both should be already available in most of the POSIX compatible systems (Linux, Solaris, OS X, etc). However, as of October 2017, the latest release of Java (SE9) introduces some breaking changes in Nextflow, and should not be used (see [here](https://github.com/nextflow-io/nextflow/issues/462) for details). 

If you are using the containerised version of the pipeline (as we strongly suggest), you will should also install [Docker](https://www.docker.com) or [Singularity](http://singularity.lbl.gov/), as explained [here](https://docs.docker.com/engine/installation/) and [here](http://singularity.lbl.gov/docs-installation), respectively.
In fact, Nextflow orchestrates, in a transparent fashion, the flow of the pipeline by wrapping and executing each step using the Docker/Singularity run command. Thus, Nextflow lies *outside* the container, that is responsible for instantiating. 
You can find more information about Docker/Singularity containers and Nextflow [here](https://www.nextflow.io/docs/latest/docker.html) and [here](https://www.nextflow.io/docs/latest/singularity.html), respectively.

Once you have either Docker or Singularity up and running, you will not need to install anything additional tools, since all the pieces of software are already available in the Docker container released with YAMP pipeline, and that you can find on [DockerHub](https://hub.docker.com/r/alesssia/viraldetect/). Please refer to [Using Docker or Singularity](#using-docker-or-singularity) for more details. 

The list of tools that should be available includes:
- Samtools v1.3.1 ([http://samtools.sourceforge.net](http://samtools.sourceforge.net))

Following the links, you will find detailed instructions on how to install them, as explained by their developers. 
Notably, MetaPhlAn2, QIIME, and HUMAnN2 are also available in [bioconda](https://anaconda.org/bioconda/). 


## Installation

Clone the pipeline repository in a directory of your choice:

```
git clone https://github.com/h3abionet/h3ameta.git
```

The repository includes:

- the Nextflow script, `main.nf`, 
- the configuration files, `nextflow.config`
- a folder (`bin`) containing two helper scripts,
- a folder (`viraldetect`) containing the Docker file used to build the Docker image (`Dockerfile`). 

**Note:** the `nextflow.config` file includes the parameters that are used in our tutorials and testing sessions.


## Other requirements

The pipeline requires a set of databases that are queried during its execution which include. the following;

- Human genome
- kraken database
- Candidate viral genomes

## Usage

1. Modify the `nextflow.config` file, specifying the necessary parameters, such as the path to the aforementioned databases.
2. From a terminal window run the `main.nf` script using the following command (when the library layout is 'paired'):
	```
	nextflow run main.nf --reads1 R1 --reads2 R2  --outdir outputdir --mode MODE
	```
	where `R1` and `R2` represent the path to the raw data (two compressed paired-end FASTQ files), `mysample` is a prefix that will be used to label all the resulting files, `outputdir` is the directory where the results will be stored, and `MODE` is any of the following: < QC, characterisation, complete >; or  the following command (when the library layout is 'single'):
	

## Using Docker or Singularity

To use the tools made available through the Docker container within both Docker, one could either pull the pre-built image from [DockerHub](https://hub.docker.com/r/alesssia/viraldetect/), using the following command:

```
docker pull alesssia/viraldetect
```

or build a local image using the file `Dokerfile` in  the `viraldetect` folder. To build a local image, one should first access the `viraldetect` folder and then run the following command (be careful to add the dot!):

```
docker build -t viraldetect .
```

In both cases, the image can be used by YAMP by running the command presented above adding `-with-docker` followed by the image name (`viraldetect`):

```
nextflow run main.nf --reads1 R1 --reads2 R2 --outdir outputdir --mode MODE -with-docker viraldetect
```

where `R1` and `R2` represent the path to the raw data (two compressed FASTQ file), `mysample` is a prefix that will be used to label all the resulting files, `outputdir` is the directory where the results will be stored, and `MODE` is any of the following: < QC, characterisation, complete >.

YAMP can also fetch the Docker container directly from DockerHub;

```
nextflow run main.nf --reads1 R1 --reads2 R2 --outdir outputdir --mode MODE -with-docker docker://alesssia/viraldetect
```

so, even simpler!

Viral detection pipeline can use a Docker image with Singularity (again without pulling the image) by adding the `-with-singularity` option followed by the image path (`--with-singularity docker://alesssia/viraldetect`), that is, the following command:

```
nextflow run main.nf --reads1 R1 --reads2 R2 --prefix mysample --outdir outputdir --mode MODE -with-singularity docker://alesssia/viraldetect
```


Please note that Nextflow is not included in the Docker container and should be installed as explained [here](https://www.nextflow.io/docs/latest/getstarted.html).


## Troubleshooting

We have listed all known issues and solutions on this [wiki page](https://github.com/alesssia/YAMP/wiki/Troubleshooting). Please report any issue using the [GitHub platform](https://github.com/alesssia/YAMP/issues).


## Changelog

## License

YAMP is licensed under GNU GPL v3.


## Acknowledgements


