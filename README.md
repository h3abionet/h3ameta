
<img src="docs/assets/images/H3ABioNetlogo2.jpg"/>

# h3ameta
`h3ameta` is a collection of **H3ABionNet Metagenomics Workflows** written in `nextflow`:
1. `TaxonomicClassification`
2. `StrainComp`
3. `ViralDetectLong`
4. `ViralDetectShort`
5. ~~`StrainSifter`~~ (in development)

## 1. Setting up the `h3ameta` workflow:
### 1.1. Download the worklfow:
```console
netflow pull h3abionet/h3ameta
```
To get the help menu:
```console
nextflow run h3ameta -r phele --help
```

### 1.2. Download the test datasets:
**NB:** *I havent found a place to put the test dataset (~3.4GB)*
```console
wget <link>
```

### 1.3. Downolad workflow Singularity containers:
```console
nextflow run h3ameta -r phele -profile slurm --mode prep.Containers 
```

### 1.4. Generate `kraken2` database:
**NB:** Installation of the `kraken2` database is time consuming, computationaly intensive and requires a lot of space. Before installing the `kraken2` database, please check with your system administrator if there is an instance of `kraken2` database already installed on the system. If the `krakend2` database is installed, skip this step and use the `--kraken_db` to specify the location of the installed `kraken2` instance throughout this walk-through of the `h3ameta` pipeline.

Should you wish to install your own instance of the database, use the `--kraken_db` in the command below to specify the location where you would like to install the database; if not specified, the database will be installed in the current working directory.
```console
nextflow run h3ameta -r phele -profile slurm --mode prep.KrakenDB
```

### 1.5. Generate `braken` database:
**NB:** Installation of the `braken` database is dependent on `kraken2` library installation files. Again, please check with your system administrator if the dattabase already exists (in principle, this **SHOULD** be in the `kraken2` database location). If the installattion of `braken` database exist (and in the correct `kraken2` database location), please skip this step.

If you ran the command above to build the `kraken2` database, you may proceed witth the command below, specifying the same `--kraken_db` path as in tthe command above.
```console
nextflow run h3ameta -r phele -profile slurm --mode prep.BrakenDB
```

## 2. Using the `h3ameta` workflow:
Once the `h3ameta` workflow has been setup (dataset, `Singularity` containers/images, `kraken2` database and `braken` database), you're ready to use the workflow. There are two ways parameters can be passed to the pipelin: (1) using a config file, which is passed to the pipeline using the `-c` option; and (2) passing options directly on the command line. This walk-through uses both options. The config files are located in the `data/confs/` directory of the test dataset downloaded in **1.2.** (have a look at them).

### 2.1. Data QC (optional)
The data QC step is optional. It is for assessing quality of your reads and removing low quality bases and contaminating adapters.
#### 2.1.1. Read QC with `fastqc`:
```console
## Using a configuration file
nextflow run h3ameta -r phele -profile slurm --mode run.ReadQC -c data/confs/read_qc.conf

## Passing command-line arguements
nextflow run h3ameta -r phele -profile slurm --mode run.ReadQC \
    --data $PWD/data/reads/short \
    --out output
```

#### 2.1.2. Read trimming with `trimmomatic`
```console
## Using a configuration file
nextflow run h3ameta -r phele -profile slurm --mode run.ReadTrimming -c data/confs/read_trimming.conf

## Passing command-line arguements
nextflow run h3ameta -r phele -profile slurm --mode run.ReadTrimming \
    --data $PWD/data/reads/short \
    --out output \
    --trim "ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:8:keepBothReads TRAILING:28 MINLEN:4"
```

### 2.2. Workflow 1: `TaxonomicClassification`
```console
## Using a configuration file
nextflow run h3ameta -r phele -profile slurm --mode run.Classification -c data/confs/classification.conf

## Passing command-line arguements
nextflow run h3ameta -r phele -profile slurm --mode run.Classification \
    --data $PWD/data/reads/short \
    --out output
```

### 2.3. Workflow 2: `StrainComp`
```console
## Using a configuration file
nextflow run h3ameta -r phele -profile slurm --mode run.StrainComp -c data/confs/strain_comp.conf

## Passing command-line arguements
nextflow run h3ameta -r phele -profile slurm --mode run.StrainComp \
    --data $PWD/data/reads/straincomp \
    --out output \
    --kraken_db /global/blast/KDB \
    --dataset_table $PWD/data/reads/straincomp/datasets.tsv \
    --readlen 100 \
    --tax_level S \
    --annot_db $PWD/data/dbs/argannot_db/ARGannot_r2.fasta \
    --singleEnd true
```

### 2.4. Workflow 3: `ViralDetectLong`
```console
## Using a configuration file
nextflow run h3ameta -r phele -profile slurm --mode run.ViralDetectLong -c data/confs/viral_detect_long.conf

## Passing command-line arguements
nextflow run h3ameta -r phele -profile slurm --mode run.ViralDetectLong \
    --data $PWD/data/reads/long/nanopore \
    --out output \
    --kraken_db /global/blast/KDB \
    --read_type nanopore \
    --viral_db $PWD/data/dbs/viral_db/all.fasta \
    --decontam_db $PWD/data/dbs/decontam_db/GCF_000001405.38_GRCh38.p12_genomic.fna.gz \
    --acc_2_tax $PWD/data/dbs/ncbi_db/test_data/nucl_gb.accession2taxid \
    --names_dmp $PWD/data/dbs/ncbi_db/names.dmp \
    --singleEnd true
```

### 2.5. Workflow 4: `ViralDetectShort`
```console
## Using a configuration file
nextflow run h3ameta -r phele -profile slurm --mode run.ViralDetectShort -c data/confs/viral_detect_short.conf

## Passing command-line arguements
nextflow run h3ameta -r phele -profile slurm --mode run.ViralDetectShort \
    --data $PWD/data/reads/short \
    --out output \
    --kraken_db /global/blast/KDB \
    --viral_db $PWD/data/dbs/viral_db/all.fasta
```

## 3. `h3ameta` workflow output:
### 3.1. Data QC (optional)
#### 3.1.1. Read QC with `fastqc`:
#### 3.1.2. Read trimming with `trimmomatic`
### 3.2. Workflow 1: `TaxonomicClassification`
### 3.3. Workflow 2: `StrainComp`
### 3.4. Workflow 3: `ViralDetectLong`
### 3.5. Workflow 4: `ViralDetectShort`

<!-- Note: other workshop materials can be found [in our Google Drive folder](https://drive.google.com/drive/u/1/folders/1g3iyBbbD0fq2TIYz3MungaOiSu4DAm8X) -->

<!-- ## Running the model workflow -->

<!-- ### 1. Set up conda, nextflow, clone the Git repository. -->

<!-- Note: this requires Singularity to be set up on your system or cluster. -->

<!-- ``` -->
<!-- cd ~ -->
<!-- mkdir -p ~/local/bin -->
<!-- export PATH="$PATH:~/local/bin" -->

<!-- wget -qO- https://get.nextflow.io | bash -->
<!-- wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -->
<!-- cd -->
<!-- git clone https://github.com/h3abionet/h3ameta.git -->
<!-- ``` -->


<!-- ### 2. Running the workflow -->

<!-- ``` -->
<!-- cd ~ -->
<!-- mkdir test_run; cd test_run -->
<!-- nextflow h3ameta/examples/taxonomic_classification/taxonomic_classification.nf  --tax_level S -resume --in h3ameta/examples/test_data/*.fq \ -->
<!-- --dataset_table h3ameta/examples/test_data/datasets.tsv --db /path/to/kraken_and_bracken_db -->
<!-- ``` -->

<!-- ## Docker images -->

<!-- We're assuming you're using singularity -- if using Docker it'll be a little simpler, so it's left as an exercise for the reader. Of course, if you're  using Nextflow this will generally be taken care of by the appropriate config file and should be transparent. -->

<!-- ### kraken2 -->

<!-- Download the latest image -->

<!-- `singularity pull docker://quay.io/h3abionet_org/kraken2 ` -->

<!-- This will create an image `kraken2.img` which contains the kraken2 suite plus auxiliary programs like dustmasker -->

<!-- Note that we do not have any databases inside the image to keep the image small. You need to download and build the databases. Here's an example: Assume that you have a directory `/local/kraken` and you're going to bulild the database inside that -->

<!-- ``` -->
<!-- singularity exec -B /local/kraken/:/mnt kraken2.simg kraken2-build --standard --threads 8 --db /mnt/krakdb -->
<!-- ``` -->
<!-- This binds the directory `/local/kraken` on the host to the `/mnt` directory in the singularity image. The directory `/mnt` is passed to the `kraken2-build` program to use for the data and the database will be called `krakdb`. -->
