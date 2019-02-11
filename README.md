
# h3ameta
H3ABionNet Metagenomics Workflows

Note: other workshop materials can be found [in our Google Drive folder](https://drive.google.com/drive/u/1/folders/1g3iyBbbD0fq2TIYz3MungaOiSu4DAm8X)

## Running the model workflow

### 1. Set up conda, nextflow, clone the Git repository.

Note: this require Singularity to be set up on your system or cluster.

```
cd ~
mkdir -p ~/local/bin
export PATH="$PATH:~/local/bin"

wget -qO- https://get.nextflow.io | bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
cd
git clone https://github.com/h3abionet/h3ameta.git
```

### Note: if singularity isn't supported on your compute cluster, set up environment manually instead.
```
cd ~/local/bin
bash Miniconda3*.sh #accept the defaults
conda install -y -c conda-forge -c bioconda -c r \
kraken2 krona kraken ncurses datrie r-ggplot2 r-doby r-rcolorbrewer r-scales r-plyr r-stringi
mkdir ~/miniconda/bin/taxonomy
ktUpdateTaxonomy.sh

git clone https://github.com/jenniferlu717/Bracken.git
cd Bracken
bash install_bracken.sh
cp bracken ~/local/bin/
cp bracken-build ~/local/bin/
```

### 2. Running the workflow

To use the singularity image, uncomment the last flag below.


```
cd ~
mkdir test_run; cd test_run
nextflow h3ameta/examples/taxonomic_classification/taxonomic_classification.nf  --tax_level S -resume --in h3ameta/test_data/*.fq #-with-singularity shub://bhattlab/wits_workshop:classification
```



## Docker images

We're assuming you're using singularity -- if using Docker it'll be a little simpler, so it's left as an exercise for the reader. Of course, if you're  using Nextflow this will generally be taken care of by the appropriate config file and should be transparent.

### kraken2

Download the latest image

`singularity pull docker://quay.io/h3abionet_org/kraken2 `

This will create an image `kraken2.img` which contains the kraken2 suite plus auxiliary programs like dustmasker

Note that we do not have any databases inside the image to keep the image small. You need to download and build the databases. Here's an example: Assume that you have a directory `/local/kraken` and you're going to bulild the database inside that

```
singularity exec -B /local/kraken/:/mnt kraken2.simg kraken2-build --standard --threads 8 --db /mnt/krakdb
```
This binds the directory `/local/kraken` on the host to the `/mnt` directory in the singularity image. The directory `/mnt` is passed to the `kraken2-build` program to use for the data and the database will be called `krakdb`.
