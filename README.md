# h3ameta
H3ABionNet Metagenomics Workflow




## Docker images

We're assuming you're using singularity -- if using Docker it'll be a little simpler, so it's left as an exercise for the reader.

### kraken2

Download the latest image 

`singularity pull docker://quay.io/h3abionet_org/kraken2 `

This will create an image `kraken2.img` which contains the kraken2 suite plus auxiliary programs like dustmasker

Note that we do not have any databases inside the image to keep the image small. You need to download and build the databases. Here's an example: Assume that you have a directory `/local/kraken` and you're going to bulild the database inside that

```
singularity exec -B /local/kraken/:/mnt kraken2.simg kraken2-build --standard --threads 8 --db /mnt/krakdb
```
This binds the directory `/local/kraken` on the host to the `/mnt` directory in the singularity image. The directory `/mnt` is passed to the `kraken2-build` program to use for the data and the database will be called `krakdb`.
