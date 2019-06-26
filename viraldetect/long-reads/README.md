## Detection of viral sequences in long reads (Nanopore or PacBio)
The workflow runs viral sequences through a process to remove human / "contaminant" reads. The cleaned set is run through Kraken2 and is also mapped against all available viral genomes using minimap2. A script parses the results, finds the highest hits and tries to match identities between the highest hits in the Kraken2 and minimap2 results. A final HTML output report is generated from the results generated during the workflow.

### Workflow diagram
![workflow](https://raw.githubusercontent.com/h3abionet/h3ameta/master/viraldetect/long-reads/main.png "Workflow")

## Software requirements
* [Nextflow](https://www.nextflow.io/)
* [Singularity](https://www.sylabs.io/guides/3.0/user-guide/installation.html)

### Reference datasets needed
1) Complete human genome (build 38) - `wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.38_GRCh38.p12/GCF_000001405.38_GRCh38.p12_genomic.fna.gz`
2) Complete database of viral genomes - `ftp://ftp.ncbi.nih.gov/genomes/Viruses/all.fna.tar.gz`
3) Kraken complete viral, bacterial and human DB "minikraken2_v2_8GB" - `wget ftp://ftp.ccb.jhu.edu/pub/data/kraken2_dbs/minikraken2_v2_8GB_201904_UPDATE.tgz`
4) NCBI accession id to taxonomy - `wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/gi_taxid_nucl.dmp.gz`
5) NCBI taxonomy to names - `wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz`

### Test datasets
#### Nanopore
We tested on small subsets of Nanopore reads downloaded from SRA. We did a search like this: https://www.ncbi.nlm.nih.gov/sra/?term=virus+nanopore+wgs, selected MinION or GridION and then selected viral sequences. Subsets of 100 sequences were created from the downloads below to use in testing.
1) https://www.ncbi.nlm.nih.gov/sra/ERX3280332[accn] - ERR3253560 - Hepatitis B virus. The test data can be downloaded [here](http://web.cbio.uct.ac.za/~gerrit/downloads/viraldetect/nanopore/hepatitis-b-virus/in/ERR3253560.100.fastq.gz)
2) https://www.ncbi.nlm.nih.gov/sra/SRX5772425[accn] - SRR8993485 - Camelpox virus Negev2016. The test data can be downloaded [here](http://web.cbio.uct.ac.za/~gerrit/downloads/viraldetect/nanopore/camelpox_virus_negev2016/in/SRR8993485.100.fastq.gz)
3) https://www.ncbi.nlm.nih.gov/sra/ERX1201052[accn] - ERR1121624 - Human alphaherpesvirus 1 strain 17. The test data can be downloaded [here](http://web.cbio.uct.ac.za/~gerrit/downloads/viraldetect/nanopore/human_alphaherpesvirus_1_strain_17/in/ERR1121624.100.fastq.gz)
### PacBio
We tested on one subset of PacBio reads downloaded from SRA. We did a search like this: https://www.ncbi.nlm.nih.gov/sra/?term=virus+pacbio+wgs, and then selected a set that contained viral sequences. A subsets of 100 sequences were created from the download below to use in testing.
1) https://www.ncbi.nlm.nih.gov/sra/SRX2829016[accn] - SRR5569827 - Human alphaherpesvirus 1 - HSV1 KOS-YA. The test data can be downloaded [here](http://web.cbio.uct.ac.za/~gerrit/downloads/viraldetect/pacbio/human_alphaherpesvirus_HSV1_KOS-YA/in/SRR5569827.100.fastq.gz)
### Running tests
The tests are run on a local machine using the `standard` profile. The `nextflow.config` also contains an entry for a SLURM setup as an example.

#### Nanopore

##### Hepatitis B virus
```
nextflow -log nextflow.log run -w /home/gerrit/scratch/metagenomics-hackathon-2019/nextflow/work -c nextflow.config main.nf -profile standard \
--in_dir /home/gerrit/scratch/metagenomics-hackathon-2019/datasets/test/nanopore/hepatitis-b-virus/100/
--out_dir /home/gerrit/scratch/metagenomics-hackathon-2019/nextflow/100/out \
--kraken_db /home/gerrit/scratch/metagenomics-hackathon-2019/dbs/kraken2/minikraken2_v2_8GB_201904_UPDATE \
--minimap2_decontaminant_db /home/gerrit/scratch/metagenomics-hackathon-2019/dbs/minimap2/human/GCF_000001405.38_GRCh38.p12_genomic.fna.gz \
--minimap2_viral_db /home/gerrit/scratch/metagenomics-hackathon-2019/dbs/minimap2/complete-viral/all.fasta \
--nucl_gb_accession2taxid_file /home/gerrit/scratch/metagenomics-hackathon-2019/dbs/ncbi-taxonomy/nucl_gb.accession2taxid  \
--names_dmp_file /home/gerrit/scratch/metagenomics-hackathon-2019/dbs/ncbi-taxonomy/names.dmp \
--read_type nanopore
```

Output can be accessed [here](http://web.cbio.uct.ac.za/~gerrit/downloads/viraldetect/nanopore/hepatitis-b-virus/out). Here we do get a match between Kraken2 and minimap2 hits.

##### Camelpox virus Negev2016
```
nextflow -log nextflow.log run -w /home/gerrit/scratch/metagenomics-hackathon-2019/nextflow/work -c nextflow.config main.nf -profile standard  \
--in_dir /home/gerrit/scratch/metagenomics-hackathon-2019/datasets/test/nanopore/camelpox_virus_negev2016/100/ \
--out_dir /home/gerrit/scratch/metagenomics-hackathon-2019/nextflow/camelpox_virus_negev2016/100/out \
--kraken_db  /home/gerrit/scratch/metagenomics-hackathon-2019/dbs/kraken2/minikraken2_v2_8GB_201904_UPDATE \
--minimap2_decontaminant_db /home/gerrit/scratch/metagenomics-hackathon-2019/dbs/minimap2/human/GCF_000001405.38_GRCh38.p12_genomic.fna.gz \
--minimap2_viral_db /home/gerrit/scratch/metagenomics-hackathon-2019/dbs/minimap2/complete-viral/all.fasta \
--nucl_gb_accession2taxid_file /home/gerrit/scratch/metagenomics-hackathon-2019/dbs/ncbi-taxonomy/nucl_gb.accession2taxid \
--names_dmp_file /home/gerrit/scratch/metagenomics-hackathon-2019/dbs/ncbi-taxonomy/names.dmp \
--read_type nanopore
```

Output can be accessed [here](http://web.cbio.uct.ac.za/~gerrit/downloads/viraldetect/nanopore/camelpox_virus_negev2016/out) . Though many human reads have initially been removed, Kraken2 still picks up human reads in the cleaned set. Camelpox virus are picked up as highest hit in minmap2 results but secondly ranked in Kraken2 results.

##### Human alphaherpesvirus 1 strain 17
```
nextflow -log nextflow.log run -w /home/gerrit/scratch/metagenomics-hackathon-2019/nextflow/work -c nextflow.config main.nf -profile standard \
--in_dir /home/gerrit/scratch/metagenomics-hackathon-2019/datasets/test/nanopore/human_alphaherpesvirus_1_strain_17/100 \
--out_dir /home/gerrit/scratch/metagenomics-hackathon-2019/nextflow/human_alphaherpesvirus_1_strain_17/100/out \
--kraken_db /home/gerrit/scratch/metagenomics-hackathon-2019/dbs/kraken2/minikraken2_v2_8GB_201904_UPDATE \
--minimap2_decontaminant_db /home/gerrit/scratch/metagenomics-hackathon-2019/dbs/minimap2/human/GCF_000001405.38_GRCh38.p12_genomic.fna.gz \
--minimap2_viral_db /home/gerrit/scratch/metagenomics-hackathon-2019/dbs/minimap2/complete-viral/all.fasta \
--nucl_gb_accession2taxid_file /home/gerrit/scratch/metagenomics-hackathon-2019/dbs/ncbi-taxonomy/nucl_gb.accession2taxid \
--names_dmp_file /home/gerrit/scratch/metagenomics-hackathon-2019/dbs/ncbi-taxonomy/names.dmp \--read_type nanopore
--read_type nanopore
```

Output can be accessed [here](http://web.cbio.uct.ac.za/~gerrit/downloads/viraldetect/nanopore/human_alphaherpesvirus_1_strain_17/out). Again many human reads have initially been removed however Kraken2 still picks up human reads in the cleaned set. Human alphaherpesvirus 1 is ranked first in the Kraken2 results. Escherichia phage hk630 is ranked first in the in the minimap2 results.#### pacbio
#### PacBio

##### Human alphaherpesvirus - HSV1 KOS-YA
```
nextflow -log nextflow.log run -w /home/gerrit/scratch/metagenomics-hackathon-2019/nextflow/work -c nextflow.config main.nf -profile standard \
--in_dir /home/gerrit/scratch/metagenomics-hackathon-2019/datasets/test/pacbio/human_alphaherpesvirus_HSV1_KOS-YA/100 \
--out_dir /home/gerrit/scratch/metagenomics-hackathon-2019/nextflow/human_alphaherpesvirus_HSV1_KOS-YA/100/out \
--kraken_db /home/gerrit/scratch/metagenomics-hackathon-2019/dbs/kraken2/minikraken2_v2_8GB_201904_UPDATE \
--minimap2_decontaminant_db /home/gerrit/scratch/metagenomics-hackathon-2019/dbs/minimap2/human/GCF_000001405.38_GRCh38.p12_genomic.fna.gz \
--minimap2_viral_db /home/gerrit/scratch/metagenomics-hackathon-2019/dbs/minimap2/complete-viral/all.fasta \
--nucl_gb_accession2taxid_file /home/gerrit/scratch/metagenomics-hackathon-2019/dbs/ncbi-taxonomy/nucl_gb.accession2taxid \
--names_dmp_file /home/gerrit/scratch/metagenomics-hackathon-2019/dbs/ncbi-taxonomy/names.dmp \--read_type nanopore
--read_type pacbio
```
Output can be accessed [here](http://web.cbio.uct.ac.za/~gerrit/downloads/viraldetect/pacbio/human_alphaherpesvirus_HSV1_KOS-YA/out). In this case Kraken2 picks up human reads in the cleaned set and ranks human first where minimap2 ranks human alphaherpesvirus 1 as the top hit.

### Run times
For a 100 read set all processes except runMinimap2Decontaminate and getIdentity runs under 10 seconds. getIdentity
does an exhaustive search and can probably be improved, however the search times will not increase as the dataset size grow. runMinimap2Decontaminate will scale in time as the dataset size grow. Minimap2 can be parallelised therefor speeds can be increased where more resources are available.
![run-times](https://raw.githubusercontent.com/h3abionet/h3ameta/master/viraldetect/long-reads/run-times.png "Run times")

### Output folder structure
* **sample_name**
  * `cleaned.fastq` - Contaminant (e.g. human) reads removed from Fastq file.
  * `minimapdc.sam` - Minimap2 SAM output from aligning original reads against contaminant read set (e.g. human)  
* **cleaned.fastq** All of the stats/reports below are on the decontaminated/cleaned read set
  * `identity-stats.txt`- Checks done for highest Kraken2 can minmap2 matches. All of this are on the decontaminated/cleaned read set.
  * `final-repoprt.html` - Final report. Currently rendering Kraken2 and Krona results.
  * `krona.html` - Krona HTML report generated from Kraken2 results.
  * `kraken.html` - Kraken2 HTML report generated from Kraken2 TSV report.
  * `cleaned.fastq` - Contaminant (e.g. human) reads removed from Fastq file. (Duplicate file)
  * `minimap2viral.sam` - Minimap2 SAM output from aligning cleaned reads against all viral genomes.

### Container setup
All Nextflow processes make use of Singularity containers to run. Singularity containers are converted from Docker containers store on quay.io . All of the containers are stored [here](https://quay.io/user/grbot) except for the the minimap2 container that is stored [here](https://quay.io/biocontainers).

## To do
* Add additional reports to `final-report.html`.
  * These are the minimap2 SAM results after aligning against the human genome / "contaminant" reads
  * These are the minimap2 SAM results after aligning against the all viral genomes
  * The identity check results from `identity-stats.txt`
* Revise human contaminant checking.
* Revise matching of Kraken2 and minimap2 results.
