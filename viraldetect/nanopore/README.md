## Detection of viral sequences in Nanopore reads
The pipeline runs viral sequences trough Kraken2 and also maps against all available viral genomes using minimap2. 'n Script parses the results, finds the highest hits and tries to match identities between the higest hits in the Kraken2 and minimap2 results.

Datasets that you will need
1) NCBI taxonomy to names - `wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz`
2) NCBI accession id to taxonomy - `wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/gi_taxid_nucl.dmp.gz`
3) Kraken DB - `kraken2-build --standard --db kraken2-full`
4) Krona DB - `/opt/exp_soft/bioinf/KronaTools/updateTaxonomy.sh .`
5) Complete database of viral genomes - `ftp://ftp.ncbi.nih.gov/genomes/Viruses/all.fna.tar.gz`

We tested on small subsets of Nanopore reads downloaded from SRA. We did a search like this: https://www.ncbi.nlm.nih.gov/sra/?term=virus+nanopore+wgs (then select MinION or GridION and then selected viral sequences)
1) https://www.ncbi.nlm.nih.gov/sra/SRX5772425[accn] - SRR8993485 - Camelpox virus Negev2016
2) https://www.ncbi.nlm.nih.gov/sra/ERX1201052[accn] - ERR1121624 - Human alphaherpesvirus 1 strain 17
3) https://www.ncbi.nlm.nih.gov/sra/ERX3280332[accn] - ERR3253560 - Hepatitis B virus

### Hepatitis B virus

```
nextflow -log nextflow.log run -w /spaces/gerrit/projects/metagenomics-hackathon-2019/nextflow/work \
-c nextflow.config main.nf -profile wits \
--in_dir /spaces/gerrit/projects/metagenomics-hackathon-2019/datasets/test/nanopore/hepatitis-b-virus/small/ \
--out_dir /spaces/gerrit/projects/metagenomics-hackathon-2019/nextflow/out \
--kraken_db /spaces/gerrit/projects/metagenomics-hackathon-2019/dbs/kraken2/full/ \
--minimap2_db /spaces/gerrit/projects/metagenomics-hackathon-2019/dbs/minimap2/complete-viral/all.fasta \
--krona_db /opt/exp_soft/bioinf/KronaTools/taxonomy \
--nucl_gb_accession2taxid_file /spaces/gerrit/projects/metagenomics-hackathon-2019/dbs/ncbi-taxonomy/nucl_gb.accession2taxid  \
--names_dmp_file /spaces/gerrit/projects/metagenomics-hackathon-2019/dbs/ncbi-taxonomy/names.dmp
```

```
cat identity-stats.txt
Kraken2 max hit is: 10407 (taxonomy id) with a 90.0% of overall hits
Minimap2 max hit is: gi|944326253|ref|NC_028129.1| (GenBank/Refseq accession id) with a 88.88888888888889% of overall hits
The match between kraken2 and minimap2 is hepatitis b virus
Kraken2 top name hit is hepatitis b virus
Minimap2 top name hit is woolly monkey hepatitis b virus
```
Here we do get a match between Kraken2 and minimap2 hits.

### Camelpox virus Negev2016

```
nextflow -log nextflow.log run -w /spaces/gerrit/projects/metagenomics-hackathon-2019/nextflow/work \
-c nextflow.config main.nf -profile wits \
--in_dir /spaces/gerrit/projects/metagenomics-hackathon-2019/datasets/test/nanopore/camelpox_virus_negev2016/small \
--out_dir /spaces/gerrit/projects/metagenomics-hackathon-2019/nextflow/out \
--kraken_db /spaces/gerrit/projects/metagenomics-hackathon-2019/dbs/kraken2/full/ \
--minimap2_db /spaces/gerrit/projects/metagenomics-hackathon-2019/dbs/minimap2/complete-viral/all.fasta \
--krona_db /opt/exp_soft/bioinf/KronaTools/taxonomy \
--nucl_gb_accession2taxid_file /spaces/gerrit/projects/metagenomics-hackathon-2019/dbs/ncbi-taxonomy/nucl_gb.accession2taxid  \
--names_dmp_file /spaces/gerrit/projects/metagenomics-hackathon-2019/dbs/ncbi-taxonomy/names.dmp
```

```
cat identity-stats.txt
Kraken2 max hit is: 9606 (taxonomy id) with a 80.0% of overall hits
Minimap2 max hit is: gi|18640237|ref|NC_003391.1| (GenBank/Refseq accession id) with a 60.0% of overall hits
Kraken2 top name hit is homo sapiens
Minimap2 top name hit is camelpox virus
```
Humans are picked up as the majority hit in the Kraken2 results.

### Human alphaherpesvirus 1 strain 17

```
nextflow -log nextflow.log run -w /spaces/gerrit/projects/metagenomics-hackathon-2019/nextflow/work \
-c nextflow.config main.nf -profile wits \
--in_dir /spaces/gerrit/projects/metagenomics-hackathon-2019/datasets/test/nanopore/human_alphaherpesvirus_1_strain_17/small \
--out_dir /spaces/gerrit/projects/metagenomics-hackathon-2019/nextflow/out \
--kraken_db /spaces/gerrit/projects/metagenomics-hackathon-2019/dbs/kraken2/full/ \
--minimap2_db /spaces/gerrit/projects/metagenomics-hackathon-2019/dbs/minimap2/complete-viral/all.fasta \
--krona_db /opt/exp_soft/bioinf/KronaTools/taxonomy \
--nucl_gb_accession2taxid_file /spaces/gerrit/projects/metagenomics-hackathon-2019/dbs/ncbi-taxonomy/nucl_gb.accession2taxid  \
--names_dmp_file /spaces/gerrit/projects/metagenomics-hackathon-2019/dbs/ncbi-taxonomy/names.dmp
```

```
cat identity-stats
Kraken2 max hit is: 9606 (taxonomy id) with a 85.71428571428571% of overall hits
Minimap2 max hit is: gi|820945227|ref|NC_001806.2| (GenBank/Refseq accession id) with a 25.0% of overall hits
Kraken2 top name hit is homo sapiens
Minimap2 top name hit is human alphaherpesvirus 1
```
Humans are picked up as the majority hit in the Kraken2 results.

To do
1) Filter out human reads or build a Kraken2 database only on viral genomes.
2) Sort out and include the report generator into the pipeline.

### Workflow diagram
![workflow](https://raw.githubusercontent.com/h3abionet/h3ameta/master/viraldetect/nanopore/main.png "Workflow")

### Run times
All processes run ~10s. The getIdentity step takes about 4 mins. This is an exhaustive search and can probably be improved. getIdentity search times will however not increase as the datset size increases.  
![run-times](https://raw.githubusercontent.com/h3abionet/h3ameta/master/viraldetect/nanopore/run-times.png "Run times")

### On the current setup
Currently the pipeline runs on the Wits cluster. All steps run from Singularity containers, except Kraken2 and Krona which uses modules. The reason why we could not get singularity containers for this is:

1) Nextflow times out on the Kraken2 build
2) We could not find a correct Krona database that suited the Krona container.
