
# h3ameta - straincomp

## This pipeline is designed for detecting bacterial species in metagenomics data, counts of single-nucleotide variants (SNVs) between input samples, building phylogenetic tree to show strains relatedness and identification of antibiotic resistance genes.

![straincomp_flowchart](https://github.com/h3abionet/h3ameta/blob/master/straincomp/documentation/README/straincomp_flowchart.png "straincomp_flowchart.png")

# Softwares:

- [Nextflow](https://www.nextflow.io/)
- [Singularity](https://sylabs.io/guides/3.0/user-guide/installation.html)

# Dependencies:

- [Kraken](http://ccb.jhu.edu/software/kraken/)
- [Bracken](https://github.com/jenniferlu717/Bracken)
- [Krona](https://github.com/marbl/Krona/wiki) 
- [Relative abundance](https://github.com/h3abionet/h3ameta/blob/master/examples/taxonomic_classification/bin/composition_barplot.R)
- [MetaPhlAn2](http://huttenhower.sph.harvard.edu/metaphlan2) 
- [SRST2](https://github.com/katholt/srst2) + [ariba](https://github.com/sanger-pathogens/ariba) (Antimicrobial Resistance)
- [Strainsifter](https://github.com/tamburinif/StrainSifter) + A specified reference genome (Relatedness)

# Input reads format:

## supports short reads either single-ends or paired-ends in fastq format.

# Output files:

## a- [Taxonomic_classification]
 kraken report (TSV) file displaying the assigned taxonomic labels to each metagenomics read.
     ![kraken_report](https://github.com/h3abionet/h3ameta/blob/master/straincomp/documentation/README/ test_kraken_bracken.tsv.png "test_kraken_bracken.tsv.png")
   - an interactive krona pie chart showing taxonomic distribution of reads with their abundancy 
   -[krona_html](file:///home/mimzy/metagenomics/krona_test1_kraken_bracken.tsv.html).
   ![krona_chart](https://github.com/h3abionet/h3ameta/blob/master/straincomp/documentation/README/krona_chart.png "krona_chart.png ")
   ![class-long.tsv](https://github.com/h3abionet/h3ameta/blob/master/straincomp/documentation/README/class-long.tsv.png "class-long.tsv.png ")
   ![barplot_pdf](https://github.com/h3abionet/h3ameta/blob/master/straincomp/documentation/README/barplot_pdf.png "barplot_pdf.png ")

### metaphlan2: waiting the results

## b- [Antimicrobial Resistance]: waiting results 
### During this step SRST2_AMR will screen all genes in the samples against ARGannot_db to detect the presence of any resistance genes. This  will output two files; a detailed report containing gene resistance results and a summary report of samples versus detected genes. If there is no any detected gene, only a detailed report will be outputted.

## c- [Relatedness]: waiting results 
During this step metagenomics reads will be mapped to a certain bacterial reference genome, the analysis will output two files; a pdf file display the phylogenetic tree and a file showing single-nucleotide variants (SNVs) counts between each pair of isolates.










```bash
#slurm @wits_cluster
                                       
