Here goes the code for the viral outbreak monitoring. Code is still in developing phase and not ready yet

## Detection of viral sequences in Illumina reads

## Detection of viral sequences in Nanopore reads

To test pipeline as is know. E.g.

```
nextflow -log nextflow.log run -w /spaces/gerrit/projects/metagenomics-hackathon-2019/nextflow/work -c nextflow.nanopore.config main.nanapore.nf -profile wits
```

![workflow](https://raw.githubusercontent.com/h3abionet/h3ameta/master/viraldetect/main.nanapore.png "Workflow")



