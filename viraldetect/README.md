Here goes the code for the viral outbreak monitoring. Code is still in developing phase and not ready yet.

## Development plan
a. Processes

b. Implementation assignments
- Voke to work on setting up custom kraken database
- Gerrit to set up the minimap2 database
- Stanford to work on the bwa alignment process for the Illumina process
- Alfred to work on the Krona process
- Mariem to work on the reports process

c. Timeline
- We should give each other updates weekly

d. Communication plan
- we will be communicate via our metagenomics-stream-3 slack channel


## Detection of viral sequences in Illumina reads
## Detection of viral sequences in Nanopore reads

To test pipeline as is know. E.g.

```
nextflow -log nextflow.log run -w /spaces/gerrit/projects/metagenomics-hackathon-2019/nextflow/work -c nextflow.nanopore.config main.nanapore.nf -profile wits
```

![workflow](https://raw.githubusercontent.com/h3abionet/h3ameta/master/viraldetect/main.nanapore.png "Workflow")

