#!/usr/bin/env bash
module load python36
multiqc sample2_1.fq_fastqc.zip sample2_1.fq_fastqc.html sample1_1.fq_fastqc.zip sample1_1.fq_fastqc.html sample2_2.fq_fastqc.zip sample2_2.fq_fastqc.html sample1_2.fq_fastqc.zip sample1_2.fq_fastqc.html -o results
