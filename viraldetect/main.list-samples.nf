#!/usr/bin/env nextflow

// Create a nextflow channel from fastq files
samples = Channel.fromPath ("${params.in_dir}/*.fastq.gz")

// List samples
process listSamples {
    tag { "${sample}.listSamples" }
    memory { 4.GB * task.attempt }
    cpus { 1 }
    publishDir "${params.out_dir}/${sample}", mode: 'copy', overwrite: false
    echo true

    input:
    file(sample) from samples
    
    script:
    """
    printf "sample: ${sample}\n"
    """
}
