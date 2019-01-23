#!/usr/bin/env nextflow
// Create a nextflow channel from fastq files
samples = Channel.fromPath ("${params.in_dir}/*.fastq.gz")


// Run Kraken
process runKraken {
    """
    tag { "${sample}.runKraken" }
    memory { 4.GB * task.attempt }
    cpus { 1 }
    publishDir "${params.out_dir}/${sample}", mode: 'copy', overwrite: false
    echo true

    input:
    file(sample) from samples
    
    output:
    // kraken standard ouput
    
    script:
    """
    // run Kraken
    """
}}

// Filter out reads that match
process filterHumanReads {

    input:
   // Kraken file

   output:
   file "clean.fastq" into clean_fq
   // Original Kraken output file witout human reads

   script:
   """
   // Remove possible human reads from Kraken output
   """
}

process runMinimap2 {

    input:
    file fastq from clean_fq

    output:
    file "minimap2.bam" into bam

    script:
    """
    // Run minimap2
    """

}


// Here we should pull in Kraken, Mappingstats and Krona visualisation. Are we still running Krona
process generateReport {

    input:
    // Kraken output without human reads
    // minimap2 bam

    output:
    // Report.tsv

    script:
    """
    /// Need to run a custom script here that validates the minimap2 against the Kraken results and create a report with high confident hits.
    """

}
