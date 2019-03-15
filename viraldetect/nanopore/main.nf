#!/usr/bin/env nextflow
// Create a nextflow channel from fastq files
samples = Channel.fromPath ("${params.in_dir}/*.fastq.gz")

samples.into { samples_1; samples_2 }


// Run Kraken
process runKraken {
    tag { "${sample}.runKraken" }
    label 'kraken'
    memory { 4.GB * task.attempt }
    cpus { 8 }
    publishDir "${params.out_dir}/${sample}", mode: 'copy', overwrite: false

    input:
    file(sample) from samples_1

    output:
    set val(sample), file("kraken-hits.tsv") into kraken_hits

    script:
    """
    kraken2 --memory-mapping --quick \
    --db ${params.kraken_db} \
    --threads ${task.cpus} \
    --report kraken-hits.tsv \
    ${sample}
    """
}

/*
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
*/

process runMinimap2 {
    tag { "${sample}.runMinimap2" }
    label 'minimap2'
    memory { 4.GB * task.attempt }
    cpus { 8 }
    publishDir "${params.out_dir}/${sample}", mode: 'copy', overwrite: false

    input:
    file(sample) from samples_2

    output:
    file "minimap2.sam" into minimap2

    script:
    """
    minimap2 -a  ${params.minimap2_db} -t ${task.cpus} ${sample} > minimap2.sam
    """
}

process runKrona {
    tag { "${sample}.runKrona" }
    label 'krona'
    memory { 4.GB * task.attempt }
    cpus { 1 }
    publishDir "${params.out_dir}/${sample}", mode: 'copy', overwrite: false

    input:
    set val(sample), file(hits) from kraken_hits

    output:
    file "krona.htm" into krona_reports

    script:
    """
    ktImportTaxonomy -m 3 -s 0 -q 0 -t 5 -i ${hits} -tax ${params.krona_db} -o krona.htm
    """
}

// Here we should pull in Kraken, mapping stats and Krona visualisation. 
process generateReport {

    input: file kro from krona_ch
	   file kra from kraken_classified
	   file statCh from mappingStats
	   file statEb from mappingStats
	   file statHp from mappingStats


    output: "report.html" into report

    script :
	"""
	python2.7 final_report.py kra statCh statEb statHP kro report.html
	"""
}
