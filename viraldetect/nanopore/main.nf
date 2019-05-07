#!/usr/bin/env nextflow
// Create a nextflow channel from fastq files
samples = Channel.fromPath ("${params.in_dir}/*.fastq.gz")

samples.into { samples_1; samples_2 }

process runKraken {
    tag { "${sample}.runKraken" }
    label 'kraken'
    memory { 4.GB * task.attempt }
    cpus { 8 }
    module 'bioinf/kraken2'
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
kraken_hits.into { kraken_hits_1; kraken_hits_2 }

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

process getIdentity {
    tag { "${sample}.getIdentity" }
    label 'getidentity'
    memory { 4.GB * task.attempt }
    cpus { 8 }
    publishDir "${params.out_dir}/${sample}", mode: 'copy', overwrite: false

    input:
    set val(sample), file(khits) from kraken_hits_1
    file(mmap) from minimap2
    file(a2t) from file(params.nucl_gb_accession2taxid_file)
    file(t2n) from file(params.names_dmp_file)
   
    output:
    file "identity-stats.txt" into identity_stats

    script:
    """
    get-identity.py -k ${khits} -s ${mmap} -a ${a2t} -n ${t2n} > identity-stats.txt
    """
}

process runKrona {
    tag { "${sample}.runKrona" }
    label 'krona'
    memory { 4.GB * task.attempt }
    cpus { 1 }
    module 'bioinf/krona'
    publishDir "${params.out_dir}/${sample}", mode: 'copy', overwrite: false

    input:
    set val(sample), file(hits) from kraken_hits_2

    output:
    file "krona.htm" into krona_reports

    script:
    """
    ktImportTaxonomy -m 3 -s 0 -q 0 -t 5 -i ${hits} -tax ${params.krona_db} -o krona.htm
    """
}
