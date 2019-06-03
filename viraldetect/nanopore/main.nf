#!/usr/bin/env nextflow
// Create a nextflow channel from fastq files
samples = Channel.fromPath ("${params.in_dir}/*.fastq.gz")

samples.into { samples_1; samples_2 }

process runMinimap2Decontaminate {
    tag { "${sample}.runMinimap2Decontaminate" }
    label 'minimap2'
    memory { 4.GB * task.attempt }
    cpus { 2 }
    publishDir "${params.out_dir}/${sample}", mode: 'copy', overwrite: false

    input:
    file(sample) from samples_1

    output:
    file "minimap2dc.sam" into minimap2dc

    script:
    """
    minimap2 -a  ${params.minimap2_decontaminant_db} -t ${task.cpus} ${sample} > minimap2dc.sam
    """
}

process removeReads {
    tag { "${sample}.removeReads" }
    label 'removereads'
    memory { 4.GB * task.attempt }
    cpus { 1 }
    publishDir "${params.out_dir}/${sample}", mode: 'copy', overwrite: false

    input:
    file(sample) from samples_2
    file(mmap) from minimap2dc

    output:
    file "cleaned.fastq" into cleaned_samples
    script:
    """
    remove-reads.py -s ${mmap} -i ${sample} > cleaned.fastq
    """
}

cleaned_samples.into { cleaned_samples_1; cleaned_samples_2 }

process runKraken {
    tag { "${sample}.runKraken" }
    label 'kraken'
    memory { 4.GB * task.attempt }
    cpus { 2 }
    publishDir "${params.out_dir}/${sample}", mode: 'copy', overwrite: false

    input:
    file(sample) from cleaned_samples_1

    output:
    set file(sample), file("kraken-hits.tsv") into kraken_hits

    script:
    """
    kraken2 --memory-mapping --quick \
    --db ${params.kraken_db} \
    --threads ${task.cpus} \
    --report kraken-hits.tsv \
    ${sample}
    """
}

kraken_hits.into { kraken_hits_1; kraken_hits_2; kraken_hits_3 }

process runMinimap2ViralCheck {
    tag { "${sample}.runMinimap2ViralCheck" }
    label 'minimap2'
    memory { 4.GB * task.attempt }
    cpus { 2 }
    publishDir "${params.out_dir}/${sample}", mode: 'copy', overwrite: false

    input:
    file(sample) from cleaned_samples_2

    output:
    file "minimap2viral.sam" into minimap2

    script:
    """
    minimap2 -a  ${params.minimap2_viral_db} -t ${task.cpus} ${sample} > minimap2viral.sam
    """
}

process getIdentity {
    tag { "${sample}.getIdentity" }
    label 'getidentity'
    memory { 4.GB * task.attempt }
    cpus { 1 }
    publishDir "${params.out_dir}/${sample}", mode: 'copy', overwrite: false

    input:
    set file(sample), file(khits) from kraken_hits_1
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
    publishDir "${params.out_dir}/${sample}", mode: 'copy', overwrite: false

    input:
    set file(sample), file(hits) from kraken_hits_2

    output:
    file "krona.html" into krona_report

    script:
    """
    ktImportTaxonomy -m 3 -s 0 -q 0 -t 5 -i ${hits} -o krona.html
    """
}

process runFinalReport {
    tag { "${sample}.runReport" }
    label 'finalreport'
    memory { 4.GB * task.attempt }
    cpus { 1 }
    publishDir "${params.out_dir}/${sample}", mode: 'copy', overwrite: false

    input:
    set file(sample), file(hits) from kraken_hits_3
    file (krona) from krona_report

    output:
    set file("final-report.html"), file("kraken.html") into final_report

    script:
    """
    final-report.py ${hits} ${krona} final-report.html
    """
}
