// channel to paired end reads
Channel
    .fromFilePairs( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .into {data; data2}

// channel for kraken DB
krakenDB = file(params.krakenDB)

// channel for host genome
host_genome = file(params.genome)

// channel for host genome
viral_genomes = file(params.viralGenomes)

process runBowtie2{
    publishDir "output", mode: "copy", overwrite: false

    input:
    set val(pairId), file(in_fastq) from data2

    output:
    file("${pairId}_bowtie2.sam") into bowtie2_aligned

    """
    bowtie2 -x $viral_genomes -1  ${in_fastq.get(0)} -2 ${in_fastq.get(1)} -S ${pairId}_bowtie2.sam
    """
}

process getMappingstats {
    tag { "${aligned}.runSAM_STATS" }
    label 'mappingStats'
    memory { 8.GB * task.attempt }
    cpus { 20 }

    publishDir "output", mode: "copy", overwrite: false

    input:
    file aligned from bowtie2_aligned

    output:
    file "mappingStats.tsv" into mappingStats


    script:
    """
    SAM_STATS $aligned >> mappingStats.tsv
    """
}


process alignReadstoHostgenome {
     publishDir "output", mode: "copy", overwrite: false

    input: 
	set val(pairId), file(in_fastq) from data

	output:
	file "${pairId}.sam" into aligned

	script:
	"""
	bowtie2 -x $host_genome -1 ${in_fastq.get(0)} -2 ${in_fastq.get(1)}  -S ${pairId}.sam
	"""
}


process removeHostReads {
	publishDir "output", mode: "copy", overwrite: false
	
	input:
	file fq_align from aligned

	output:
	file "${fq_align.simpleName}_?.fastq" into clean_fq
	
	script:
	"""
       /usr/local/anaconda/envs/shared_env/bin/samtools view -bS $fq_align -o tmp.bam
        /usr/local/anaconda/envs/shared_env/bin/samtools view -b tmp.bam -f 12  -o unmapped.bam
	bedtools bamtofastq -i unmapped.bam -fq ${fq_align.simpleName}_1.fastq -fq2 ${fq_align.simpleName}_2.fastq
	"""
}

clean_fq.into {clean_fq_kraken}


process runKraken {
   	publishDir "output", mode: "copy", overwrite: false

	input: 
	file clean_fq  from clean_fq_kraken
    
	output:
	set val(clean_fq), file("report.kraken.tsv") into kraken_hits
        file "kraken.out"

    	"""
	kraken2 --memory-mapping --quick --db ${krakenDB} --threads $task.cpus ${clean_fq}  --report report.kraken.tsv >> kraken.out
    	"""
}


process runKrona {
    label 'krona'
    memory { 4.GB * task.attempt }
    cpus { 1 }
    publishDir "output", mode: 'copy', overwrite: false

    input:
    file(hits) from kraken_hits

    output:
    file "krona.htm" into krona_reports

    script:
    """
    ktImportTaxonomy -m 3 -s 0 -q 0 -t 5 -i ${hits} -o krona.htm
    """
}
