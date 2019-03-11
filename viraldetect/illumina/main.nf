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
    publishDir "output/reference_mapping", mode: "copy", overwrite: false

    input:
    set pairId, file(in_fastq) from data2

    output:
    file("${pairId}.sam") into bowtie2_aligned

    """
    bowtie2 -x $viral_genomes -1 ${in_fastq.get(0)} -2 ${in_fastq.get(1)} -S ${pairId}.sam
    """
}

process getMappingstats {
    label 'mappingStats'
    memory { 8.GB * task.attempt }
    cpus { 20 }

    publishDir "output/reference_mapping/${aligned.simpleName}_stats", mode: "copy", overwrite: false

    input:
    file aligned from bowtie2_aligned

    output:
    file "${aligned.simpleName}.tsv" into mapStats_Ch

    script:
    """
    samtools view -S -b ${aligned} | samtools sort -o sample.sorted.bam
    samtools index sample.sorted.bam   
    samtools idxstats sample.sorted.bam > ${aligned.simpleName}.tsv
    """
}


process alignReadstoHostgenome {
     publishDir "output/hostfree", mode: "copy", overwrite: false

    input: 
	set pairId, file(in_fastq) from data

	output:
	file "${pairId}.sam" into aligned

	script:
	"""
	bowtie2 -x $host_genome -1 ${in_fastq.get(0)} -2 ${in_fastq.get(1)} -S ${pairId}.sam
	"""
}


process removeHostReads {
	publishDir "output/hostfree", mode: "copy", overwrite: false
	
	input:
	file fq_align from aligned

	output:
	file "${fq_align.simpleName}_?.fastq" into clean_fq_kraken
	
	script:
	"""
        samtools view -bS $fq_align -o tmp.bam
        samtools view -b tmp.bam -f 12  -o unmapped.bam
	bedtools bamtofastq -i unmapped.bam -fq ${fq_align.simpleName}_1.fastq -fq2 ${fq_align.simpleName}_2.fastq
	"""
}

process runKraken {
	label 'kraken'
   	publishDir "output/kraken", mode: "copy", overwrite: false

	input: 
	file clean_fq from clean_fq_kraken
    
	output:
	file("${clean_fq.get(0).simpleName}.tsv") into kraken_hits
        file "kraken.out"

    	"""
	kraken2 --memory-mapping --quick --db ${krakenDB} --threads $task.cpus ${clean_fq}  --report ${clean_fq.get(0).simpleName}.tsv > kraken.out
    	"""
}
kraken_hits.into {krak_hits_krona; krak_hits_report}
process runKrona {
    label 'krona'
    memory { 4.GB * task.attempt }
    cpus { 1 }
    publishDir "output/krona", mode: 'copy', overwrite: false

    input:
    file(hits) from krak_hits_krona

    output:
    file "${hits.simpleName}.htm" into krona_reports

    script:
    """
    ktImportTaxonomy -m 3 -s 0 -q 0 -t 5 -i ${hits} -o ${hits.simpleName}.htm
    """
}

process generateReport {
    publishDir "output/reports/${kro.simpleName}_report", mode: 'copy', overwrite: false
    
    input: 
    file kro from krona_reports
    file kra from krak_hits_report
    file statCh from mapStats_Ch

    output: "report.html"

    script :
    """
    #!/usr/bin python
    sample=${statCh.simpleName}
    python final_report.py $kra $statCh $kro ${kro.simpleName} report.html
    """
}

