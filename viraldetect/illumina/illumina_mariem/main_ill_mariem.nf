read_pairs = Channel.fromFilePairs(params.reads,  flat: true)
read_pairs.into {screen1; screen2 ; kraken ; report ; name} 

println "viral detect pipeline    "
println "================================="


process removehost2 {
	tag { sample_id }
    	label 'Screening_for_host_DNA'
   	memory { 4.GB * task.attempt }
    	cpus { 30 }
        publishDir "${params.out_dir}/$sample_id/htmlreport", pattern: "*.html", mode: 'copy', overwrite: false
	publishDir "${params.out_dir}/$sample_id/clean", pattern: "*.fastq", mode: 'copy', overwrite: false
        input:
	set val(sample_id), file(read1), file(read2) from  screen2

        output :
        file "*" into (clean_for_kraken1,clean_for_kraken2, clean_for_bowtie1,clean_for_bowtie2 )


        script :
        """
	perl /home/mariemhanachi/fastq_screen_v0.13.0/fastq_screen $read2 --nohits --aligner bowtie2  --outdir ./
	perl  /home/mariemhanachi/fastq_screen_v0.13.0/fastq_screen $read1 --nohits --aligner bowtie2  --outdir ./
        """
}

process bowtie2 {
    	label 'Alignement_against_viral_genomes'
   	memory { 4.GB * task.attempt }
    	cpus { 30 }
	publishDir "${params.out_dir}" 

	input:
    	set val (file2) ,file ('*2.tagged_filter.fastq') from  clean_for_bowtie1
	set val (file1) , file ('*1.tagged_filter.fastq') from  clean_for_bowtie2
	

	output:
	file "fq12_mapped.sam" into mapped

	script:
	"""
	bowtie2  -x ${params.viralGenomes} -1 $file1 -2 $file2 -S fq12_mapped.sam --no-unal
	"""
}


process getMappingstats1 {
    	label 'Generate_alignement_statistcis'
   	memory { 4.GB * task.attempt }
    	cpus { 8 }
   	publishDir "${params.out_dir}", mode: 'copy', overwrite: false

	
 	input:
	file mapped12 from mapped

   	output:
    	 file "mappingStats12.txt" into  mappingstat_for_report12


    	script:
    	"""
       samtools view -b -S -f 12 -F 256 $mapped12 > aligned12.bam
       samtools sort aligned12.bam  -o sample12.sorted.bam 
       samtools index  sample12.sorted.bam
       samtools idxstats sample12.sorted.bam >> mappingStats12.txt
      """
}


process kraken {
    	label 'Taxonomic_classification_by_kraken'
   	memory { 4.GB * task.attempt }
    	cpus { 8 }
        publishDir "${params.out_dir}", mode: 'copy', overwrite: false

        input:
	set val (samplekr1) ,file ('*1.tagged_filter.fastq')    from  clean_for_kraken1
	set val (samplekr2) ,file ('*2.tagged_filter.fastq')     from  clean_for_kraken2
	output:
	
	file("kraken-hits.tsv") into (kraken_hits , kraken_report)
	
	script:
	"""
    	/usr/local/kraken2-master/kraken2 --memory-mapping --quick --db ${params.krakenDB2} --threads 5 --report kraken-hits.tsv --paired $samplekr1 $samplekr2
   	 """
}

process krona {
    	label 'Krona'
   	memory { 4.GB * task.attempt }
    	cpus { 8 }
	publishDir "${params.out_dir}", mode: 'copy', overwrite: false

	input: 
	file (krak) from kraken_hits

	output:
	file ("krona.html")  into krona_result

	script :
	"""
	ktImportTaxonomy -m 3 -s 0 -q 0 -t 5 -i $krak -o krona.html
	"""
}

process report{
    	label 'Generate_final_report'
   	memory { 4.GB * task.attempt }
    	cpus { 8 }
	publishDir "${params.out_dir}", mode: 'copy', overwrite: false
	
	input:
	
	file (kro) from krona_result
	file (kra) from kraken_report
	file (stat) from mappingstat_for_report12
	set pair_id, file(read1: "*1") from  report
	
	
	output:
	file ("${read1}_final_report.html") into final_rep

	script:
	"""
	python2.7 $baseDir/add_script/final_report_v5.py  $kra $stat $kro $read1"_final_report.html"
	"""
}

	

