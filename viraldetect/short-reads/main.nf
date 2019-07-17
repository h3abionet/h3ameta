read_pairs = Channel.fromFilePairs(params.reads)
read_pairs.into {screen1; screen2 ; kraken ; report ; name; other}

other.subscribe{ println it }

println "viral detect pipeline    "
println "================================="

process removehost2 {
		tag { "{$sample_id}/removehost2" }
    memory { 4.GB * task.attempt }
		publishDir "${params.out_dir}/$sample_id/htmlreport", pattern: "*.html", mode: 'copy', overwrite: false
	  publishDir "${params.out_dir}/$sample_id/clean", pattern: "*.fastq", mode: 'copy', overwrite: false
		label 'fastq_screen'

		input:
	  set val(sample_id), file(reads) from screen2

    output :
		set val(sample_id), file ("${reads[0].simpleName}.100.tagged_filter.fastq.gz" )into (clean_for_pairfq1)
		set val(sample_id), file ("${reads[1].simpleName}.100.tagged_filter.fastq.gz") into (clean_for_pairfq2)

    script :
    	"""
 			fastq_screen ${reads[0]} --conf ${params.fastq_screen_conf} --nohits --aligner bowtie2 --outdir ./
   	  fastq_screen ${reads[1]} --conf ${params.fastq_screen_conf} --nohits --aligner bowtie2 --outdir ./
      """
}

process pairfq {
	tag { "{$sample_id}/pairfq" }
  memory { 4.GB * task.attempt }
	publishDir "${params.out_dir}"
	label 'pairfq'

	input:
  	set val (sample_id), file (file1) from clean_for_pairfq1
		set val (sample_id), file (file2) from clean_for_pairfq2

	output:
	set val(sample_id), file ("${file1.simpleName}.paired.fastq" )into (clean_for_kraken1, clean_for_bowtie1)
	set val(sample_id), file ("${file2.simpleName}.paired.fastq") into (clean_for_kraken2, clean_for_bowtie2)

  // Need to check fastq extension here if we need to gunzip
	script:
	  if ( file1.extension == "gz" && file2.extension == "gz" ) {
		"""
			gunzip -f ${file1}
			gunzip -f ${file2}
			pairfq_lite addinfo -i ${file1.baseName} -o ${file1.simpleName}.addinfo.fastq -p 1
			pairfq_lite addinfo -i ${file2.baseName} -o ${file2.simpleName}.addinfo.fastq -p 2
			pairfq_lite makepairs -f ${file1.simpleName}.addinfo.fastq -r ${file2.simpleName}.addinfo.fastq -fp ${file1.simpleName}.paired.fastq -rp ${file2.simpleName}.paired.fastq -fs ${file1.simpleName}.singleton.fastq -rs ${file2.simpleName}.singleton.fastq
		"""
		} else if ( file1.extension == "fastq" && file2.extension == "fastq" ) {
		"""
			pairfq_lite addinfo -i ${file1.baseName} -o ${file1.simpleName}.addinfo.fastq -p 1
			pairfq_lite addinfo -i ${file2.baseName} -o ${file2.simpleName}.addinfo.fastq -p 2
			pairfq_lite makepairs -f ${file1.simpleName}.addinfo.fastq -r ${file2.simpleName}.addinfo.fastq -fp ${file1.simpleName}.paired.fastq -rp ${file2.simpleName}.paired.fastq -fs ${file1.simpleName}.singleton.fastq -rs ${file2.simpleName}.singleton.fastq
		"""
		} else {
		"""
		echo "File type not supported"
		exit 1
		"""
		}
}

process bowtie2 {
	tag { "{$sample_id}/bowtie2" }
  memory { 4.GB * task.attempt }
	publishDir "${params.out_dir}"
	label 'bowtie2'

	input:
  	set val (sample_id), file (file1) from clean_for_bowtie1
		set val (sample_id), file (file2) from clean_for_bowtie2

	output:
		file "fq12_mapped.sam" into mapped

	script:
		"""
		bowtie2  -x ${params.viralGenomes} -1 $file1 -2 $file2 -S fq12_mapped.sam --no-unal
		"""
}

process getMappingstats1 {
  tag { "getMappingstats1" }
 	memory { 4.GB * task.attempt }
	publishDir "${params.out_dir}", mode: 'copy', overwrite: false
	label 'samtools'

 	input:
	file mapped12 from mapped

 	output:
  file "mappingStats12.txt" into  mappingstat_for_report12

  script:
	  """
	  samtools view -S -b -f 12 -F 256 ${mapped12} | samtools sort - > sample12.sorted.bam
		samtools index sample12.sorted.bam
		samtools idxstats sample12.sorted.bam > mappingStats12.txt
  	"""
}

process kraken {
    tag { "kraken" }
   	memory { 4.GB * task.attempt }
    publishDir "${params.out_dir}", mode: 'copy', overwrite: false
		label 'kraken'

    input:
			set val (samplekr1) ,file (file1) from clean_for_kraken1
			set val (samplekr2) ,file (file2) from clean_for_kraken2
		output:
			file("kraken-hits.tsv") into (kraken_hits , kraken_report)

		script:
			"""
		  kraken2 --memory-mapping --quick --db ${params.kraken_db} --threads 1 --report kraken-hits.tsv --paired $file1 $file2
	   	"""
}

process krona {
		tag { "krona" }
  	memory { 4.GB * task.attempt }
		publishDir "${params.out_dir}", mode: 'copy', overwrite: false
		label 'krona'

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
  tag { "final_report_short_reads" }
  memory { 4.GB * task.attempt }
	publishDir "${params.out_dir}", mode: 'copy', overwrite: false
	label 'final_report_short_reads'

	input:
		file (kro) from krona_result
		file (kra) from kraken_report
		file (stat) from mappingstat_for_report12

	output:
		file ("final-report.html") into final_rep

	script:
	"""
	final-report-short-reads.py  $kra $stat $kro final-report.html ${params.search_pattern}
	"""
}
