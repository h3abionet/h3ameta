//nextflow main_ill_dynamic2.nf --reads 'input/*_{1,2}*.fastq' -c config.config


log.info "********************************************"
log.info " H3Bionet Viral Detect pipeline short reads"
log.info "*******************************************"
log.info "Reads        : ${params.reads}"
log.info "Current home : $HOME"
log.info "Current user : $USER"
log.info "Current path : $PWD"
log.info "Working dir  : $workDir"
log.info "Output dir   : ${params.outdir}"

log.info "====================================================================================="

read_pairs = Channel.fromFilePairs(params.reads,  flat: true)
read_pairs.into {info; screen ; kraken ; report ; name}
//name.println() 


process print_sample_info {
    tag { "${sample_id}" }
    echo true
        
    input:
    set val(sample_id), file(fastq_r1_file), file(fastq_r2_file) from info
    script:
    """
    printf "[sample_info] sample: ${sample_id}\tFastqR1: ${fastq_r1_file}\tFastqR2: ${fastq_r2_file}\n"
    """
}

process rename {
	echo true
	tag { "${sample_id}" }
	
	input: 	
	set val(sample_id), file(fastq_R1_file), file(fastq_R2_file) from name

	output:
	file  '*.{1,2}*.fastq' into renamed 

	script :
	"""
	cp ${fastq_R1_file} \$(echo ${fastq_R1_file} | sed 's/_/./')
	cp ${fastq_R2_file}  \$(echo ${fastq_R2_file} | sed 's/_/./')
	"""
}

process removehost {
	tag { "${read}" }
    	label 'Screening_host_DNA'
   	cpus { 30 }
	publishDir "${params.outdir}/${pair_id}/fastq_screen_results", mode: 'copy'
   	
	input:
	 val read  from  renamed

	output:
	file '*{1,2}*.tagged_filter.fastq' into clean_for_fastqlite1

	file '*.tagged.fastq' into discared
	file '*html' into html_fastqscreen
	file '*txt' into info_fastqscreen
	
    script :
    	"""
 	echo ${read[0]}
	echo ${read[1]} 
	${params.screen} ${read[0]} --nohits --aligner bowtie2  --outdir ./
	${params.screen} ${read[1]}  --nohits --aligner bowtie2  --outdir ./
      """
}

process pairing {
        label 'Pairing'
        cpus {30}
	tag {cleanfastq[1].simpleName}
	publishDir "${params.outdir}/${cleanfastq[1].simpleName}/paired", mode: 'copy'
        input :
        val cleanfastq from clean_for_fastqlite1
       

        output :
        file  '*.paired.fastq' into (clean_for_kraken, clean_for_bowtie)
	file '*singleton.fastq' into single
	
        script :

"""
${params.lite} addinfo -i ${cleanfastq[0]} -o ${cleanfastq[0]}.addinfo.fastq -p 1
${pairfq_lite} addinfo -i ${cleanfastq[1]} -o ${cleanfastq[1]}.addinfo.fastq -p 2
${params.lite} makepairs -f ${cleanfastq[0]}.addinfo.fastq -r ${cleanfastq[1]}.addinfo.fastq -fp ${cleanfastq[0].simpleName}.1.paired.fastq -rp ${cleanfastq[1].simpleName}.2.paired.fastq -fs ${cleanfastq[0].simpleName}.1.singleton.fastq -rs ${cleanfastq[1].simpleName}.2.singleton.fastq

echo ${cleanfastq[0]}.addinfo.fastq
echo ${cleanfastq[1]}.addinfo.fastq
echo ${cleanfastq[1].simpleName}_2.paired.fastq
echo ${cleanfastq[1]}.getSimpleName_2.paired.fastq 


"""
}

process bowtie {
 	label 'bowtie2_step'
	cpus {30}
	echo true
	tag {"${pairedfastqbow[1].simpleName}"}
	publishDir "${params.outdir}/${pairedfastqbow[1].simpleName}/bowtie", mode: 'copy'
	
	input:
	val pairedfastqbow from  clean_for_bowtie
   	

 	output : 
	file '*.sam'  into mapped
	
	script :	

	"""
	bowtie2  -x ${params.viralGenomes} -1 ${pairedfastqbow[0]} -2 ${pairedfastqbow[1]} -S ${pairedfastqbow[1].simpleName}.sam  --no-unal
	 """
}


process kraken {
	label 'kraken_step'
	cpus {20}
	tag {"${pairedfastqkrk[1].simpleName}"}
	publishDir "${params.outdir}/${pairedfastqkrk[1].simpleName}/kraken", mode: 'copy'

	input :
	val pairedfastqkrk from clean_for_kraken

	output :
	file '*.tsv' into (kraken_for_report, kraken_for_krona) 

	script :
	"""
	echo ${pairedfastqkrk[0]}
	echo ${pairedfastqkrk[1]}
 	${params.kraken} --memory-mapping --quick --db ${params.krakenDB2} --threads 1 --report  ${pairedfastqkrk[1].simpleName}.tsv  --paired ${pairedfastqkrk[0]}  ${pairedfastqkrk[1]}
	   	"""
}

process krona {
	label 'krona_step'
	cpus {20}
	tag {"${krakenfile.simpleName}"}
	publishDir "${params.outdir}/${krakenfile.simpleName}/krona",  pattern:"*.html", mode: 'copy'
	publishDir "${params.outdir}/${krakenfile.simpleName}/final_reporting",  pattern:"*.html", mode: 'copy'	
	input : 
	set val (krakenfile), file ('*tsv') from kraken_for_krona

	output:
	file '*.krona.html'  into krona_result

	script :
	"""
	echo
	ktImportTaxonomy -m 3 -s 0 -q 0 -t 5 -i $krakenfile -o ${krakenfile.simpleName}.krona.html
	"""
}


process stat {
	tag {"${samfile.simpleName}"}
	label 'samtools_step'
	cpus {20}
	publishDir "${params.outdir}/${samfile.simpleName}/SAM", mode: 'copy'
	input : 
	set val(samfile) , file (sam) from mapped
	
	output:
	file '*txt' into mappingstat

	script:
	"""
	samtools view -b -S -f 12 -F 256 $samfile > ${samfile}.bam
	samtools sort ${samfile}.bam  -o ${samfile}.sorted.bam 
	samtools index  ${samfile}.sorted.bam
	samtools idxstats ${samfile}.sorted.bam >> ${samfile.simpleName}.txt
	"""

}


process report{
    	label 'Generate_final_report'
   	memory { 4.GB * task.attempt }
    	cpus { 8 }
	tag {"${kronafile.simpleName}"}
	publishDir "${params.outdir}/${kronafile.simpleName}/final_reporting",  pattern:"*.html", overwrite: true, mode : 'move'
	
	input: 
	set val (kronafile), file ('*html') from krona_result
	file (kra) from kraken_for_report	
	file (stat) from mappingstat
	
	
	output:
	
	file '*html' into final_rep

	script:
	
	"""
	echo $kronafile	
	grep  ">"  ${params.fastavir}  | cut -d "|" -f4,5 >  pattern.txt
	python2.7 ${params.finalreport} $kra $stat $kronafile pattern.txt ${kronafile.simpleName}_finalreport.html
	"""
}

