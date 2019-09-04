//nextflow mdr4.nf --reads 'input/*_{1,2}.fastq' -c illumina3.config 
params.outdir = "/home/hanachi/Hack_M"
params.screen = "/home/hanachi/shotgun3/shotg_software/fastq_screen_v0.11.1/fastq_screen"
params.reads = "$baseDir/input/*{1,2}.fastq"
params.bowtie2 ="/usr/bin/bowtie2"
params.kraken = "/usr/local/kraken2-master/kraken2"
params.viralGenomes ="/home/hanachi/Hack_M/REF/all_viral/all_"
params.krakenDB2="/home/hanachi/Hack_M/REF/kraken/minikraken2_v2_8GB_201904_UPDATE/"
params.kronadb="/home/hanachi/Hack_M/REF/krona/taxonomy/taxonomy.tab"
params.finalreport = "/home/hanachi/Hack_M/final_5.py"
params.fastavir = "/home/hanachi/Hack_M/REF/all_viral/all.fasta"


read_pairs = Channel.fromFilePairs(params.reads,  flat: true)



read_pairs.into {info; screen ; kraken ; report ; name}
//name.println() 

println "viral detect pipeline    "
println "************************************"

//println "Files in process: "
println "===================================="

println "===================================="

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

process removehost2 {
	tag {pair_id }
    	label 'Screening_for_host_DNA'
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
        label 'Screening_for_host_DNA'
        cpus {30}
	tag {a}
	publishDir "${params.outdir}/${a[1].simpleName}/paired", mode: 'copy'
        input :
        val a  from clean_for_fastqlite1
       

        output :
        file  '*.paired.fastq' into (clean_for_kraken1, clean_for_bowtie1)
	file '*singleton.fastq' into single
	
        script :

"""
/usr/local/pairfq_lite addinfo -i ${a[0]} -o ${a[0]}.addinfo.fastq -p 1
/usr/local/pairfq_lite addinfo -i ${a[1]} -o ${a[1]}.addinfo.fastq -p 2

/usr/local/pairfq_lite makepairs -f ${a[0]}.addinfo.fastq -r ${a[1]}.addinfo.fastq -fp ${a[0].simpleName}.1.paired.fastq -rp ${a[1].simpleName}.2.paired.fastq -fs ${a[0].simpleName}.1.singleton.fastq -rs ${a[1].simpleName}.2.singleton.fastq

echo ${a[0]}.addinfo.fastq
echo ${a[1]}.addinfo.fastq
echo ${a[1].simpleName}_2.paired.fastq
echo ${a[1]}.getSimpleName_2.paired.fastq 


"""
}

process bowtie {
 	label 'bow'
	cpus {30}
	echo true	
	publishDir "${params.outdir}/${b[1].simpleName}/bowtie", mode: 'copy'
	input:

val b from  clean_for_bowtie1
   	

 	output : 
	file '*.sam'  into mapped
	
	script :	

	"""
	bowtie2  -x ${params.viralGenomes} -1 ${b[0]} -2 ${b[1]} -S ${b[1].simpleName}.sam  --no-unal
	 """
}


process kraken {
	label 'kraken'
	cpus {20}
publishDir "${params.outdir}/${c[1].simpleName}/kraken", mode: 'copy'

input :
val c from clean_for_kraken1

output :

file '*.tsv' into (kraken_for_report, kraken_for_krona) 


script :
"""
echo ${c[0]}
echo ${c[1]}
 ${params.kraken} --memory-mapping --quick --db ${params.krakenDB2} --threads 1 --report  ${c[1].simpleName}.tsv  --paired ${c[0]}  ${c[1]}
	   	"""
}

process krona {
	label 'krona'
	cpus {20}
publishDir "${params.outdir}/${krakenfile.simpleName}/krona",  pattern:"*.html", mode: 'copy'
publishDir "${params.outdir}/${krakenfile.simpleName}/final_reporting",  pattern:"*.html", mode: 'copy'	
input : 
	set val(krakenfile) , file (krak) from kraken_for_krona

	output:
	file '*.html'  into krona_result

	script :
	"""
	ktImportTaxonomy -m 3 -s 0 -q 0 -t 5 -i $krak -o ${krakenfile.simpleName}.html
	"""
}


process stat {

	label 'bowtie_stat'
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
	tag { "${id_sample.baseName}" }
	publishDir "${params.outdir}/${id_sample.simpleName}/final_reporting",  mode: 'copy', pattern:"*.html", overwrite: false
	
	input:
	
	set val (id_sample), file (kro) from krona_result
	file (kra) from kraken_for_report
	file (stat) from mappingstat
	
	
	output:
	
	file '*html' into final_rep

	script:
	"""
	
	grep  ">" ${params.fastavir} > ${params.fastavir}_pattern.txt
	python2.7 ${params.finalreport} $kra $stat  ${kro} ${params.fastavir}_pattern.txt ${id_sample.simpleName}_report.html
	"""
}

