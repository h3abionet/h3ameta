// define all input files
params.fq1="data/*1.fastq.gz"	// input directory for fastq files
params.fq2="data/*2.fastq.gz"
params.krakenDB="krakenDB"	 // Path to kraken DB
params.genome="hostDB/host" // path to human genome
params.viralGenomes="viralDB/viruses"
params.katDB="katDB"	// # will symlink to Gerrit's /spaces/gerrit ...

// create a nextflow channel
input_fq1 = Channel.fromPath("${params.fq1}")
input_fq2 = Channel.fromPath("${params.fq2}")

// channel for kraken DB
krakenDB = file(params.krakenDB)

// channel for KAT kmer DB
kat_db = params.katDB

// channel for host genome
host_genome = file(params.genome)

// channel for host genome
viral_genomes = file(params.viralGenomes)

/* process alignReadstoHostgenome {
    input: 
	file fq1 from input_fq1
	file fq2 from input_fq2  

	output:
	file "fq12.sam" into aligned

	script:
	"""
	bowtie2 -x $host_genome -1 $fq1 -2 $fq2 -S fq12.sam
	"""
} */

process removeHostReads {
	
	input:
	//file fq_align from aligned
	file fq1 from input_fq1
        file fq2 from input_fq2
	
	output:
	//file "clean_f1.fastq" into clean_fq1
	//file "clean_f2.fastq" into clean_fq2
	file "clean*1*.fq" into clean_fq1
	file "clean*2*.fq" into clean_fq2
	
	//clean_fq1.into {clean_fq1_kraken; clean_fq1_bowtie2}
        //clean_fq2.into {clean_fq2_kraken; clean_fq2_bowtie2}

	script:
	"""
	# samtools view -bS $fq_align | samtools view -b -f 12 -F 256 | samtools sort -n > unmapped.bam
	# bedtools bamtofastq -i unmapped.bam -fq clean_f1.fastq -fq2 clean_f2.fastq
	kat filter seq --output_prefix clean --invert --seq $fq1 --seq2 $fq2 $kat_db
	"""
}

clean_fq1.into {clean_fq1_kraken; clean_fq1_bowtie2}
clean_fq2.into {clean_fq2_kraken; clean_fq2_bowtie2}

process runBowtie2{
    publishDir "output", mode: "copy", overwrite: false

    input:
    file clean_fq1  from clean_fq1_bowtie2
    file clean_fq2  from clean_fq2_bowtie2

    output:
    file "bowtie2Out.sam" into bowtie2_aligned

    script:
    """
    bowtie2 -x $viral_genomes -1 ${clean_fq1} -2 ${clean_fq2} -S bowtie2Out.sam
    """
}

process getMappingstats {
    publishDir "output", mode: "copy", overwrite: false

    input:
    set val(sample), file(aligned) from bowtie2_aligned

    output:
    file "mappingStats.txt" into mappingStats


    script:
    """
    samtools view -S -b ${aligned} | samtools sort -n -o sample.sorted.bam
    samtools index sample.sorted.bam
    samtools idxstats sample.sorted.bam>>mappingStats.txt
    #samtools flagstat ${aligned}>>mappingStats.txt
    #samtools view $aligned | cut -f3 | sort | uniq -c | \
    #awk -v i=$aligned '{print "\t"\$1"\t"\$2}' >> mappingStats.txt
    """


}


process runKraken {
   	 publishDir "output", mode: "copy", overwrite: false

	input: 
	file fq1 from clean_fq1_kraken
	file fq2 from clean_fq2_kraken  
    
	output:
	file "report.kraken.tsv" into kraken_classified

	script:
    	"""
	kraken --db $krakenDB --threads $task.cpus --fastq-input --paired $fq1 $fq2 > kraken.output

	kraken-report --db $krakenDB kraken.output > report.kraken.tsv
    	"""
}
