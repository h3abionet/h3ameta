// define all input files
params.fq1="data/*1.fastq.gz"	// input directory for fastq files
params.fq2="data/*2.fastq.gz"
params.krakenDB="krakenDB"	 // Path to kraken DB
params.genome="hostDB/host" // path to human genome
params.viralGenomes="viralDB/viruses"

// create a nextflow channel
input_fq1 = Channel.fromPath("${params.fq1}")
input_fq2 = Channel.fromPath("${params.fq2}")

// channel for kraken DB
krakenDB = file(params.krakenDB)

// channel for host genome
host_genome = file(params.genome)

// channel for host genome
viral_genomes = file(params.viralGenomes)

process alignReadstoHostgenome {
    input: 
	file fq1 from input_fq1
	file fq2 from input_fq2  

	output:
	file "fq12.sam" into aligned

	script:
	"""
	bowtie2 -x $host_genome -1 $fq1 -2 $fq2 -S fq12.sam
	"""
}

process removeHostReads {
	
	input:
	file fq_align from aligned

	output:
	file "clean_f1.fastq" into clean_fq1
	file "clean_f2.fastq" into clean_fq2
	
	//clean_fq1.into {clean_fq1_kraken; clean_fq1_bowtie2}
        //clean_fq2.into {clean_fq2_kraken; clean_fq2_bowtie2}

	script:
	"""
	samtools view -bS $fq_align | samtools view -b -f 12 -F 256 | samtools sort -n > unmapped.bam
	bedtools bamtofastq -i unmapped.bam -fq clean_f1.fastq -fq2 clean_f2.fastq
	"""
}

clean_fq1.into {clean_fq1_kraken; clean_fq1_bowtie2}
clean_fq2.into {clean_fq2_kraken; clean_fq2_bowtie2}

process runBowtie2{
    publishDir "stats", mode: "copy", overwrite: false

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
    publishDir "stats", mode: "copy", overwrite: false

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

/*
process runKraken {

	input: 
	file fq1 from clean_fq1_kraken
	file fq2 from clean_fq2_kraken  
    
	output:
    file "report.kraken.tsv" into kraken_classified

	script:
    """
	kraken --db $krakenDB --threads $task.cpus --report report.kraken.tsv \
	--fastq-input --paired $fq1 $fq2
    """
}
*/

// Need to check if we have short reads then run this.
/*
process runBowtie2{

    input:
    file clean_fq1  from clean_fq1_bowtie2
    file clean_fq2  from clean_fq2_bowtie2

    output:
    file "bwaOut.sam" into bwa_aligned

    script:
    """
    bowtie2 -x $viral_genomes ${clean_fq1} ${clean_fq2} > bowtie2Out.sam
    """
}


// Here we can look at https://github.com/IARCbioinfo/bametrics-nf and https://github.com/IARCbioinfo/mpileup-nf
//Using samtools for now
process getMappingstats {

    input:
    file aligned from bwa_aligned

    output:
    file "mappingStats.txt" into mappingStats


    script:
    """
    samtools view $aligned | cut -f3 | sort | uniq -c | \
    awk -v i=$aligned '{print i"\t"$1"\t"$2}' >> mappingStats.txt
    """


}

*/
