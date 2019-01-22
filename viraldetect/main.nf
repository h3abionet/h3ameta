
#!/usr/bin/env nexflow

# define all input files
params.fq1="data/cleanReads/*1.fq"	# input directory for fastq files
params.fq2="data/cleanReads/*2.fq"
params.krakenDB=“data/krakenDB”	# Path to kraken DB
params.genome=“data/HumanGenome” #path to human genome
params.bitmask=”data/HumanGenome_index” #bmtagger index file
params.sampleNames=		#Names (prefixes) for each sample (and/or fastq files)

# create a nextflow channel
input_fq1 = Channel.fromPath("${params.fq1}”)
input_fq2 = Channel.fromPath("${params.fq2}”)

# Create a bm index of human genome if not exist. Add nextflow code to check if ${params.bitmask} ELSE
process buildIndex {

   script:
		"""
    bmtool -d ${params.genome} -o ${params.bitmask} -A 0 -w 18
		"""
}


# Identify matching reads
process indentifyHostReads {

    input:
	  file input_fq1
	  file input_fq2

    output:
		file "host.matches" into host_matches

    script:
		"""
		bmtagger.sh -b ${params.bitmask} –x reference.srprism –T tmp -q0 - 1 ${input_fq1} -2 ${input_fq2} -o host.matches
		"""
}

# Filter out reads that match
process removeHostReadsF1 {
    input:
		file matches from host_matches

		output:
		file "clean_f1.fastq" into clean_fq1

		script:
		"""
		seqtk subseq ${matches} > clean_f1.fastq
		"""
}

process removeHostReadsF2 {
    input:
		file matches from host_matches

		output:
		file "clean_f2.fastq" into clean_fq2

		script:
		"""
		seqtk subseq ${matches} > clean_f2.fastq
		"""
}

# Maybe need to check if F1 and F2 read ids matches / paired before continuing further. Some tools might have an issue.

# Need to properly collate clean_fq1 and clean_fq2 so that they are used correctly downstream.

# Running Kraken, need to know how we can use the input?
process runKraken {

    input: file ${params.bmtOut}

    file seqs from host_free_reads             output:
    file “*.kraken” into kraken_classified

    script:
    “””
    Kraken  —db ${params.krakenDB} ${seqs} > ${}.kraken
    “””
}

## Need to check if we have short reads then run this.
process runBwa{

    input:
    file seqs from host_free_reads

    output:
    file “*.sam” into bwa_aligned

    script:
    “””
    bwa mem ${params.bwa_ref} read1.fq read2.fq > aln-pe.sam
    “””
}

## Need to check if we have long reads then run this.
process runMinimap2 {

    input:


    output:


    script:

}

# Here we can look at https://github.com/IARCbioinfo/bametrics-nf and https://github.com/IARCbioinfo/mpileup-nf
process getMappingstats {

    input:


    output:


    script:


}

# Here we should pull in Kraken, Mappingstats and Krona visualisation. Are we still running Krona
process generateReport {

    input:


    output:


    script:

}
