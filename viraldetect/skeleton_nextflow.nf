
#!/usr/bin/env nexflow

# define all input files
params.fq1="data/cleanReads/*1.fq"	# input directory for fastq files
params.fq2="data/cleanReads/*2.fq"
params.krakenDB=“data/krakenDB”	# Path to kraken DB
params.genome=“data/HumanGenome” #path to human genome
params.bitmask=”data/HumanGenome_index” #bmtagger index file
params.bmtOut=”data/bmtagger_out”
params.sampleNames=		#Names (prefixes) for each sample (and/or fastq files)



# create a nextflow channel
input_fq1 = Channel.fromPath("${params.fq1}”)
input_fq2 = Channel.fromPath("${params.fq2}”)




# four nextflow processes for each step

# process 1: get processed reads and filter out 

process removeHostReads {
Input: 
	file input_fq1
	file input_fq2
Output: file params.bmtOut into kraken_ch

Script:

"""
# Run bmtagger
# Make index for bmfilter
bmtool -d ${params.genome} -o ${params.bitmask} -A 0 -w 18 

#Run BMTtagger
Bmtagger.sh -b ${params.bitmask} –x reference.srprism –T tmp -q0 - 1 ${input_fq1} -2 ${input_fq2} -o ${params.bmtOut}

"""
}
  
process runKraken {
            
             input: file ${params.bmtOut}

             file seqs from host_free_reads             output: 
             file “*.kraken” into kraken_classified

             script:
             “””
             Kraken  —db ${params.krakenDB} ${seqs} > ${}.kraken
             “””


}

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


process runSNAP {
            
             input: 


             output: 


             script:


}


process runMinimap2 {
            
             input: 


             output: 


             script:


}


process getMappingstats {
            
             input: 


             output: 


             script:


}

process generateReport {
            
             input: 


             output: 


             script:


}



