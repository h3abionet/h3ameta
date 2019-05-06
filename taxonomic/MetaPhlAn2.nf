#!/usr/bin/env nextflow

/*
comments here about what this pipeline does*/
//The parameters below can all be overridden with --parametername on the commandline (e.g. --in or --dataset_table)
params.in = "/home/ansieyssel/h3ameta/test_datasets/taxonomic_classification/*.f*q*" /*here I should repalce it wit a generic thing so that its not "hard coded"
/* change this to our own parameters if needed: params.db = "/labs/asbhatt/data/program_indices/kraken2/kraken_unmod/standard/"
/* note to ansie, here i should define parameters that will be used my MetaHlAn, like params.outfile or pramas.format*/
data = file(params.in)
sequencing_data = file(params.in)

//this process details still have to be modified so that it can accept fasta or fq or fastq files
//perhaps interleaving the paired reads (merging?) 
/* Add a process for unzipping sipped files if needed.... */



process MetaPhlAn2 {

	input:
	//file d from sequencing_data //input channel is a file, as declared above
	file (infile) from sequencing_data// input channel is a file as declared above
	output: 
	file "${infile.baseName}_MetaPhlAn2_profile.txt" into MetaPhlAn2_ch
	//file "${infile.baseName}_MetaPhlAn2_output.biom" into AlphaDiversity_ch
	file "${infile.baseName}_MetaPhlAn2_microbes_list.tsv" into FunctionalProfiling_ch
	file "${infile.baseName}_MetaPhlAn2_bt2out.txt" into Bowtie2Output_ch
	//resource requirements are specified in this way:
	cpus 2
	time '4h'
	memory '4GB' //4

	
	script:
	"""	
	#!/usr/bin/env bash
	metaphlan2.py --input_type fastq --bowtie2out ${infile.baseName}_bt2out.txt --nproc ${task.cpus} $infile > ${infile.baseName}_MetaPhlAn2_microbes_list.tsv
	"""
}



