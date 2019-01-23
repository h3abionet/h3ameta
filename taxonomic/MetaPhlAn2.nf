#!/usr/bin/env nextflow

/*
comments here about what this pipeline does*/

//The parameters below can all be overridden with --parametername on the commandline (e.g. --in or --dataset_table)
params.in = "/home/ansieyssel/h3ameta/test_datasets/taxonomic_classification/*.f*q*" /*here I should repalce it wit a generic thing so that its not "hard coded"
/* change this to our own parameters if needed: params.db = "/labs/asbhatt/data/program_indices/kraken2/kraken_unmod/standard/"
params.readlen = 150
params.tax_level = 'S'
params.dataset_table = 'test_data/datasets.tsv'*/
/* note to ansie, here i should define parameters that will be used my MetaHlAn, like params.outfile or pramas.format*/
sequencing_data = file(params.in)
//krakdb = file(params.db)
//dataset_table = file(params.dataset_table)

//this process details still have to be modified so that it can accept fasta or fq or fastq files
//perhaps interleaving the paired reads (merging?) with SeqTK before running this prcess is a good idea
//if we have paired end data....
process MetaPhlAn2 {
	//publishDir 'outs/' //, mode: symlink, overwrite: true
	//don't publish to the outs folder, as this is an intermediate
        //publishing the results can be done by symlinking or copying the outputs
	input:
		file d from sequencing_data //input channel is a file, as declared above
	output: 
		//file "${d}_MetaPhlAn2_profile.biom" into MetaPhlAn2_ch //output channel consists of *MetaPhlAn_prfile.txt files
		file "${d}_MetaPhlAn2_output.biom" into MetaPhlAn2_ch
	//resource requirements are specified in this way:
	cpus 2
	time '4h'
	memory '20GB' //20

	script:
	"""
	#!/usr/bin/env bash
	#metaphlan2.py $d --input_type fastq --nproc $task.cpus > ${d}_MetaPhlAn2_profile.txt
	metaphlan2.py $d --input_type fastq --nproc $task.cpus --o ${d}_MetaPhlAn2_profile.txt --biom ${d}_MetaPhlAn2_output.biom

	
	"""
}




/*process krona {
	publishDir 'outs/'
	input: file k from bracken_ch2
	output: file "krona_${k}.html" into krona_ch
	cpus 1
	time '1h'
	memory '1GB'


	"""
	ktImportTaxonomy -m 3 -s 0 -q 0 -t 5 -i ${k} -o krona_${k}.html \
	 -tax \$(which ktImportTaxonomy | sed 's/\\/ktImportTaxonomy//g')/taxonomy
	"""
}*/
