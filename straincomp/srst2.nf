#!/usr/bin/env nextflow

/*
srst2-Nextflow 
Written by Penistacia, Heyam and Mushal
*/

//The parameters below can all be overridden with --parametername on the commandline (e.g. --in or --dataset_table)
params.in = "/home/mushalali/mehabo-h3ameta/test_datasets/strain_comparison/test"
params.ResFinder.db = "/home/mushalali/mehabo-h3ameta/straincomp/data/ResFinder.fasta"

data = file(params.in)
ResFinder = file(params.ResFinder.db)


process srst2 {
	input:
		file d from data //input channel is a file, as declared above
		file ResFinder
	output: file "${d.baseName}_srst2_resistant.txt"


	script:
	"""
	#!/usr/bin/env bash
	srst2 --forward R1 --reverse R2 --output ${d.baseName}_srst2_resistant.txt --log --gene_db $ResFinder
	"""
}




