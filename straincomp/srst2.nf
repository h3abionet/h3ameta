#!/usr/bin/env nextflow

/*
srst2-Nextflow 
Written by Penistacia, Heyam and Mushal
*/

//The parameters below can all be overridden with --parametername on the commandline (e.g. --in or --dataset_table)
params.in = ""
params.ResFinder.db = "/home/mehabo/h3ameta/straincomp/srst2/data"

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
	srst2 --input_pe $d --output ${d.baseName}_srst2_resistant.txt --log --gene_db $ResFinder
	"""
}




