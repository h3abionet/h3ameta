#!/usr/bin/env nextflow

/*
Kraken-Nextflow - a NextFlow port of the Bhatt lab's Snakemake workflow for taxonomic short read classification using
Kraken2.  Written by Eli Moss.  Original workflow written by Ben Siranosian and Eli Moss, and can be found at
github.com/bhattlab/metagenomics_workflows

This workflow performs initial classification, then refinement with Bracken, table aggregation
and finally visualization with a custom bargraph script written in R, as well as Krona.

The accompanying nextflow.config file allows this workflow to be run on the Stanford SCG cluster with the
SLURM scheduler.

 ~/nextflow ../wits_workshop/nextflow/taxonomic_classification/taxonomic_classification.nf  --tax_level S -resume \
 -profile wits --in ../wits_workshop/nextflow/test_data/*.fq

*/

//The parameters below can all be overridden with --parametername on the commandline (e.g. --in or --dataset_table)
params.in = "test_data/*.f*q" //please note that asterisks must be escaped on the command line
params.db = "/labs/asbhatt/data/program_indices/kraken2/kraken_custom_oct2018/genbank_bacteria/"
params.readlen = 150
params.tax_level = 'S'
params.dataset_table = 'h3ameta/examples/test_data/datasets.tsv'

data = Channel.fromPath(params.in).flatten()
krakdb = file(params.db)
dataset_table = file(params.dataset_table)


//data.flatten().subscribe{println it}


process kraken {
	//publishDir 'outs/' //, mode: symlink, overwrite: true
		//don't publish to the outs folder, as this is an intermediate
		//publishing the results can be done by symlinking or copying the outputs
	input:
		file d from data.flatten() //input channel is a file, as declared above
		file krakdb
	output: file "${d.baseName}_kraken.tsv" into kraken_ch //output channel consists of *kraken.tsv files

	//resource requirements are specified in this way:
	cpus 1
	time '4h'
	memory '20GB' //20
	container = 'quay.io/biocontainers/kraken2:2.0.7_beta--pl526h6bb024c_1' //''shub://bhattlab/wits_workshop:classification' //

	script:
	"""
	#!/usr/bin/env bash

	#note that the code below can be written in any language, so long as the interpreter is named
	#in the above shebang line.

	#also note: variables from this or higher scopes can be named inside strings with a dollar sign, as below.
	#when part of a larger name, as in the tsv output filename below, the variable must be enclosed with {}
	# the final $d names the input.

	kraken2 --db $krakdb/ --threads $task.cpus \
	--report ${d.baseName}_kraken.tsv \
	--quick --memory-mapping \
	$d
	"""
}

process bracken {
	publishDir 'outs/' //the output from this process will be published to the outs folder
  input: file f from kraken_ch.flatten() //the input comes from the kraken process output channel above, and will be referenced with the varaible f
  output: file "${f.baseName}_bracken.tsv" into bracken_ch //the output filenames depend on the input filenames.
		//this is done to avoid name collisions in the outs publication folder; naming collisions are not possible during
		//process execution due to nextflow's built-in encapsulation

	//resource requirements
	cpus 1
	time '1h'
	memory { 8.GB * task.attempt }
	container = 'quay.io/biocontainers/bracken:2.2--py27h2d50403_1' //'shub://bhattlab/wits_workshop:classification' //

	script:
	"""
	#!/usr/bin/env bash

	bracken -d $params.db -i $f \
	-o ${f.baseName}_bracken.tsv -r $params.readlen -l $params.tax_level \
	"""
}

//the bracken output channel is used by multiple downstream processes (collect_results and krona) and so must be split into
//two channels (one per downstream process)
bracken_ch.into{bracken_ch1; bracken_ch2}

process collect_results {
	publishDir 'outs/'
	input:
		file data from dataset_table
		file f from bracken_ch1.collect() //this process runs on the outputs of all bracken instances, so .collect() is added
	output: file 'class_long.tsv' into collect_results_ch
	cpus 1
	time '1h'
	memory {2.GB * task.attempt } //dynamic resource allocation!

	//this script uses only python builtins, so no image is needed
	//container = 'shub://bhattlab/wits_workshop:classification'

	script:
	"""
	collate_results.py $params.tax_level $data class_long.tsv $f
	"""
}

process barplot {
	publishDir 'outs/'
	input: file f from collect_results_ch
	output: file 'barplot.pdf' into barplot_ch

	//a container is probably overkill for this, but several R libraries are used
	container = 'shub://bhattlab/wits_workshop:classification'

	script:
	"""
	composition_barplot.R $f barplot.pdf
	"""
}


process krona {
	publishDir 'outs/'
	input: file k from bracken_ch2
	output: file "krona_${k}.html" into krona_ch
	cpus 1
	time '1h'
	memory '1GB'
	container = 'shub://bhattlab/wits_workshop:classification'

	"""
	ktImportTaxonomy -m 3 -s 0 -q 0 -t 5 -i ${k} -o krona_${k}.html \
	 -tax \$(which ktImportTaxonomy | sed 's/\\/ktImportTaxonomy//g')/taxonomy
	"""
}
