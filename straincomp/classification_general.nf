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
params.in = "/home/mushalali/h3ameta/straincomp/nature_fastq/*.gz"
params.db = "/local/kraken/mini8GB"
params.readlen = 150
params.tax_level = 'S'
params.dataset_table = '/home/mushalali/h3ameta/examples/test_data/datasets.tsv'
params.ResFinder = '/home/mushalali/h3ameta/straincomp/ResFinder_db.fasta'
params.ARGannot = '/home/mushalali/h3ameta/straincomp/ARGannot_r3_db.fasta'
data = file(params.in)
kraken_db = file(params.db)
dataset_table = file(params.dataset_table)
ResFinder_db = file(params.ResFinder)
ARGannot_db = file(params.ARGannot)
process kraken {
	//publishDir 'general_outs/', mode: 'copy', overwrite: true
		//don't publish to the outs folder, as this is an intermediate
		//publishing the results can be done by symlinking or copying the outputs
	input:
		file sample_fq from data //input channel is a file, as declared above
		file kraken_db
	output: file "${sample_fq.baseName}_kraken.tsv" into kraken_out //output channel consists of *kraken.tsv files

	//resource requirements are specified in this way:
	cpus 2
	time '4h'
	memory '10GB' //20

	script:
	"""
	#!/usr/bin/env bash
	#note that the code below can be written in any language, so long as the interpreter is named
	#in the above shebang line.
	#also note: variables from this or higher scopes can be named inside strings with a dollar sign, as below.
	#when part of a larger name, as in the tsv output filename below, the variable must be enclosed with {}
	# the final $sample_fq names the input.
	kraken2 --db $kraken_db/ --threads $task.cpus \
	--report ${sample_fq.baseName}_kraken.tsv \
	--quick --memory-mapping \
	$sample_fq
	"""
}

process bracken {
	publishDir 'general_outs/', mode: 'copy', overwrite: true //the output from this process will be published to the outs folder
  input: file f from kraken_out //the input comes from the kraken process output channel above, and will be referenced with the varaible f
  output: file "${f.baseName}_bracken.tsv" into bracken_out //the output filenames depend on the input filenames.
		//this is done to avoid name collisions in the outs publication folder; naming collisions are not possible during
		//process execution due to nextflow's built-in encapsulation

	//resource requirements
	cpus 1
	time '1h'
	memory { 8.GB * task.attempt }

	script:
	"""
	#!/usr/bin/env bash
	bracken -d $params.db -i $f \
	-o ${f.baseName}_bracken.tsv -r $params.readlen -l $params.tax_level \
	"""
}
bracken_out.into{bracken_out1; bracken_out2}

process krona {
	publishDir 'general_outs/', mode: 'copy', overwrite: true
	input: file k from bracken_out2
	output: file "krona_${k}.html" into krona_out
	cpus 1
	time '1h'
	memory '1GB'

	"""
	ktImportTaxonomy -m 3 -s 0 -q 0 -t 5 -i ${k} -o krona_${k}.html \
	 -tax /global/taxonomy
	"""
}
