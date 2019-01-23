params.input = "quality_control/*.fq.gz" //this can be overridden by providing --in_file on the command line.

in_file = file(params.input) //create a channel that can be used as input to the first step

process runFastQCOriginal {
	input: file in_file //one input
	output: 
           set file("$base/*.zip"), file("$base/*.html") into step1_ch //many outputs
	publishDir "results/FastQCOriginal" //where should the results get linked to from the work folder?
        cpus 4
	script:
        base = in_file.baseName
	"""
	#!/usr/bin/env bash

        mkdir $base
	fastqc -t 4 $in_file -o FastQCOriginal/$base
	
	"""
	//notice in the above:
	//any interpreter can be used (e.g. bash, python, ruby, perl), so long as it's specified in the shebang (#!) line.
	//variables can be used with a dollar sign ($).  When they are used not to indicate a variable, they must be escaped.
}


/*
process Cutadapt{
	input: file in_file  //one input
	output: file "output2.txt" into step2_ch //one output
	publishDir "results/"

	script:
	"""
	(cat $f; echo 'this is the output from step 2') > output2.txt
	"""
	//notice in the above:
	//the default interpreter is bash.  This is used when there is no shebang line.
	//the input file is referenced according to a variable defined in the input clause within this process
}

process step3{
	input: file f from step2_ch.collect() //many inputs.  Note the .collect()
	output: file "output3.txt" into step3_ch //one output
	publishDir "results/"

	script:
	"""
	(cat $f; echo 'this is the output from step 3') > output3.txt
	"""
}
*/
