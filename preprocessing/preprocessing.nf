params.input = "input.txt" //this can be overridden by providing --in_file on the command line.

in_file = file(params.input) //create a channel that can be used as input to the first step

process step1{
	input: file in_file //one input
	output: file "output1_*.txt" into step1_ch //many outputs
	publishDir "results/" //where should the results get linked to from the work folder?

	script: //what do?!
	"""
	#!/usr/bin/env bash

	for OUTNUM in `seq 1 10`; do
		(cat $in_file; echo 'this is the output from step 1') > output1_\$OUTNUM.txt
	done
	"""
	//notice in the above:
	//any interpreter can be used (e.g. bash, python, ruby, perl), so long as it's specified in the shebang (#!) line.
	//variables can be used with a dollar sign ($).  When they are used not to indicate a variable, they must be escaped.

}

process step2{
	input: file f from step1_ch //one input
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
