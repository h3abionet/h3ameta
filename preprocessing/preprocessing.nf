params.input = "quality_control/*.fq.gz" //this can be overridden by providing --in_file on the command line.

in_file = file(params.input) //create a channel that can be used as input to the first step

process runFastQCOriginal {
	input: file in_file //one input
	output: 
<<<<<<< HEAD
           set file("$base/*.zip"), file("$base/*.html") into step1_ch //many outputs
	publishDir "results/FastQCOriginal" //where should the results get linked to from the work folder?
=======
        set file("$base/*.zip"), file("$base/*.html") into step1_ch //many outputs
	publishDir "results/" //where should the results get linked to from the work folder?
>>>>>>> 7fe7644b5e66dd0e20320521cc855db062f7639b
        cpus 4
	script:
        base = in_file.baseName
	"""
	#!/usr/bin/env bash

        mkdir $base
<<<<<<< HEAD
	fastqc -t 4 $in_file -o FastQCOriginal/$base
=======
	//fastqc -t 4 $in_file -o $base
>>>>>>> 7fe7644b5e66dd0e20320521cc855db062f7639b
	
	"""
	//notice in the above:
	//any interpreter can be used (e.g. bash, python, ruby, perl), so long as it's specified in the shebang (#!) line.
	//variables can be used with a dollar sign ($).  When they are used not to indicate a variable, they must be escaped.
<<<<<<< HEAD
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
=======
>>>>>>> 7fe7644b5e66dd0e20320521cc855db062f7639b
}

process runMultiQc {
	input: file f from step1_ch.collect() //one input
        output:
        file("results") into step2_ch //many outputs
	publishDir "results/"
"""	
#!/usr/bin/env bash

multiqc . -o results
"""
}
*/
