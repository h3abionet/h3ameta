params.input = "quality_control/*.fq.gz" //this can be overridden by providing --in_file on the command line.

in_file = file(params.input) //create a channel that can be used as input to the first step

process runFastQCOriginal {
	input: file in_file //one input
	output: 
        set file("$base/*.zip"), file("$base/*.html") into step1_ch //many outputs
	publishDir "results/" //where should the results get linked to from the work folder?
        cpus 4
	script:
        base = in_file.baseName
	"""
	#!/usr/bin/env bash

        mkdir $base
	fastqc -t 4 $in_file -o $base
	
	"""
	//notice in the above:
	//any interpreter can be used (e.g. bash, python, ruby, perl), so long as it's specified in the shebang (#!) line.
	//variables can be used with a dollar sign ($).  When they are used not to indicate a variable, they must be escaped.
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
