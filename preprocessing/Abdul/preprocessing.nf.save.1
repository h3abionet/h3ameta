params.input = "quality_control/*.fq.gz"

in_file = file(params.input)

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

	fastqc -t 4 $in_file -o $base

	"""
	 
}


process runMultiQc {
        input: file f from step1_ch.collect() //one input
        output:
        file("results") into step2_ch //many outputs
        publishDir "results/MultiQcOriginal"
"""
#!/usr/bin/env bash
module load python36
multiqc $f -o results
"""
}



process runTrimmomatic{

        //input: file in_file //one input
        output:
        file("results") into step3_ch
         publishDir "results/TrimmedData"

        script:
        """
         #!/bin/bash
         
         input="/home/abdulrahman/h3ameta/preprocessing/Abdul/quality_control"
         for i in $input/*_1.fq.gz; 
         do
         withpath="${i}"
         filename=${withpath##*/}
         base="${filename%*_*.fq.gz}"
         sample_name=`echo "${base}" | awk -F ".fastq.gz" '{print $1}'` 
         java -jar /home/abdulrahman/h3ameta/preprocessing/Abdul/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 6 -trimlog $output/"${base}".log $input/"${base}"_1.fq.gz $input/"${base}"_R2.fq.gz 

         done
        """
        //notice in the above:
        //the default interpreter is bash.  This is used when there is no shebang line.
        //the input file is referenced according to a variable defined in the input clause within this process
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
