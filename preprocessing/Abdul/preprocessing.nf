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

/*

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
         filename=${withpath##/}
         base="${filename%*_*.fq.gz}"
         sample_name=`echo "${base}" | awk -F ".fastq.gz" '{print $1}'` 
         java -jar /home/abdulrahman/h3ameta/preprocessing/Abdul/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 6 -trimlog $output/"${base}".log $input/"${base}"_1.fq.gz $input/"${base}"_2.fq.gz $output/"${base}"_R1.trimmed_PE.fastq $output/"${base}"_R1.trimmed_SE.fastq $output/"${base}"_2.trimmed_PE.fastq $output/"${base}"_2.trimmed_SE.fastq LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20

         done
        """
        }


process runFastQCtrimmeddata {
        input: file step3_ch //one input
        output: 
        file("results") into step4_ch
        publishDir "results/TrimmedData"
        publishDir "results/FastQCtrimmeddata" //where should the results get linked to from the work folder?

        cpus 4
        script:
        base = in_file.baseName
        """
        #!/usr/bin/env bash

        mkdir $base

        fastqc -t 4 $step3_ch/*trimmed_PE.fq -o $base

        """

}

process runMultiQc {
	input: file f from step4_ch.collect() //one input
        output:
        file("results") into step5_ch //many outputs
	publishDir "results/MultiQCtrimmed"
"""	
#!/usr/bin/env bash

multiqc . -o results
"""
}

*/
