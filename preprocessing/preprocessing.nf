

org.apache.commons.lang3.StringUtils

params.input = "quality_control/*.fq.gz"


if (!params.paired)  {
   println "not implemented"
   System.exit(0)
}


Channel.fromFilePairs(params.input)
       .into { qc_input; data_input }



process runFastQCOriginal {
   cpus 4
   input: 
     set val(base), file(fq) from qc_input
   output: 
     set file("$base/*.zip"), file("$base/*.html") into step1_ch //many outputs
   publishDir "${params.out_dir}/fastqc"
   script:
       base = in_file.baseName
       """
          mkdir $base
  	  fastqc -t 4 $fq -o $base
       """
}


process runMultiQc {
   input: 
      file f from step1_ch.collect() //one input
   output:
     file("results") into step2_ch 
   publishDir "${params.out_dir}/MultiQcOriginal"
    script:
    """
     #!/usr/bin/env bash
     multiqc $f -o results
   """
}

process runTrimmomatic{
   cpus 6
   input: 
      set val(base), file(fq) from data_input
   output:
     file("results") into trimmed
   script:
    """
    java -jar $trimmo_jar $paired -threads 6 -trimlog ${base}.log ${fq[0]} ${fq[1]} \
            ${base}_trmP_1.fqz ${base}_trmP_2.fqz ${base}_trmU_1.fqz ${base}_trmU_2.fqz \
            LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20
    """
}


process runFastQCtrimmeddata {
   input: 
      file step3_ch //one input
        output: 
          file("results") into step4_ch
        publishDir "results/TrimmedData"
        publishDir "results/FastQCtrimmeddata" 

        cpus 4
        script:
        base = in_file
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
