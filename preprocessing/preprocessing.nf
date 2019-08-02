

import org.apache.commons.lang3.StringUtils

params.input = "quality_control/*_{1,2}.fq.gz"

trimmo_jar="/opt/exp_soft/bioinf/Trimmomatic-0.36/trimmomatic-0.36.jar"

if (!params.paired)  {
   println "not implemented"
   System.exit(0)
}  else {
   paired = "PE"
}


inp_src = Channel.fromFilePairs(params.input)  

inp1 = Channel.create()
inp2 = Channel.create()
inp_src.separate (inp1, inp2) { x -> [x, x] }


process runFastQCOriginal {
   cpus 4
   input: 
     set val(base), file(fq) from inp1
   output: 
     file("$base/*.{zip,html}") into step1_ch //many outputs
   publishDir "${params.out_dir}/fastqc"
   script:
    """
      mkdir $base
      fastqc -t 4 $fq -o $base
    """
}


process runMultiQc {
   input: 
      file f from step1_ch.collect() //one input
   output:
     file("results/*") into step2_ch 
   publishDir "${params.out_dir}/MultiQcOriginal"
    script:
    """
     #!/usr/bin/env bash
     mkdir results
     multiqc $f -o results
   """
}

process runTrimmomatic{

   input: 
     set val(base), file(fq) from  inp2
   output:
     set val(base), file("${base}_trm*fqz") into step3_ch
   script:
    """
    java -jar $trimmo_jar $paired -threads 6 -trimlog ${base}.log \
            ${fq[0]} ${fq[1]} \
            ${base}_trmP_1.fqz ${base}_trmP_2.fqz\
            ${base}_trmU_1.fqz ${base}_trmU_2.fqz \
            LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20
    """
}


process runFastQCtrimmeddata {
     input: 
	 set val(base), file(f) from step3_ch
     output: 
     file("${base}/*") into step4_ch
     publishDir "${params.out_dir}/FastQCtrimmeddata" 
     cpus 4
     script:
     """
     #!/usr/bin/env bash
     mkdir $base
     fastqc -t 4 $f -o $base
     """
}


