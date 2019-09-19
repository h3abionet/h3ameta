

org.apache.commons.lang3.StringUtils

params.input = "quality_control/*.fq.gz"


if (!params.paired)  {
   println "not implemented"
   System.exit(0)
}


trimmo_env = System.getenv("TRIMMO")

if  (params.trimmo.size()>0)
  trimmo = file(params.trimmo)
else 
  if ((trimmo_env.getClass() != org.codehaus.groovy.runtime.NullObject) && (trimmo_env.size()>0))
       trimmo = file(trimmo_env)
   else {
     println "I don't know where Trimmomatic is";
     println "Please specify in the config file or set the TRIMMO environment variable"
     System.exit(-1)
  }


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
     file trimmo
   output:
     set val(base), file("${base}_trm*fqz") into step3_ch
   script:
    """
    java -jar $trimmo $paired -threads 6 -trimlog ${base}.log \
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



process multiQcTrimmed {
  input: 
      file f from step4_ch.collect() //one input
  output:
     file("results") into step5_ch //many outputs
     publishDir "results/MultiQCtrimmed"
     """	
     #!/usr/bin/env bash
       multiqc . -o results
     """
}


