#!/usr/bin/env nextflow



/* Copyright 2019 University of the Witwatersrand, Johannesburg on behalf of the Pan-African Bioinformatics Network for H3Africa.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 * Contributors:
 *  Eli Moss, Mushal Allam, Ansie Yssel, Heyam Mohammed, Scott Hazelhurst, Mariem Hanachi, 
 *  Stanford Kwenda, Ami Bhatt, Gerrit Botha, Alfred Ssekagiri, Ovokeraye Oduaran
 *
 */


import org.apache.commons.lang3.StringUtils

params.input = "quality_control/*.fq.gz"

out_dir  = params.out_dir




if (!params.paired)  {
   paired="SE"
   inp_src = Channel.fromPath(params.input).map { file -> [file.simpleName, file] }
}  else {
   paired = "PE"
   inp_src = Channel.fromFilePairs(params.input)  
}



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
      fastqc -t ${params.num_cpus} $fq -o $base
    """
}


process runMultiQc {
   input: 
      file f from step1_ch.collect() //one input
   output:
     file("*") into step2_ch 
   publishDir "${out_dir}/multiqc"
    script:
    """
     #!/usr/bin/env bash
     multiqc $f -o .
   """
}

process runTrimmomatic{

   input: 
     set val(base), file(fq) from  inp2
   output:
     set val(base), file("${base}_trm*fqz") into step3_ch
   script:
     if (params.paired) {
       data="${fq[0]} ${fq[1]}"
       out ="${base}_trmP_1.fqz ${base}_trmP_2.fqz ${base}_trmU_1.fqz ${base}_trmU_2.fqz"
     } else {
        out = "${base}_trm.fqz"
        data = fq
     }
    """
    java org.usadellab.trimmomatic.Trimmomatic  $paired -threads 6 -trimlog ${base}.log \
            $data $out \
            ${params.leading} ${params.trailing} ${params.window} ${params.adapter} \
            ${params.minlen}
    """
}


process runFastQCtrimmeddata {
   input: 
     set val(base), file(f) from step3_ch
   output: 
     file("${base}/*") into step4_ch
   publishDir "${out_dir}/fastqc_trim"
     cpus params.num_cpus
     script:
   """
     #!/usr/bin/env bash
     mkdir $base
     fastqc -t ${params.num_cpus} $f -o $base
   """
}



process multiQcTrimmed {
  input: 
      file f from step4_ch.collect() //one input
  output:
     file("*") into step5_ch //many outputs
     publishDir "${out_dir}/multiqc_trim"
     """	
     #!/usr/bin/env bash
       multiqc . -o .
     """
}


