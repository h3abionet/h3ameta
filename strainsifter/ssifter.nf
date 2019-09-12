

prefix = params.out_dir


reference = file(params.reference)

def getFQBase = { 
  file ->
    base = file.baseName;
    m =  base =~ /(.*)_[12].fastq/
    return m[0][1]
}


samples_ch  = Channel.fromFilePairs(params.input_pat) 
     {file -> getFQBase(file) }



process bwa_index {
 input: 
   file (ref) from reference_ch
 output:
   file	"${ref}.{amb,ann,bwt,pac,sa}" into ref_ind_ch
 script:
   "bwa index $ref"
}


process bwa_align {
   memory '6GB'
   time  '6h'
   cpus 10   // should be a param and the thread arguments below dynamically chosen
   input:
     file(ref) from reference_ch
     file("*") from ref_ind_ch
     set val(sample), file(seqs_12) from samples_ch
   output:
     set val(sample), file(outname) into (filtered_ch, filtered_cov_ch)
   script:
     outname = "${sample}_filtered.bam"
     """
     bwa mem -t 9  $ref ${seqs_12} | samtools view -b -q ${params.mapq} |\
         bamtools filter -tag 'NM:<=${params.n_mismatches}' | \
         samtools sort --threads 9 -o $outname
     """
}


process genomecov {
     label 'bigmem'
     input:
       set val(sample), file(filtered) from filtered_cov_ch
     output:
       set val(sample), file(outname) into coverage_ch
     script:
	outname = filtered.baseName+".tsv"
	"bedtools genomecov -ibam $filtered > ${outname}"
}



process calc_coverage {
    label = 'bigmem'
    input:
       set val(sample), file(coverage) from coverage_ch	
    output:
       set val(sample), file(cvg_calc) into calc_coverage_ch
    script:
	cvg_calc = coverage.baseName+".cvg"
	"getCoverage.py $coverage ${params.min_cvg} $cvg_calc"
}



// Now we filter out the BAMs that do not meet our coverage requirement
// we join the BAM file and coverage scores and then do the filter
passed_ch = filtered_ch.join(calc_coverage_ch).filter { 
  label, bam, cvg -> 
    (s_cvg, s_pct) = cvg.text.split("\t")
    println "${label}, ${s_cvg}"
    return (s_cvg.toFloat() > params.min_cvg) && (s_pct.toFloat()>=params.min_genome_percent)
}
                                       


process faidx {
   cpus 2
   time '1h'
   input:
     file(reference)
   output: 
     set file(reference), file("${reference}.fai") into fai_ch
   script:
      "samtools faidx $reference"
}


process pileup {
    memory "32GB"
    cpus 16
    input:
       set file(ref), file(ref_idx) fom fai_ch
       set val (label), file(passed), file(coverage) from passed_ch
    output: 
	file("${output}") into pileup
    script:
      output = passed.baseName+".pileup"
      "samtools mpileup -f ${ref} -B -aa -o ${output} $passed"
}



process call_snps {
   time '2h'
   mem  '32GB'
   cpus  16
   input: 
    file(pileup)
   output: 
     file call
   script:
     call  = "${pileup.baseName}.tsv"
     script:
	"""
        callSNPs.py  $pileup $call
        """
}


// get consensus sequence from pileup
process snp_consensus {
  mem '2GB'
  time '2h'
  cpus 1
  input: 
     file(call)
  output: 
     file(consensus) into (consensus1, consensus2, consensus3)
  script:
      out = call.baseName+".txt"
      "echo $call > ${output}; cut -f4 ${call} >> ${output}"
}




// combine consensus sequences into one file
process combine {
  input:
    file(consensi) from consensus1.toList()
  output: 
    file output into consensus_combine
  script:
    output = "${params.label}.tsv"
    "paste ${consensi}} > ${output}"
}


// find positions that have a base call in each input genome and at least
// one variant in the set of input genomes
process core_snps {
   mem '16GB'
   input: 
     file(consensus_combine)
   output: 
     file(out) into core
   script:
    out = "${consensus_combine.baseName}_core_snps.tsv"     
    """
    findCoreSNPs.py $consensus_combine $out
    """
}
// convert core SNPs file to fasta format
process core_snps_to_fasta {
  mem '16GB'
  input: 
    file(core)
  output:
    file(core_fa) 
  script:
    core_fa = core.baseName + ".fa"
    "coreSNPs2fasta.py $core $core_fa"
}

// perform multiple sequence alignment of fasta file
process multi_align {
  mem '200GB'
  input: 
    file(core_fa)
  output: 
    file(aligned)
  script:
     aligned = "${core_fa.baseName}.aln"
     "muscle -in ${core_fa} -out ${aligned}"
}

// calculate phylogenetic tree from multiple sequence alignment
process build_tree {
  mem '8GB'
  input:
   file(aligned)
  output:
   file(tree)
  script:
   tree = "${aligned.baseName}.tree"
   "fasttree -nt ${aligned} > ${tree}"
}

// plot phylogenetic tree
process plot_tree {
  input: 
    file(tree)
  output:
    file(pdf)
  script:
    pdf = "${tree.baseName}.pdf"
   script:
     "renderTree.R $tree $pdf"
}

pairs = consensus2.merge(consensus2).filter( a, b -> a.baseName == b.baseName )

// count pairwise SNVs between input samples
process pairwise_snvs_indiv {
   input:
     set file(a), file(b) from pairs
  output:
     stdout into result
  shell:
     '''
        totalBases=`wc -l !{a}`
        totalBases=$( totalBases - 1 )
        totalPos=`paste  !{a} !{b} | sed '1d' | grep -v N | wc -l`
        diffPos=`paste  !{a} !{b} | sed '1d' | grep -v N | \
          awk '$1 != $2 {print $0}' | wc -l`
        echo "!{a}\t!{b}\t!{diffPos}\t!{totalPos}\t!{totalBases}" 
     '''
}

process combine_dists {
  input:
    val(combine) from  pairs.toSortedList()
  output:
    file("result_table.tsv")
  script:
    table = combine.join("\n")  
    out = new File("result_table.tsv")
    out.write(table)
}
