

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
   cpu  26
   input: 
    file(pileup)
   output: 
     file "${sample}.tsv"
   script:
     sample = pileup.baseName
		min_cvg=5,
		min_freq=0.8,
		min_qual=20
	script:
		"scripts/callSNPs.py"

# get consensus sequence from pileup
process snp_consensus {
	input: rules.call_snps.output
	output: "consensus/{sample}.txt"
	resources:
		mem=2,
		time=2
	threads: 1
	shell:
		"echo {wildcards.sample} > {output}; cut -f4 {input} >> {output}"

# combine consensus sequences into one file
process combine {
	input:
		dynamic("consensus/{sample}.txt")
	output: "{name}.cns.tsv".format(name = prefix)
	resources:
		mem=2,
		time=1
	threads: 1
	shell:
		"paste {input} > {output}"

# find positions that have a base call in each input genome and at least
# one variant in the set of input genomes
process core_snps {
	input: rules.combine.output
	output: "{name}.core_snps.tsv".format(name = prefix)
	resources:
		mem=16,
		time=1
	threads: 1
	script:
		"scripts/findCoreSNPs.py"

# convert core SNPs file to fasta format
process core_snps_to_fasta {
	input: rules.core_snps.output
	output: "{name}.fasta".format(name = prefix)
	resources:
		mem=16,
		time=1
	threads: 1
	script:
		"scripts/coreSNPs2fasta.py"

# perform multiple sequence alignment of fasta file
process multi_align {
	input: rules.core_snps_to_fasta.output
	output: "{name}.afa".format(name = prefix)
	resources:
		mem=200,
		time=12
	threads: 1
	shell:
		"muscle -in {input} -out {output}"

# calculate phylogenetic tree from multiple sequence alignment
process build_tree {
	input: rules.multi_align.output
	output: "{name}.tree".format(name = prefix)
	resources:
		mem=8,
		time=1
	threads: 1
	shell:
		"fasttree -nt {input} > {output}"

# plot phylogenetic tree
process plot_tree {
	input: rules.build_tree.output
	output: "{name}.tree.pdf".format(name = prefix)
	resources:
		mem=8,
		time=1
	threads: 1
	script:
		"scripts/renderTree.R"

# count pairwise SNVs between input samples
process pairwise_snvs {
	input: dynamic("consensus/{sample}.txt")
	output: "{name}.dist.tsv".format(name = prefix)
	resources:
		mem=8,
		time=1
	threads: 1
	script:
		"scripts/pairwiseDist.py"
*/
