


params.fq="data/*{1,2}.fastq.gz"	// input directory for fastq files

params.krakenDB="krakenDB"	 // Path to kraken DB
params.genome="hostDB/host" // path to human genome
params.viral_genome_dir="viralDB/viruses"
params.katDB="katDB"	// # will symlink to Gerrit's /spaces/gerrit ...

// create a nextflow channel
input_fq = Channel.fromFilePairs("${params.fq}*_{1,2}.fastq.gz")
// Note we use fromFilePairs rather than two separate fromPath because it's slightly
// simpler and also allows us to do many at the same time -- if you have separate fromPath
// then you need to be sure that the two channels synchronise


// channel for kraken DB
krakenDB = file(params.krakenDB)

// channel for KAT kmer DB
kat_db = params.katDB

// channel for host genome
host_genome = file(params.genome)

// channel for host genome
viral_genomes = file(params.viral_genome_dir)



// kat requires an uncompressed file so check whether we need to 
// uncompress first -- return command to uncompress and resultant file name
// in theory, it would be better to do this uncompressing in a separate process because
// one could then parallelise, resume, but the nextflow gymnastics are complex and
// the cost is actually very low compared to running kat so we sick with the simpler
// approach
def check_input_compress = { file ->
  if (['gz','fqz'].contains(file.getExtension()))
    return ["gunzip -c $file > ${file.baseName}", file.baseName]
  else
    return ["", file]
}






process removeHostReads {
   input:
     set val(name), file(fq) from input_fq
   output:
     set val(name), file("clean*1*.*q"),file("clean*2*.*q") into (clean_fq_bt, clean_fq_kr)
   script:
     (uncompressA, fa) = check_input_compress(fq[0])
     (uncompressB, fb) = check_input_compress(fq[1])
     """
     $uncompressA
     $uncompressB
     kat filter seq --output_prefix clean --invert --seq $fa --seq2 $fb $kat_db
     """
}



process runBowtie2{
    publishDir params.out_dir, mode: "copy", overwrite: false

    input:
      set val(name), file(cfq1), file(cfq2) from clean_fq_bt
      file(viral_genomes)
    output:
      set val(outname), val(outfname) into bowtie2_aligned
    script:
     outfname = "${name}.sam"
     outname  = "${name}"
     """
       bowtie2 -x $viral_genomes/${params.viral_genome} \
               -1 ${cfq1} -2 ${cfq2} -S $outfname
     """
}

process getMappingstats {
    publishDir params.out_dir, mode: "copy", overwrite: false

    input:
       set val(sample), file(aligned) from bowtie2_aligned
    output:
       file(out) into mappingStats
    script:
    out = "${sample}_mapstats.txt"
    """
     samtools view -S -b ${aligned} | samtools sort -n -o sample.sorted.bam
     samtools index sample.sorted.bam
     samtools idxstats sample.sorted.bam >> $out
     #samtools flagstat ${aligned}>>mappingStats.txt
     #samtools view $aligned | cut -f3 | sort | uniq -c | \
     #awk -v i=$aligned '{print "\t"\$1"\t"\$2}' >> mappingStats.txt
    """


}


process runKraken {
   cpus 2

   input: 
      set val(name),  file(fq1), file(fq2) from clean_fq_kr
      file (krakenDB)

   output:
      set file(report), file(kraken_out) into kraken_classified

   publishDir params.out_dir, mode: "copy", overwrite: false
   script:
     report     = "report.tsv"
     kraken_out = "kraken_out.txt"
     """
       kraken2 --report $report --db $krakenDB --threads $task.cpus\
               --fastq-input --paired $fq1 $fq2 >  $kraken_out
    """
}

