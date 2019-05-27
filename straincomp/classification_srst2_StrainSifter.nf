#!/usr/bin/env nextflow
/*
(c) University of the Witwatersrand, Johannesburg on behalf of the H3A Bioinformatics Network, 2019
Contributors:
   Eli Moss, Mushal Allam, Ansie Yssel, Heyam Mohammed, Scott Hazelhurst

and finally visualization with a custom bargraph script written in R, as well as Krona.
The accompanying nextflow.config file allows this workflow to be run on the Stanford SCG cluster with the
SLURM scheduler.
 ~/nextflow ../wits_workshop/nextflow/taxonomic_classification/taxonomic_classification.nf  --tax_level S -resume \
 -profile wits --in ../wits_workshop/nextflow/test_data/*.fq
*/

//The parameters below can all be overridden with --parametername on the commandline (e.g. --in or --dataset_table)

Channel.fromPath(params.inp).into {data1; data2}
kraken_db     = file(params.db)
dataset_table = file(params.dataset_table)
annot_db      = file(params.annot_db)
strainSifter_config = file(params.strainSifter_config)





if (params.paired) {
    paired   = "--paired"
    in_srst2 = "--input_pe"
} else {
    paired = " "
    in_srst2 = "--input_se"
}


process kraken {
    cpus 2
    time '48h'
    memory '20GB' //20
    input:
       file sample_fq from data1
       file kraken_db
    output: 
       file "${sample_fq.baseName}_kraken.tsv" into kraken_out 
    script:
    """
    #!/usr/bin/env bash
    #note that the code below can be written in any language, so long as the interpreter is named
    #in the above shebang line.
    #also note: variables from this or higher scopes can be named inside strings with a dollar sign, as below.
    #when part of a larger name, as in the tsv output filename below, the variable must be enclosed with {}
    # the final $sample_fq names the input.
    kraken2 --db $kraken_db/ --threads $task.cpus $paired\
    --report ${sample_fq.baseName}_kraken.tsv \
    --quick --memory-mapping \
    $sample_fq
    """
}

process bracken {
   cpus 1
   time '48h'
   memory { 8.GB * task.attempt }
   input: 
      file f from kraken_out 
   output: 
      file "${f.baseName}_bracken.tsv" into bracken_out 
   publishDir '${params.out}', mode: 'copy', overwrite: true 
   script:
   """
	#!/usr/bin/env bash
	bracken -d $params.db -i $f \
	-o ${f.baseName}_bracken.tsv -r $params.readlen -l $params.tax_level \
   """
}

bracken_out.into{bracken_out1; bracken_out2}


process krona {
    cpus 1
    time '1h'
    memory '1GB'
    input: 
      file k from bracken_out2
    output: 
      file "krona_${k}.html" into krona_out
    publishDir "${params.out}", mode: 'copy', overwrite: true
    """
    ktImportTaxonomy -m 3 -s 0 -q 0 -t 5 -i ${k} -o krona_${k}.html \
     -tax /global/taxonomy
    """
}

process collect_results {
   cpus 1
   time '1h'
   memory {2.GB * task.attempt } //dynamic resource allocation!
   publishDir "${params.out}", mode: 'copy', overwrite: true
   input:
       file data from dataset_table
       file f from bracken_out1.collect() 
   output: file 'class_long.tsv' into collect_results_out
   script:
   """
   collate_results.py $params.tax_level $data class_long.tsv $f
   """
}

process barplot {
   input: 
     file f from collect_results_out
   output: 
     file out_pdf into barplot_out
   publishDir "${params.out}", mode: 'copy', overwrite: true
   script:
      out_pdf = 'relative_abundance_barplot.pdf'
      """
      composition_barplot.R $f $out_pdf
      """
}

process srst2{
   cpus 1
   time '1h'
   memory '4GB'
   module 'samtools18'
   input: 
    file sample_fq from data2 //input channel is a file, as declared above
    file annot_db
   output: 
     file outfname into srst2_annot_db_out
   publishDir "${params.out}", mode: 'copy', overwrite: true
   script:
   outfname = "*_results.txt"
   """
   #!/usr/bin/env bash
   srst2 ${in_srst2} $sample_fq --output ${sample_fq.baseName}_srst2 --gene_db $annot_db
   """   
}

process StrainSifter {
        cpus 1
        time '96h'
        label 'bigmem'
        input:             
            file strainSifter_config
        output: 
           file pdf  into StrainSifter_out
        publishDir "${params.out}", mode: 'copy', overwrite: true
        script:
	pdf = "nature.tree.pdf"
        """
        echo Need to handle this >  $pdf
        """
}
