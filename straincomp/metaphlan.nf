#!/usr/bin/env nextflow

Channel.fromPath(params.inp).into { sequencing_data; report_in_data_ch }

mpa_db = file(params.mpa)


process MetaPhlAn2 {
   cpus 2
   time '4h'
   memory '4GB'
   input:
      file (infile) from sequencing_data
   output: 
      file "$out" into (bowtie_out_ch,metaphl_tables_ch)
   publishDir "${params.out}", mode: 'copy', overwrite: true
   script:
       out="${infile.baseName}_MetaPhlAn2_microbes_list.tsv"
       """ 
       #!/usr/bin/env bash
	metaphlan2.py --input_type fastq --mpa_pkl ${mpa_db}/mpa_v20_m200.pkl \
		      --bowtie2db ${mpa_db}/mpa_v20_m200 \
		      --bowtie2out ${infile.baseName}_bt2out.txt \
		      --nproc ${task.cpus} $infile\
			> $out
       """
}



process mergeTables {
  
  input:
    file (table) from bowtie_out_ch.toList()
  output:
    file (out) into (cluster_ch, report_merged_abundance_ch)
  publishDir "${params.out}", mode: 'copy', overwrite: true
  script:
    out = "merged_abundance_table.tsv"
   """
     merge_metaphlan_tables.py $table > $out
   """
}


process heatMap {
  input:
     file(abundance) from cluster_ch
  publishDir "${params.out}", mode: 'copy', overwrite: false
  output:
     file(heatmap) into heatmap_ch
  publishDir "${params.out}", mode: 'copy', overwrite: true
  script:
    heatmap = "abundance_heatmap.png"
    """
     metaphlan_hclust_heatmap.py -c bbcry --top 25 --minv 0.1 -s log \
              --in $abundance --out $heatmap
   """
}


process makeReport {
   publishDir "${params.out}", mode: 'copy', overwrite: false
   input:
      file (infile) from report_in_data_ch.toList()
      file (table) from report_merged_abundance_ch
      file (heatmap) from heatmap_ch
      file (metatabs) from metaphl_tables_ch.toList()
   output:
     file "report.txt" into metaphlan_report_ch
     cpus 1

   script:
   """
   #!/usr/bin/env bash
   echo "The input files were:" > report.txt
   for fn in ${infile}; do
       echo "  \$fn " >> report.txt
   done
   echo " " >> report.txt # blank like
   echo "The abundance table is $table:" >> report.txt
   echo "The heatmap file is $heatmap " >> report.txt 
   echo "#############" >> report.txt
   echo " " >> report.txt # blank like
   echo "The individual Metaphlan taxonomy files were:" >> report.txt
   for fn in ${metatabs}; do
       echo "  \$fn " >> report.txt
   done
   """

}

