#!/usr/bin/env nextflow


params.in = '/home/ansieyssel/h3ameta/test_datasets/taxonomic_classification/*.f*q*' 
/*modified that it can take varois formats for input also unzipped, also * added to path so it can be run by various users*/
params.mpa = '/opt/exp_soft/bioinf/metaphlan2/db_v20'

Channel.fromPath(params.in).into { sequencing_data; report_in_data_ch }

mpa_db = file(params.mpa)


process MetaPhlAn2 {
    cpus 2
    time '4h'
    memory '4GB'
    publishDir 'metaphlan2_outs/', mode: 'copy', overwrite: true
    input:
      file (infile) from sequencing_data
    output: 
       file "$out" into bowtie_out_ch
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
  script:
    out = "merged_abundance_table.tsv"
   """
     merge_metaphlan_tables.py $table > $out
   """
}


process heatMap {
  input:
     file(abundance) from cluster_ch
  publishDir 'metaphlan2_outs/', mode: 'copy', overwrite: false
  output:
     file(heatmap) into heatmap_ch
  script:
    heatmap = "abundance_heatmap.png"
    """
     metaphlan_hclust_heatmap.py -c bbcry --top 25 --minv 0.1 -s log \
              --in $abundance --out $heatmap
   """
}


process makeReport {
	publishDir 'metaphlan2_outs/', mode: 'copy', overwrite: false
	input:
           file (infile) from report_in_data_ch.toList()
	   file (table) from report_merged_abundance_ch
	   file (heatmap) from heatmap_ch
	output:
  	  file "report.txt" into Metaphlan_report_ch
	  cpus 1
	
	script:
	"""
	#!/usr/bin/env bash
	echo "The input files were:" > report.txt
	echo $infile >> report.txt
	echo >> report.ch # blank like
	echo "The abundance table is $table:" >> report.txt
	echo "The heatmap file is $heatmap " >> report.txt 
	echo "#############" >> report_ch
	"""
 
}

