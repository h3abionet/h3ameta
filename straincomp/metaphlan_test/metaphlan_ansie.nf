#!/usr/bin/env nextflow


params.in = '/home/ansieyssel/h3ameta/test_datasets/taxonomic_classification/*.f*q*' /*modified that it can take varois formats for input also unzipped, also * added to path so it can be run by various users*/
params.mpa = '/home/mushalali/bin/miniconda3/bin/db_v20'

sequencing_data = file(params.in)
mpa_db = file(params.mpa)


process MetaPhlAn2 {
    publishDir 'metaphlan2_outs/', mode: 'copy', overwrite: true
    input:
    file (infile) from sequencing_data
    output: 
    file "${infile}" into report_in_data_ch 
    file "${infile.baseName}_MetaPhlAn2_microbes_list.tsv" into FunctionalProfiling_ch, report_ch1
    file "${infile.baseName}_bt2out.txt" into Bowtie2Output_ch, report_ch2
    cpus 2
    time '4h'
    memory '4GB'
    container = 'quay.io/biocontainers/metaphlan2:2.7.8--py27h24bf2e0_0'
    
    script:
    """ 
    #!/usr/bin/env bash
     metaphlan2.py --input_type fastq --mpa_pkl ${mpa_db}/mpa_v20_m200.pkl --bowtie2db ${mpa_db}/mpa_v20_m200 --bowtie2out ${infile.baseName}_bt2out.txt --nproc ${task.cpus} $infile > ${infile.baseName}_MetaPhl#An2_microbes_list.tsv
     merge_metaphlan_tables.py /home/ansieyssel/h3ameta/straincomp/metaphlan_test/metaphlan2_outs/*_bt2out.txt > merged_abundance_table.txt
     metaphlan_hclust_heatmap.py -c bbcry --top 25 --minv 0.1 -s log --in merged_abundance_table.txt --out abundance_heatmap.png
    
    
    """
}

process Make_report {
	publishDir 'metaphlan2_outs/', mode: 'copy', overwrite: false
	input:
	file (infile) from report_in_data_ch
	file (outfile1) from report_ch1
	file (outfile) from report_ch2
	output:
	file "report.txt" into Metaphlan_report_ch
	cpus 1
	
	script:
	"""
	#!/usr/bin env bash
	echo "the input file was:" > report.txt
	echo infile > report.txt
	echo "the output files are:" > report.txt
        echo outfile1 > report.txt
	echo outfile2 > report.txt 
	echo "#############"
	"""
 
}
