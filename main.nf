#!/usr/bin/env nextflow
println "clear".execute().text

/*  ======================================================================================================
 *  HELP MENU
 *  ======================================================================================================
 */
line="=".multiply(100)
ver="nf-h3ameta v0.1"

if (params.help) {
    println "\n${line}"
    println "#".multiply(48 - ("${ver}".size() / 2 )) + "  ${ver}   " + "#".multiply(48 - ("${ver}".size() / 2 ))
    println "${line}\n"
    println "USAGE:"
    println "nextflow run nf-rnaSeqCount -profile \"slurm\" --data \"/path/to/data\" --genome \"/path/to/genome.fa\" --genes \"/path/to/genes.gtf\"\n" 
    println "HELP:"
    println "nextflow run nf-rnaSeqCount --help\n"
    println "MANDATORY ARGUEMENTS:"
    exit 1
}

// MAIN USER INPUT ERRORS
data_error = """
${line}
Oooh no!! Looks like there's a serious issue in your command! 
I do not recognise the \'--data ${params.data}\' option you have given me, or you have not given me any \'--data\' option at all!
Please provide a valid directory with you input FASTQ reads with the \'--data\' option to run the nf-rnaSeqCount workflow! 
${line}
"""

genome_error = """
${line}
Oooh no!! Looks like there's a serious issue in your command! 
I do not recognise the \'--genome ${params.genome}\' option you have given me, or you have not given me any \'--genome\' option at all!
Please provide a valid FASTA file (.fasta or .fa) for your reference genome with the \'--genome\' option to run the nf-rnaSeqCount workflow! 
${line}
"""

kraken_db_error = """
${line}
Oooh no!! Looks like there's a serious issue in your command! 
I do not recognise the \'--kraken_db ${params.kraken_db}\' option you have given me, or you have not given me any \'--kraken_db\' option at all!
Please provide a valid directory with your Kraken database with the \'--kraken_db\' option to run the nf-rnaSeqMetagen workflow! 
${line}
"""

// genes_error = """
// ${line}
// Oooh no!! Looks like there's a serious issue in your command! 
// I do not recognise the \'--genes ${params.genes}\' option you have given me, or you have not given me any \'--genes\' option at all!
// Please provide a valid GTF annotation file (.gtf) for your reference genome with the \'--genes\' option to run the nf-rnaSeqCount workflow! 
// ${line}
// """

mode_error = """
${line}
Oooh no!! Looks like there's an serious issue in your command! 
I do not recognise the \'--mode ${params.mode}\' option you have given me, or you have not given me any \'--mode\' option at all!
The allowed options for \'--mode\' are:
\tprep.Containers\t\t: For downloading Singularity containers used in this workflow.
\tprep.STARIndex\t\t: For indexing your reference genome using STAR.
\tprep.BowtieIndex\t: For indexing your reference genome using Bowtie2.
\trun.ReadQC\t\t: For performing general QC on your reads using FastQC. 
\trun.ReadTrimming\t: For trimming low quality bases and removing adapters from your reads using Trimmmomatic.
\trun.ReadAlignment\t: For aligning your reads to your reference genome using STAR.
\trun.ReadCounting\t: For counting features in your reads using HTSeq-count and featureCounts.
\trun.MultiQC\t\t: For getting a summary of QC through the analysis using MultiQC.
\nPlease use one of the above options with \'--mode\' to run the nf-rnaSeqCount workflow!
${line}
"""
from_error = """
${line}
Oooh no!! Looks like there's an serious issue in your command! 
I do not recognise the \'--from ${params.from}\' option you have given me!
The allowed options for \'--from\' are:
\trun.ReadQC\t\t: To resume from the QC step.
\trun.ReadTrimming\t: To resume from the trimming step.
\nPlease use one of the above options with \'--from\' to run the nf-rnaSeqCount workflow!
${line}
"""
generic_error = """
${line}
OOOHHH!!!! YOU DONE FUCKED UP!!!
${line}
"""

// EMPTY LIST FOR COLLECTING ALL THE PATHS TO BIND TO SINGULARITY IMAGE
bind_dirs = []

// USER PARAMETER INPUT: DATA DIRECTORY
switch (params.data) {
    case [null]:
        data_dir = "NOT SPECIFIED!"
        break
    default:
        data_dir = file(params.data, type: 'dir')
        bind_dirs.add(data_dir)
        break
}

// USER PARAMETER INPUT: GENOME FASTA FILE
switch (params.genome) {
    case [null]:
        genome = "NOT SPECIFIED!"
        break
    default:
        genome = file(params.genome, type: 'file', checkIfExists: true)
        index_dir = genome.getParent()
        bind_dirs.add(genome.getParent())
        break
}

// USER PARAMETER INPUT: OUTPUT DIRECTORY
switch (params.out) {
    case [null]:
        out_dir = file("${PWD}/results_nf-h3ameta", type: 'dir')
        break
    default:
        out_dir = file(params.out, type: 'dir')
        bind_dirs.add(out_dir)
        break
}

// USER PARAMETER INPUT: KRAKEN2 DB DIRECTORY
switch (params.kraken_db) {
    case [null]:
        kraken_db = "NOT SPECIFIED!"
        break
    default:
        kraken_db = file(params.kraken_db, type: 'dir')
        kraken_db.mkdir()
        taxonomy = file("$kraken_db/taxonomy", type: 'dir')
        // names_dmp = file("${taxonomy}/names.dmp", type: 'file')
        bind_dirs.add(kraken_db)
}

// USER INPUT RESUME FROM: WHERE TO PICK UP??
if(params.from == null) {
    resume_from = "NOT SPECIFIED YET!"
} else if(params.from in ["run.ReadTrimming", "run.ReadQC"]) {
    resume_from = params.from
} else {
    exit 1, "$from_error"
}

// USER STRANDED MODE: ARE WE DOING PAIRED- OR SINGLE-END?
if (params.singleEnd == null && params.pairedEnd == null) {
    stranded = "paired-end"
} else if(params.singleEnd) {
    stranded = "single-end"
} else if(params.pairedEnd){
    stranded = "paired-end"
} else {}

// FUNCTIONS FOR CHECKING DATA
def breakIfNull(parameter,error) {
    if (parameter == null) {
        exit 1, error
    } else {}
}

def checkDataStrandedness() {
    switch (stranded) {
        case ["paired-end"]:
            read_pairs = Channel.fromFilePairs("${data_dir}/*{R,read}[1,2]*.{${ext}}", type: 'file')
                .ifEmpty { exit 1, "$main_data_error" }
            break
        case["single-end"]:
            read_pairs = Channel.fromFilePairs("${data_dir}/*.{${ext}}", type: 'file', size:1)
                .ifEmpty { exit 1, "$main_data_error" }
            break
    }
}

// USER PARAMETER INPUT: TRIMMOMATIC OPTIONS
switch (params.trim) {
    case [null]:
        switch (stranded) {
            case ["paired-end"]:
                trim_params = "ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:8:keepBothReads TRAILING:28 MINLEN:40"
                break
            case ["single-end"]:
                trim_params = "ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 TRAILING:28 MINLEN:40"
                break
        }
        break
    default:
        trim_params = params.trim
        break
}

// OUTPUT DIRECTORIES
qc_dir         = file("${out_dir}/Read_QC", type: 'dir')
trim_dir       = file("${out_dir}/Read_Trimming", type: 'dir')
tax_clas_dir   = file("${out_dir}/Taxonomy_Classification", type: 'dir')
str_comp_dir   = file("${out_dir}/Strain_Comp", type: 'dir') 
vir_det_l_dir  = file("${out_dir}/Viral_Detect_Long", type: 'dir')
vir_det_s_dir  = file("${out_dir}/Viral_Detect_Short", type: 'dir')

ext = "fastq,fastq.gz,fastq.bz2,fq,fq.gz,fq.bz2"

// DATA ABSENT ERRORS
main_data_error = """
${line}
Oooh no!! Looks like there's a serious issue in your input data! There are no FASTQ file in the directory:
\t${data_dir}
Please ensure that you have given me the correct directory for you FASTQ input reads using the \'--data\' option!
${line}
"""

// USER INPUT MODE: WHICH ANALYSIS TO RUN!
switch (params.mode) {
    case [null]:
        exit 1, "$mode_error"
    
    case ["prep.Containers", "prep.Reference", "prep.KrakenDB", "prep.BrackenDB"]:
        mode = params.mode
        switch (mode) {
            case ["prep.Containers"]:
                break
            case ["prep.Reference"]:
                breakIfNull(params.genome,"$genome_error")
                break
            case ["prep.KrakenDB","prep.BrackenDB"]:
                breakIfNull(params.kraken_db,"$kraken_db_error")
                break
        }
        break

    // ONLY REQUIRE INPUT READS
    case ["run.ReadQC", "run.ReadTrimming"]:
        mode = params.mode
        breakIfNull(params.data,"$data_error")
        switch (mode) {
            case ["run.ReadQC"]:
                read_pairs = Channel.fromFilePairs("${data_dir}/*.{${ext}}", type: 'file', size:-1) { "all_reads" }
                    .ifEmpty { exit 1, "$main_data_error" }
                break
            case ["run.ReadTrimming"]:
                checkDataStrandedness()                
                break
        }
        break

    case ["run.StrainComp","run.TaxonomicClassification","run.ViralDetectLong","run.ViralDetectShort"]:
        mode = params.mode
        breakIfNull(params.data,"$data_error")
        breakIfNull(params.kraken_db,"$kraken_db_error")
        switch (mode) {
            case ["run.StrainComp"]:
                breakIfNull(params.dataset_table,"$generic_error")
                breakIfNull(params.annot_db,"$generic_error")
                breakIfNull(params.readlen,"$generic_error")
                breakIfNull(params.tax_level,"$generic_error")

                dataset_table = file(params.dataset_table, type: 'file')
                annot_db = file(params.annot_db, type: 'file')
                readlen = params.readlen
                tax_level = params.tax_level
                checkDataStrandedness()
                break
            case ["run.TaxonomicClassification"]:
                checkDataStrandedness()
                break
            case ["run.ViralDetectLong"]:
                breakIfNull(params.read_type,"$generic_error")
                breakIfNull(params.decontam_db,"$generic_error")
                breakIfNull(params.viral_db,"$generic_error")
                breakIfNull(params.acc_2_tax,"$generic_error")
                breakIfNull(params.names_dmp,"$generic_error")
                
                viral_db = file(params.viral_db, type: 'file')
                decontam_db = file(params.decontam_db, type: 'file')
                acc_2_tax = file(params.acc_2_tax, type: 'file')
                names_dmp = file(params.names_dmp, type: 'file')
                
                switch (params.read_type) {
                    case ["nanopore"]:
                        read_type = "map-ont"
                        break
                    case ["pacbio"]:
                        read_type = "map-pb"
                        break
                    default:
                        exit 1, "$generic_error"
                }
                checkDataStrandedness()
                break
            case ["run.ViralDetectShort"]:
                breakIfNull(params.viral_db,"$generic_error")
                checkDataStrandedness()
                viral_db = file(params.viral_db, type: 'file')
                break
        }
        out_dir.mkdir()
        break
    default:
        exit 1, "$mode_error"
        break
}

// USER PARAMETER INPUT: PATHS TO BE BINDED TO THE IMAGE
bind_dirs = bind_dirs
    .unique()
    .collect { it -> "-B ${it}"}
    .join("\n" + ' '.multiply(26))
    .toString()

// ========

switch (mode) {
        
        // MODE 1 - DOWNLOAD SINGULARITY CONTAINERS
    case ['prep.Containers']:
        process run_DownloadContainers {
            label 'mini'
            tag { "Downloading h3ameta Singularity image!" }
            publishDir "$PWD/containers", mode: 'copy', overwrite: true
            
            output:
            file("*.sif") into containers
            
            """
            singularity pull nf-h3ameta.sif shub://phelelani/nf-h3ameta:h3ameta
            """
        }
        break
        
        // MODE 2 - INDEX REFERENCE GENOME
    case ['prep.Reference']:
        process run_BWA {
            label 'midi'
            tag { "Index Reference (BWA)" }
            publishDir "$index_dir", mode: 'copy', overwrite: true

            """
            bwa index ${genome}
            samtools faidx ${genome}
            """
        }

        process run_Bowtie2 {
            label 'maxi'
            tag { "Index Reference (Bowtie2" }
            publishDir "$index_dir", mode: 'copy', overwrite: true
            
            output:
            set val("Bowtie2Index"), file("*") into out_bowtie_index
            
            """
             bowtie2-build --threads ${task.cpus} ${genome} ${genome.getBaseName()}
            """
        }
        break

        // MODE 3 - GENERATE KRAKEN DATABASE
    case ['prep.KrakenDB']:
        process run_GenerateKrakenDB {
            label 'maxi'
            tag { "Generate Kraken DB" }
            publishDir "$kraken_db", mode: 'copy', overwrite: false
            
            output:
            file("*") into out_kraken_db
            file("taxonomy/taxdump.tar.gz") into taxonomy_dump
            
            """
            kraken2-build --standard --threads ${task.cpus} --db .
            """
        }

        process run_UpdateTaxonomy {
            label 'mini'
            tag { "Update NCBI Taxonomy" }
            publishDir "$taxonomy", mode: 'copy', overwrite: true
            
            input:
            file(dmp) from taxonomy_dump
            
            output:
            file("*") into taxonomy_update
            
            """
            /opt/KronaTools-2.7/updateTaxonomy.sh --only-build --preserve .
            """
        }
        break

        // MODE 4 - GENERATE BRAKEN DATABASE
    case ['prep.BrackenDB']:
        process run_GenerateBrackenDB {
            label 'maxi'
            tag { "Generate Bracken DB" }
            
            """
            bracken-build -t ${task.cpus} -d ${kraken_db}
            """
        }
        break

        // MODE 5 - PERFORM READ QC
    case ['run.ReadQC']:
        process run_FastQC {
            label 'midi'
            tag { samples }
            publishDir "${qc_dir}", mode: 'copy', overwrite: true
            
            input:
            set samples, file(reads) from read_pairs
            
            output:
            set samples, file("*.{html,zip}") into out_fastqc

            """
            fastqc ${reads.findAll().join(' ') } \
                --threads ${task.cpus} \
                --noextract
            """
        }

        process run_MultiQC {
            label 'mini'
            tag { samples }
            publishDir "${qc_dir}", mode: 'copy', overwrite: true
            
            input:
            set samples, file(qc) from out_fastqc

            output:
            set samples, file('fastqc_multiqc') into out_multiqc
            
            """
            multiqc ${qc.findAll().join(' ') } --force -o fastqc_multiqc
            """
        }
        break

        // MODE 6 - PERFORM READ TRIMMING
    case ['run.ReadTrimming']:
        process run_Trimmomatic {
            label 'maxi'
            tag { sample }
            publishDir "${trim_dir}/${sample}", mode: 'copy', overwrite: true
            
            input:
            set sample, file(reads) from read_pairs
            
            output:
            set sample, file("${sample}_trimmed*.fastq.gz") into out_trimmed_reads

            """
            ln -s /opt/Trimmomatic-0.39/adapters/*.fa .
            if [[ ${stranded} == "paired-end" ]]
            then
                java -jar /opt/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
                    -threads ${task.cpus} \
                    ${reads.findAll().join(' ')} \
                    -baseout ${sample}_trimmed.fastq.gz \
                    ${trim_params}
            elif [[ ${stranded} == "single-end" ]]
            then
                java -jar /opt/Trimmomatic-0.39/trimmomatic-0.39.jar SE \
                    -threads ${task.cpus} \
                    ${reads.findAll().join(' ')} \
                    ${sample}_trimmed.fastq.gz \
                    ${trim_params}
            else
                :
            fi
            """
        }
        break

        // MODE 7 - PERFORM TAXONOMIC CLASSIFICATION
    case ['run.TaxonomicClassification']:
        process run_MetaPhlAn2 {
            label 'midi'
            tag { sample }
            publishDir "${tax_clas_dir}/${sample}", mode: 'copy', overwrite: true 
            
            input:
            set sample, file(reads) from read_pairs

            output:
            set sample, file("${sample}*") into out_metaphlan
            
            """
            /opt/metaphlan2/metaphlan2.py ${reads.findAll().join(' ')} \
                --bowtie2out ${sample}.bowtie2.bz2 \
                --nproc ${task.cpus} \
                --input_type fastq \
                -o ${sample}_MetaPhlAn2_microbes_list.tsv
            """
        }
        break

        // MODE 8 - STRAIN-COMPUTE
    case ['run.StrainComp']:
        
        read_pairs.into { read_pairs_1; read_pairs_2; read_pairs_3 }

        process run_Kraken2SC {
            label 'maxi'
            tag { sample }
            publishDir "${str_comp_dir}/${sample}", mode: 'copy', overwrite: true

            input:
            set sample, file(reads) from read_pairs_1
            
            output: 
            set sample, file("kraken_${sample}.tsv") into out_kraken 

            """
            if [[ ${stranded} == "paired-end" ]]
            then
                kraken2 --db $kraken_db --threads $task.cpus \
                    --paired \
                    --report kraken_${sample}.tsv \
                    --quick --memory-mapping \
                    ${reads.findAll().join(' ')}
            elif [[ ${stranded} == "single-end" ]]
            then
                kraken2 --db $kraken_db --threads $task.cpus \
                    --report kraken_${sample}.tsv \
                    --quick --memory-mapping \
                    ${reads.findAll().join(' ')}
            else
                :
            fi
            """
        }

        process run_BrackenSC {
            label 'maxi'
            tag { sample }
            publishDir "${str_comp_dir}/${sample}", mode: 'copy', overwrite: true 

            input: 
            set sample, file(table) from out_kraken
            
            output: 
            set sample, file("${sample}_bracken.tsv") into out_bracken_1, out_bracken_2

            """
            bracken -d $kraken_db -i ${table} \
                -o ${sample}_bracken.tsv \
                -r $readlen \
                -l $tax_level \
           """
        }

        process run_KronaSC {
            label 'mini'
            tag { sample }
            publishDir "${str_comp_dir}/${sample}", mode: 'copy', overwrite: true
            
            input:
            set sample, file(bhits) from out_bracken_1

            output: 
            set sample, file("${sample}_krona.html") into krona_out
            
            """
            ktImportTaxonomy -m 3 -s 0 -q 0 -t 5 \
                -tax ${taxonomy} \
                -i ${bhits} \
                -o ${sample}_krona.html
            """
        }
        
        out_bracken_2.collectFile() { item -> [ 'braken_file.txt', "${item.get(1)}" + '\n' ] }
            .set { bracken_files }
        
        process run_CollectResults {
            label 'mini'
            publishDir "${str_comp_dir}/collect", mode: 'copy', overwrite: true

            input:
            file(list) from bracken_files

            output: 
            file 'class_long.tsv' into collect_results_out

            """
            collate_results.sh ${tax_level} ${dataset_table} class_long.tsv ${list}
            """           
        }

        process run_CreateBarPlot {
            label 'mini'
            tag { "Create BarPlot" }
            publishDir "${str_comp_dir}/plots", mode: 'copy', overwrite: true

            input: 
            file(combined_bracken) from collect_results_out

            output: 
            file("relative_abundance_barplot.pdf") into barplot_out
            
            """
            composition_barplot.R ${combined_bracken} relative_abundance_barplot.pdf
            """
        }

        process run_SRST2 {
            label 'midi' 
            tag { sample }
            publishDir "${str_comp_dir}/srst2", mode: 'copy', overwrite: true

            input: 
            set val(sample), file(reads) from read_pairs_2 

            output: 
            set val(sample), file("*") into out_srst2
            
            """
            if [[ ${stranded} == "paired-end" ]]
            then 
                srst2 --input_pe ${reads.findAll().join(' ')} \
                    --forward _R1 --reverse _R2 \
                    --output srst2_${sample} \
                    --threads ${task.cpus} \
                    --gene_db ${annot_db}
            elif [[ ${stranded} == "single-end" ]]
            then
                srst2 --input_se ${reads.findAll().join(' ')} \
                    --output srst2_${sample} \
                    --gene_db ${annot_db}
            else
                :
            fi
            """
        }        
        break

        // MODE 9 - STRAIN SIFTER (IN DEV)
    case ['run.StrainSifter']:
        
        // process run_BWAAlign {
        //     label 'maxi'
        //     tag { sample }
        //     publisDir "${out_dir}/filtered_bam", mode: 'copy', overwrite: true
            
        //     input:
        //     set val(sample), file(reads) from read_pairs

        //     output:
        //     set val(sample), file("${sample}.filtered.bam") into out_bwa_align
            
        //     """
        //     bwa mem -t ${task.cpus} ${genome} ${reads.findAll().join(' ')} | \
        //         samtools view -b -q {params.qual} | \
        //         bamtools filter -tag 'NM:<={params.nm}' | \
        //         samtools sort --threads ${task.cpus} -o ${sample}.filtered.bam"
        //     """
        // }

        // process run_GenomeCov {
        //     label ''
        //     tag { sample }
        //     publisDir "${out_dir}/genome_cov", mode: 'copy', overwrite: true
            
        //     input:
        //     set val(sample), file(bam) from out_bwa_align

        //     output:
        //     set val(sample), file(${sample}_genome_cov.tsv) into out_genome_cov
            
        //     """
        //     bedtools genomecov -ibam ${bam} > ${sample}_genome_cov.tsv
        //     """
        // }

        // process run_

        // process run_StrainSifter {
        //     label 'maxi'
        //     tag { "Strain Sifter" }
        //     publishDir "${out_dir}/strainsifter", mode: 'copy', overwrite: true

        //     input:
        //     file strainSifter_config

        //     output: 
        //     file pdf  into StrainSifter_out
            
        //     """

        //     """
        // }
        break

        // MODE - VIRAL DETECTION FOR LONG READS
    case ['run.ViralDetectLong']:

        // read_pairs.subscribe { println it }

        read_pairs.into { read_pairs_1; read_pairs_2 }
        
        process run_Minimap2Decontaminate {
            label 'midi'
            tag { sample }
            publishDir "${vir_det_l_dir}/${sample}", mode: 'copy', overwrite: false, pattern: "${sample}_minimap2dc.sam"

            input:
            set sample, file(reads) from read_pairs_1

            output:
            set sample, file("${sample}_minimap2dc.sam"), file(reads) into out_minimap2dc_1, out_minimap2dc_2
            
            """
            minimap2 -ax ${read_type} ${params.decontam_db} -t ${task.cpus} ${reads.findAll().join(' ')} > ${sample}_minimap2dc.sam
            """
        }

        process run_RemoveReads {
            label 'mini'
            tag { sample }
            publishDir "${vir_det_l_dir}/${sample}", mode: 'copy', overwrite: true

            input:
            set sample, file(mmap), file(reads) from out_minimap2dc_1

            output:
            set sample, file("${sample}_cleaned.fastq") into clean_samples_1, clean_samples_2
            
            """
            remove-reads.py -s ${mmap} -i ${reads} > ${sample}_cleaned.fastq
            """
        }

        process run_Kraken2VDL {
            label 'maxi'
            tag { sample }
            publishDir "${vir_det_l_dir}/${sample}", mode: 'copy', overwrite: true

            input:
            set sample, file(reads) from clean_samples_1

            output:
            set sample, file("${sample}_kraken.tsv") into out_kraken_1, out_kraken_2, out_kraken_3

            """
            kraken2 --memory-mapping --quick \
                --db ${kraken_db} \
                --threads ${task.cpus} \
                --report ${sample}_kraken.tsv \
                ${reads.findAll().join(' ')}
            """
        }

        process run_Minimap2ViralCheck {
            label 'maxi'
            tag { sample } 
            publishDir "${vir_det_l_dir}/${sample}", mode: 'copy', overwrite: false

            input:
            set sample, file(reads) from clean_samples_2

            output:
            set sample, file("${sample}_minimap2viral.sam") into out_minimap2viral

            script:
            """
            minimap2 -ax ${read_type} ${viral_db} -t ${task.cpus} ${reads} > ${sample}_minimap2viral.sam
            """
            }

         out_kraken_1.join(out_minimap2dc_2)
            .map { it -> [ it[0], [ it[1], it[2], it[3] ] ] }
            .set { kraken_minimap2dc }        

        // process run_getIdentity {
        //     label 'mini'
        //     tag { sample }
        //     publishDir "${vir_det_l_dir}/${sample}", mode: 'copy', overwrite: false

        //     input:
        //     set sample, file(km2dc) from kraken_minimap2dc
            
        //     output:
        //     set sample, file("${sample}_identity-stats.txt") into out_identity_stats

        //     script:
        //     """
        //     get-identity.py \
        //         -k ${km2dc.get(0)} \
        //         -s ${km2dc.get(1)} \
        //         -a ${acc_2_tax} \
        //         -n ${params.names_dmp} > ${sample}_identity-stats.txt
        //     """
        // }
        
        process run_KronaVDL {
            label 'mini'
            tag { sample }
            publishDir "${vir_det_l_dir}/${sample}", mode: 'copy', overwrite: true

            input:
            set sample, file(khits) from out_kraken_2

            output:
            set sample, file("${sample}_krona.html") into out_krona_report
            
            """
            ktImportTaxonomy -m 3 -s 0 -q 0 -t 5 \
                -tax ${taxonomy} \
                -i ${khits} \
                -o ${sample}_krona.html
            """
        }

         out_kraken_3.join(out_krona_report)
            .map { it -> [ it[0], [ it[1], it[2] ] ] }
            .set { kraken_krona_reports }

        process run_FinalReportVDL {
            label 'mini'
            tag { sample }
            publishDir "${vir_det_l_dir}/${sample}", mode: 'copy', overwrite: true
            
            input:
            set sample, file(krakro)from kraken_krona_reports
            
            output:
            set sample, file("${sample}_*") into out_final_report
            
            """
            final-report-long-reads.py ${krakro.get(0)} ${krakro.get(1)} ${sample}_final-report.html
            """
        }
        break
        
        // MODE - VIRAL DETECT FOR SHORT READS
    case ['run.ViralDetectShort']:        

        process run_RemoveHost2 {
            label 'midi'
            tag { sample }
            publishDir "${vir_det_s_dir}/${sample}", mode: 'copy', overwrite: true
            
            input:
            set sample, file(reads) from read_pairs
            
            output:
            set sample, file("*") into out_fastq_screen
            
            """
            cat <<EOF > myConf.conf
            DATABASE`printf '\t'`VIRALDB`printf '\t'`${viral_db.getParent()}/${viral_db.getBaseName()}
            EOF

            fastq_screen --conf myConf.conf --nohits --aligner bowtie2 --threads ${task.cpus} ${reads.findAll().join(' ') }
            """
        }

        process run_PairFQ {
            label 'midi'
            tag { sample }
            publishDir "${vir_det_s_dir}/${sample}", mode: 'copy', overwrite: true

            input:
            set sample, file(reads) from out_fastq_screen
            
            output:
            set sample, file ("*.{paired_R1.fastq.gz,paired_R2.fastq.gz}") into paired_kraken1, paired_kraken2, paired_bowtie1, paired_bowtie2

            """
            pairfq addinfo \
                -i ${reads.get(0)} \
                -o ${sample}.addinfo_R1.fastq \
                -p 1

            pairfq addinfo \
                -i ${reads.get(1)} \
                -o ${sample}.addinfo_R2.fastq \
                -p 2

            pairfq makepairs \
                -f ${sample}.addinfo_R1.fastq \
                -r ${sample}.addinfo_R2.fastq \
                -fp ${sample}.paired_R1.fastq \
                -rp ${sample}.paired_R2.fastq \
                -fs ${sample}.singleton_R1.fastq \
                -rs ${sample}.singleton_R2.fastq \
                -c gzip

            """
        }

        process run_Bowtie2VDS {
            label 'maxi'
            tag { sample }
            publishDir "${vir_det_s_dir}/${sample}", mode: 'copy', overwrite: true
            
            input:
            set sample, file(reads) from paired_bowtie1

            output:
            set sample, file("${sample}_mapped.sam") into reads_mapped
            
            """
            bowtie2 -x ${viral_db.getParent()}/${viral_db.getBaseName()} \
                -1 ${reads.get(0)} -2 ${reads.get(1)} \
                --threads ${task.cpus} \
                -S ${sample}_mapped.sam --no-unal
            """
        }

        process run_MappingsStats {
            label 'midi'
            tag { sample }
            publishDir "${vir_det_s_dir}/${sample}", mode: 'copy', overwrite: false

            input:
            set sample, file(sam) from reads_mapped

            output:
            set sample, file("${sample}_mapping_stats.txt") into reads_mapped_stats

            """
            samtools view -S -b -f 12 -F 256 ${sam} | samtools sort - > ${sample}.sorted.bam
            samtools index ${sample}.sorted.bam
            samtools idxstats ${sample}.sorted.bam > ${sample}_mapping_stats.txt
            """
        }

        process run_Kraken2VDS {
            label 'maxi'
            tag { sample }
            publishDir "${vir_det_s_dir}/${sample}", mode: 'copy', overwrite: false

            input:
            set sample, file(reads) from paired_kraken1

            output:
            set sample, file("${sample}_kraken-hits.tsv") into kraken_hits, kraken_report

            """
            kraken2 --memory-mapping --quick \
                --db ${kraken_db} \
                --threads ${task.cpus} \
                --report ${sample}_kraken-hits.tsv \
                --paired ${reads.findAll().join(' ')}
            """
        }

        process run_KronaVDS {
            label 'mini'
            tag { sample }
            publishDir "${vir_det_s_dir}/${sample}", mode: 'copy', overwrite: false
            
            input:
            set sample, file(khits) from kraken_hits

            output:
            set sample, file("${sample}_krona.html")  into krona_charts

            """
            ktImportTaxonomy -m 3 -s 0 -q 0 -t 5 \
                -tax ${taxonomy} \
                -i ${khits} \
                -o ${sample}_krona.html
            """
        }


        kraken_report.join(reads_mapped_stats)
            .join(krona_charts)
            .map { it -> [ it[0], [ it[1], it[2], it[3] ] ] }
            .set { combined_reports }
 
        process run_FinalReportVDS {
            label 'mini'
            tag { sample }
            publishDir "${vir_det_s_dir}/${sample}", mode: 'copy', overwrite: false
            
            input:
            set sample, file(reports) from combined_reports

            output:
            set sample, file ("final-report.html") into out_final_report_vds

            """
            grep "^>" ${viral_db} > search_pattern.txt
            final-report-short-reads.py ${reports.get(0)} ${reports.get(1)} ${reports.get(2)} final-report.html search_pattern.txt
            """
        }
        
        break
}

workflow.onComplete {
    println "\n${line}"
    println "#".multiply(48 - ("${summary}".size() / 2 )) + "  ${summary}  " + "#".multiply(48 - ("${summary}".size() / 2 ))    
    println "${line}\n"
    println "Execution command   : ${workflow.commandLine}"
    println "Execution name      : ${workflow.runName}"
    println "Workflow start      : ${workflow.start}"
    println "Workflow end        : ${workflow.complete}"
    println "Workflow duration   : ${workflow.duration}"
    println "Workflow completed? : ${workflow.success}"
    println "Work directory      : ${workflow.workDir}"
    println "Project directory   : ${workflow.projectDir}"
    println "Execution directory : ${workflow.launchDir}"
    println "Configuration files : ${workflow.configFiles}"
    println "Workflow containers : ${workflow.container}"
    println "exit status         : ${workflow.exitStatus}"
    println "Error report        : ${workflow.errorReport ?: '-'}"
    println "${line}\n"
    println "\n"
}
