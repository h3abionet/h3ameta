#!/bin/bash
output=/home/abdulrahman/h3ameta/preprocessing/Abdul/Trim
input="/home/abdulrahman/h3ameta/preprocessing/Abdul/quality_control"
for i in $input/*_1.fq.gz;
do
withpath="${i}" filename=${withpath##*/}
base="${filename%*_*.fq.gz}"
sample_name=`echo "${base}" | awk -F ".fastq.gz" '{print $1}'`
java -jar /home/abdulrahman/h3ameta/preprocessing/Abdul/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 6 -trimlog $output/"${base}".log $input/"${base}"_1.fq.gz $input/"${base}"_2.fq.gz $output/"${base}"_R1.trimmed_PE.fastq $output/"${base}"_R1.trimmed_SE.fastq $output/"${base}"_2.trimmed_PE.fastq $output/"${base}"_2.trimmed_SE.fastq LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20



done
