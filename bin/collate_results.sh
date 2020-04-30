#!/usr/bin/env bash

tax=$1
table=$2
out=$3
list=$4

while IFS=$'\t' read -r sample timepoint condition read
do
    if [[ $sample =~ ^Sample ]]
    then
        :
    else
        while read -r braken_file
        do
            if [[ $sample == `basename $braken_file | sed 's/_bracken.tsv//' ` ]]
            then 
                while IFS=$'\t' read -r name taxonomy_id taxonomy_lvl kraken_assigned_reads added_reads new_est_reads fraction_total_reads
                do
                    if [[ $name == "name" ]]
                    then
                        :
                    elif [[ $taxonomy_lvl == "$tax" && $kraken_assigned_reads > 1 ]]
                    then
                        name=`sed 's/ \{1,\}/_/g' <<< "$name"`
                        echo -e $sample"\t"$condition"\t"$timepoint"\t"$fraction_total_reads"\t"$kraken_assigned_reads"\t"$new_est_reads"\t"$name"\t"$taxonomy_lvl
                    else
                        :
                    fi
                done < $braken_file
            else
                :
            fi
        done < $list
    fi
done < $table > $out