#!/usr/bin/env python3.6
# Takes Kraken2 and Minimap2 results and check what the percentages are of the top most hits. It also tries to check if the top hits in the Krakan2 and Minimap2 results matches.
# Input files (1):
#    kraken_file    minimap2 file (sam file)    NCBI's accession to taxonomic id mapping file   NCBI's taxonomic to naming mapping file
# Stdout (2)
#   Stdout containing taxonomies in both Kraken2 and Minimap2 with highest hits. Also the species / strain match between the two reports.
import sys
import re
import string
from optparse import OptionParser
import subprocess
import tempfile
import csv
import pysam
import re

def main():
    usage = "usage: %prog -k kraken2_file -s sam_file -a nucl_gb_accession2taxid_file -n names_dmp_file"
    parser = OptionParser(usage=usage)
    parser.add_option("-k", "--kraken2", dest="kraken2_file", help="Kraken2 TSV input file.")
    parser.add_option("-s", "--sam", dest="sam_file", help="SAM input file.")
    parser.add_option("-a", "--accession2taxid", dest="nucl_gb_accession2taxid_file", help="NCBIs GenBank accession to taxonomic id mapping file.")
    parser.add_option("-n", "--names_dmp", dest="names_dmp_file", help="NCBIs taxonomic id to naming mapping file.")
    (options, args) = parser.parse_args()

    if not options.kraken2_file:
        print ("Please specify the kraken2 TSV input file. (-k kraken2_file)")
        return - 1
    if not options.sam_file:
        print ("Please specify the SAM input file (-s sam_file)")
        return - 2
    if not options.nucl_gb_accession2taxid_file:
        print ("Please specify the NCBIs GenBank accession to taxonomic id mapping file (-a nucl_gb_accession2taxid_file)")
        return - 3
    if not options.names_dmp_file:
        print ("Please specify the NCBIs taxonomic id to naming mapping file (-n names_dmp_file)")
        return - 4
    if (len(args) > 0):
        print ("Too many input arguments")
        return - 5

    kraken2_file = options.kraken2_file
    sam_file = options.sam_file
    a2t_file = options.nucl_gb_accession2taxid_file
    t2n_file  = options.names_dmp_file

    # Get results from kraken2 report. Need to make an external Linux call.
    tmp_file = tempfile.mkstemp()[1]
    p = subprocess.Popen(["kraken-biom",kraken2_file, "--fmt","tsv","-o",tmp_file], shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (stdout, stderr) = p.communicate()
    ## Now we get the kraken2 results in a more parsable format
    kraken_file_ids_only_tsv = csv.reader(open(tmp_file), delimiter='\t')
    next(kraken_file_ids_only_tsv) # Skip first two lines, those are headers
    next(kraken_file_ids_only_tsv)
    kraken2_results = {}
    for row in kraken_file_ids_only_tsv:
        if row[0] in kraken2_results:
            kraken2_results[row[0]] = kraken2_results[row[0]] + row[1]
        else:
            kraken2_results[row[0]] = row[1]
    ## Now we get hit statistics
    kraken2_max_hit = 0.0
    kraken2_max_taxonomy = ""
    kraken2_total_hits = 0.0
    for kraken2_taxonomy in kraken2_results:
        if (float(kraken2_results[kraken2_taxonomy]) >= kraken2_max_hit):
            kraken2_max_hit = float(kraken2_results[kraken2_taxonomy])
            kraken2_max_taxonomy = kraken2_taxonomy
        kraken2_total_hits = kraken2_total_hits + float(kraken2_results[kraken2_taxonomy])
    kraken2_percentage = 100*(float(kraken2_max_hit) / kraken2_total_hits)
    print ("Kraken2 max hit is: " + kraken2_max_taxonomy + " (taxonomy id) with a " + str(kraken2_percentage) + "% of overall hits")
    ## Now we try do find the hit details in the NCBI databases
    t2n_file_tsv = csv.reader(open(t2n_file), delimiter='\t')
    kraken2_taxonomy_name_hit = ""
    for row in t2n_file_tsv:
        if (row[0] == kraken2_max_taxonomy) and (row[6] == "scientific name"):
            kraken2_taxonomy_name_hit = row[2]

    # Now work on minimap2/SAM hits
    sam_file_fd = pysam.AlignmentFile(sam_file, "r")
    minimap2_results = {}
    for read in sam_file_fd.fetch():
        if read.reference_name in minimap2_results:
            minimap2_results[read.reference_name] = minimap2_results[read.reference_name] + 1
        else:
            if read.reference_name is not None:
                minimap2_results[read.reference_name] = 1
    ## Now we get hit statistics
    minimap2_max_hit = 0.0
    minimap2_max_accession = ""
    minimap2_total_hits = 0.0
    for minimap2_accesssion in minimap2_results:
        if (float(minimap2_results[minimap2_accesssion]) >= minimap2_max_hit):
            minimap2_max_hit = float(minimap2_results[minimap2_accesssion])
            minimap2_max_accession = minimap2_accesssion
        minimap2_total_hits = minimap2_total_hits + float(minimap2_results[minimap2_accesssion])
    minimap2_percentage = 100*(float(minimap2_max_hit) / minimap2_total_hits)
    print ("Minimap2 max hit is: " + minimap2_max_accession + " (GenBank/Refseq accession id) with a " + str(minimap2_percentage) + "% of overall hits")
    ## Now we try do find the hit details in the NCBI's accession to taxonomy databases
    ### First we need to find the correct refseq id from the complete bacterial / viral dbs dowloaded from NCBI
    m = re.match(r'.*ref\|(.*)\|',minimap2_max_accession)
    tmp_accession = m.group(1)
    minimap2_max_accession = tmp_accession
    a2t_file_tsv = csv.reader(open(a2t_file), delimiter='\t')
    minimap2_max_taxonomy = ""
    for row in a2t_file_tsv:
        if row[1] == minimap2_max_accession:
            minimap2_max_taxonomy = row[2] # We assume we will only get one hit here
    ## Now we try do find the hit details in the NCBI taxonomy to name databases
    t2n_file_tsv = csv.reader(open(t2n_file), delimiter='\t')
    minimap2_taxonomy_name_hit = ""
    for row in t2n_file_tsv:
        if (row[0] == minimap2_max_taxonomy) and (row[6] == "scientific name"):
            minimap2_taxonomy_name_hit = row[2]

    # Now let us check matches
    ## First check if we have a match in the kraken2 string with that of the minimap2 string. It might be that we have the species in the one and strain in the other
    ## Lets convert to lower case migth help with searching
    minimap2_taxonomy_name_hit = minimap2_taxonomy_name_hit.lower()
    kraken2_taxonomy_name_hit = kraken2_taxonomy_name_hit.lower()

    match = re.findall(minimap2_taxonomy_name_hit, kraken2_taxonomy_name_hit)
    if (match):
        print ("The match between kraken2 and minimap2 is " + match[0])
    ## Now we check the other way around
    match = re.findall(kraken2_taxonomy_name_hit, minimap2_taxonomy_name_hit)
    if (match):
        print ("The match between kraken2 and minimap2 is " + match[0])

    print ("Kraken2 top name hit is " + kraken2_taxonomy_name_hit)
    print ("Minimap2 top name hit is " + minimap2_taxonomy_name_hit)

    return 0

if __name__ == "__main__":
    sys.exit(main())
