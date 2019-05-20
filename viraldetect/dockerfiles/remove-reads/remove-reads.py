#!/usr/bin/env python3.6
# Find queries in SAM file that did not have matches to the reference. These reads are then selected from the original read set and those that did match are removed. Can be used as a decontamination tool.
# Input files (1):
#    minimap2 file (sam file)   fastq/fastq.gz input file
# Sdtoutfiles (2)
#   Fastq file with reads that did not match.
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
    usage = "usage: %prog -s sam_file -i fastq_file_in"
    parser = OptionParser(usage=usage)
    parser.add_option("-s", "--sam", dest="sam_file", help="SAM input file.")
    parser.add_option("-i", "--fastq", dest="fastq_file", help="Fastq/Fastq.gz input file.")

    (options, args) = parser.parse_args()

    if not options.sam_file:
        print ("Please specify the SAM input file (-s sam_file)")
        return - 1
    if not options.fastq_file_in:
        print ("Please specify the Fastq/Fastq.gz input file (-i fastq_file_in)")
        return - 2
    if (len(args) > 0):
        print ("Too many input arguments")
        return - 3

    sam_file = options.sam_file
    fastq_file_in = options.fastq_file_in

    # Get hits that did not hit a reference
    tmp_file = tempfile.mkstemp()[1]
    tmp_file_fd = open(tmp_file, "w")
    sam_file_fd = pysam.AlignmentFile(sam_file, "r")
    for read in sam_file_fd.fetch():
        if read.reference_name is None:
            tmp_file_fd.write ("%s\n" % read.query_name)
    tmp_file_fd.close()
    # Now we select only those no hits for our downstream processing
    p = subprocess.Popen(["seqtk","subseq",fastq_file, tmp_file], shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (stdout, stderr) = p.communicate()
    print (stdout.decode('ascii'))

    return 0

if __name__ == "__main__":
    sys.exit(main())
