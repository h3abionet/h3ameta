#!/usr/bin/env python3

import sys
import itertools
import subprocess
import re

# input and output files from snakemake
in_files = snakemake.input
out_file = snakemake.output[0]

cmd = " ".join(["cat", in_files[0], "| wc -l"])
totalBases = subprocess.check_output(cmd, stdin=subprocess.PIPE, shell=True ).decode('ascii').rstrip('\n')
totalBases = int(totalBases) - 1

# for every pairwise combination of files, check SNV distance
with open(out_file, "w") as sys.stdout:

    # print header
    print('\t'.join(["Sample1", "Sample2", "SNVs", "BasesCompared", "TotalBases"]))

    # get all pairwise combinations of input files
    for element in itertools.combinations(in_files, 2):
        file1, file2 = element

        cmd1 = " ".join(["paste", file1, file2, "| sed '1d' | grep -v N | wc -l"])
        totalPos = subprocess.check_output(cmd1, stdin=subprocess.PIPE, shell=True ).decode('ascii').rstrip('\n')

        cmd2 = " ".join(["paste", file1, file2, "| sed '1d' | grep -v N | awk '$1 != $2 {print $0}' | wc -l"])
        diffPos = subprocess.check_output(cmd2, stdin=subprocess.PIPE, shell=True ).decode('ascii').rstrip('\n')

        fname1 = re.findall('consensus/(.+)\.', file1)[0]
        fname2 = re.findall('consensus/(.+)\.', file2)[0]

        # dist = subprocess.check_output(['./count_snvs.sh', file1, file2]).decode('ascii').rstrip('\n')
        print('\t'.join([fname1, fname2, str(diffPos), str(totalPos), str(totalBases)]))
