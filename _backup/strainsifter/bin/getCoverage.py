#!/usr/bin/env python3
# Fiona Tamburini
# Jan 26 2018
# Modified for NextFlow by SH, March 2019

# calculate average coverage and % bases covered above min threshold from bedtools genomecov output
# usage python3 getCoverage.py 5 sample < depth.tsv

import sys


# minimum coverage threshold -- report percentage of bases covered at or
# beyond this depth


sample_file = sys.argv[1]
minCvg      = float(sys.argv[2])
out_file    = sys.argv[3]



totalBases = 0
coveredBases = 0
weightedAvg = 0
with open(sample_file, 'r') as sample:
    for line in sample:
        if line.startswith("genome"):
            chr, depth, numBases, size, fraction = line.rstrip('\n').split('\t')

            depth = int(depth)
            numBases = int(numBases)

            totalBases += numBases
            weightedAvg += depth * numBases

            if depth >= minCvg:
                coveredBases += numBases
avgCvg = 0
percCovered = 0

if totalBases > 0:
    avgCvg = float(weightedAvg) / float(totalBases)
    percCovered = float(coveredBases) / float(totalBases)


with open(out_file, 'w') as out:
    print("\t".join([str(round(avgCvg, 2)), str(round(percCovered, 2))]), file = out)
