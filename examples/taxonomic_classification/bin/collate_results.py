#!/usr/bin/env python

import sys
import os
import re


tax_level = sys.argv[1]
datasets_table_f = sys.argv[2]
output_f = sys.argv[3]
class_filenames = sys.argv[4:]
output = []

sample_condtp = {}
with open(datasets_table_f, 'r') as asmf:
	basedir = os.path.dirname(datasets_table_f)
	for l in asmf.readlines():
		s = l.strip().split("\t")
		if len(s) == 1 or s[0] == 'Sample':
			continue

		sample = s[0]
		timepoint = s[1]
		condition = s[2]

		sample_condtp[sample] = "\t".join([condition, timepoint])

for inf in class_filenames:
	with open(inf, 'r') as i:
		for l in i.readlines()[1:]:
			s = re.sub(r' +', "\t", l.strip()).split("\t")
			sample = os.path.splitext(os.path.basename(inf))[0].split('.')[0]
			sample = sample.replace('_kraken_bracken','')

			percent = s[0]
			reads = s[1]
			direct_reads = s[2]
			taxonomy = s[3]
			taxid = '_'.join(s[6:])
			newline = [sample, sample_condtp[sample], percent, reads, direct_reads, taxid, taxonomy]
			if taxonomy == tax_level and int(reads) > 1:
				output.append(newline)

with open(output_f, 'w') as out:
	out.write("\n".join(["\t".join(o) for o in output]))

#sys.exit('crashbangboom')
