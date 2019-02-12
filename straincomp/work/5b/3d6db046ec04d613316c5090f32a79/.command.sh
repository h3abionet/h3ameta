#!/usr/bin/env bash
#note that the code below can be written in any language, so long as the interpreter is named
#in the above shebang line.
#also note: variables from this or higher scopes can be named inside strings with a dollar sign, as below.
#when part of a larger name, as in the tsv output filename below, the variable must be enclosed with {}
# the final test4.fq names the input.
kraken2 --db mini8GB/ --threads 2 	--report test4_kraken.tsv 	--quick --memory-mapping 	test4.fq
