#!/usr/bin/env bash

       mkdir sample1_1.fq

fastqc -t 4 sample1_1.fq.gz -o sample1_1.fq
