#!/usr/bin/env bash

       mkdir sample1_2.fq

fastqc -t 4 sample1_2.fq.gz -o sample1_2.fq
