#!/usr/bin/env bash

       mkdir sample2_1.fq

fastqc -t 4 sample2_1.fq.gz -o sample2_1.fq
