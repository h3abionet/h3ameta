#!/usr/bin/env bash

       mkdir sample2_2.fq

fastqc -t 4 sample2_2.fq.gz -o sample2_2.fq
