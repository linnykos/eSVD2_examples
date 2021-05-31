#!/bin/bash
#$ -N bm_naive_pvalue
#$ -j y
#$ -o ../../../../out/writeup6/qsub/

Rscript --no-save writeup6_citeseq_bm_naive_pvalue.R
