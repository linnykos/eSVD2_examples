#!/bin/bash
#$ -N bm_glmpca_poisson
#$ -j y
#$ -o ../../../../out/writeup6/qsub/

Rscript --no-save writeup6_citeseq_bm_glmpca_poisson.R
