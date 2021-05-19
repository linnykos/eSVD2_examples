#!/bin/bash
#$ -N fasttopics
#$ -j y
#$ -o ../../../../out/writeup6/qsub/

Rscript --no-save writeup6_citeseq_bm_fasttopics.R
