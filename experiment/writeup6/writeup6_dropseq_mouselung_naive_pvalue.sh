#!/bin/bash
#$ -N lung_naive_pvalue
#$ -j y
#$ -o ../../../../out/writeup6/qsub/
#$ -l m_mem_free=12G

Rscript --no-save writeup6_dropseq_mouselung_naive_pvalue.R
