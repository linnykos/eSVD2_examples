#!/bin/bash
#$ -N lung_glmpca_nb
#$ -j y
#$ -o ../../../../out/writeup6/qsub/
#$ -l m_mem_free=12G

Rscript --no-save writeup6_dropseq_mouselung_glmpca_nb_pvalue.R
