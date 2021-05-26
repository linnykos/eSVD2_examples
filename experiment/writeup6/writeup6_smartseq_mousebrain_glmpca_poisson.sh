#!/bin/bash
#$ -N brain_glmpca_poisson
#$ -j y
#$ -o ../../../../out/writeup6/qsub/
#$ -l m_mem_free=20G

Rscript --no-save writeup6_smartseq_mousebrain_glmpca_poisson.R
