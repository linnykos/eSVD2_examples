#!/bin/bash
#$ -N brain_glmpca_nb
#$ -j y
#$ -o ../../../../out/writeup6/qsub/
#$ -l m_mem_free=12G

Rscript --no-save writeup6_smartseq_mousebrain_glmpca_nb.R
