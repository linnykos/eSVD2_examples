#!/bin/bash
#$ -N cortical_glmpca_p
#$ -j y
#$ -o ../../../../out/writeup6b/qsub/
#$ -l m_mem_free=50G

Rscript --no-save writeup6b_dropseq_humancortical_glmpca_poisson.R
