#!/bin/bash
#$ -N cortical_esvd_glmgampoi
#$ -j y
#$ -o ../../../../out/writeup8/qsub/
#$ -l m_mem_free=50G

Rscript --no-save writeup8_dropseq_humancortical_esvd_nb_glmgampoi.R