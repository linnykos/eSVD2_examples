#!/bin/bash
#$ -N retinal_esvd_glmgampoi
#$ -j y
#$ -o ../../../../out/writeup8/qsub/
#$ -l m_mem_free=50G

Rscript --no-save writeup8_10x_mouseretinal_esvd_glmgampoi.R
