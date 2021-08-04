#!/bin/bash
#$ -N retinal_esvd
#$ -j y
#$ -o ../../../../out/writeup6b/qsub/
#$ -l m_mem_free=50G

Rscript --no-save writeup6b_10x_mouseretinal_esvd.R
