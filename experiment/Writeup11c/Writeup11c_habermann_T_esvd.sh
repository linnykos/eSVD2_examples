#!/bin/bash
#$ -N habermann_T
#$ -j y
#$ -o ../../../../out/Writeup11c/qsub/
#$ -l m_mem_free=100G

Rscript --no-save Writeup11c_habermann_T_esvd.R
