#!/bin/bash
#$ -N habermann_T
#$ -j y
#$ -o ../../../../out/Writeup11d/qsub/
#$ -l m_mem_free=100G

Rscript --no-save Writeup11d_habermann_T_esvd.R
