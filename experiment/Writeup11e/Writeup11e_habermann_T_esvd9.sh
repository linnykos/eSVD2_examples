#!/bin/bash
#$ -N habermann_T_9
#$ -j y
#$ -o ../../../../out/Writeup11e/qsub/
#$ -l m_mem_free=100G

Rscript --no-save Writeup11e_habermann_T_esvd9.R
