#!/bin/bash
#$ -N habermann_T
#$ -j y
#$ -o ../../../../out/Writeup11g/qsub/
#$ -l m_mem_free=100G

Rscript --no-save Writeup11g_habermann_T_esvd.R
