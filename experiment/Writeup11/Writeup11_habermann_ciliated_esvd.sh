#!/bin/bash
#$ -N habermann_ciliated_esvd
#$ -j y
#$ -o ../../../../out/Writeup11/qsub/
#$ -l m_mem_free=100G

Rscript --no-save Writeup11_habermann_ciliated_esvd.R
