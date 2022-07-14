#!/bin/bash
#$ -N sns_layer23_2
#$ -j y
#$ -o ../../../../out/Writeup11g/qsub/
#$ -l m_mem_free=100G

Rscript --no-save Writeup11g_sns_layer23_esvd2.R
