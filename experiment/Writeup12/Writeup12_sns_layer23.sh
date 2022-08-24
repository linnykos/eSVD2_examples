#!/bin/bash
#$ -N sns_layer23
#$ -j y
#$ -o ../../../../out/Writeup12/qsub/
#$ -l m_mem_free=50G

Rscript --no-save Writeup12_sns_layer23.R
