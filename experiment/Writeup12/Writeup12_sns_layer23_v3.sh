#!/bin/bash
#$ -N sns_layer23_v3
#$ -j y
#$ -o ../../../../out/Writeup12/qsub/
#$ -l m_mem_free=15G

Rscript --no-save Writeup12_sns_layer23_v3.R
