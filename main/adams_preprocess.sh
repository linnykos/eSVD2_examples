#!/bin/bash
#$ -N adams_preprocess
#$ -j y
#$ -o ../../../out/main/qsub/
#$ -l m_mem_free=200G

Rscript --no-save adams_preprocess.R
