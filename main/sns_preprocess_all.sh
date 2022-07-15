#!/bin/bash
#$ -N sns_preprocess
#$ -j y
#$ -o ../../../out/main/qsub/
#$ -l m_mem_free=200G

Rscript --no-save sns_preprocess_all.R
