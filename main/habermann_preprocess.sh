#!/bin/bash
#$ -N habermann_preprocess
#$ -j y
#$ -o ../../../out/main/qsub/
#$ -l m_mem_free=200G

Rscript --no-save habermann_preprocess.R
