#!/bin/bash
#$ -N habermann_T_sctransform
#$ -j y
#$ -o ../../../out/main/qsub/
#$ -l m_mem_free=100G

Rscript --no-save habermann_T_sctransform.R
