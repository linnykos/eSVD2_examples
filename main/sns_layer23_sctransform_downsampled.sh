#!/bin/bash
#$ -N sns_layer23_sctransform_downsampled
#$ -j y
#$ -o ../../../out/main/qsub/
#$ -l m_mem_free=15G

Rscript --no-save sns_layer23_sctransform_downsampled.R
