#!/bin/bash
#$ -N sns_layer23_de_seurat_mast
#$ -j y
#$ -o ../../../../out/writeup8c/qsub/
#$ -l m_mem_free=50G

Rscript --no-save writeup8c_sns_layer23_de_seurat_mast.R
