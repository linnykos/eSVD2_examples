#!/bin/bash
#$ -N lung_glmpca_nb
#$ -j y
#$ -o ../../../../out/writeup6/qsub/

Rscript --no-save writeup6_dropseq_mouselung_glmpca_nb.R
