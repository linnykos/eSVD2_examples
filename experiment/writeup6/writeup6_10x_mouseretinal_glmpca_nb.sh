#!/bin/bash
#$ -N retinal_glmpca_nb
#$ -j y
#$ -o ../../../../out/writeup6/qsub/
#$ -l m_mem_free=12G

Rscript --no-save writeup6_10x_mouseretinal_glmpca_nb.R
