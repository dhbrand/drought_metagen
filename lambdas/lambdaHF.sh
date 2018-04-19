#!/bin/bash
#
#
#SBATCH -J HF
#SBATCH -A iPlant-Collabs
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -t 12:00:00
#SBATCH -p skx-normal
#SBATCH -o job_%j_%N.out

Rscript --vanilla --verbose ./lambdaHF.R > outputHF.Rout
