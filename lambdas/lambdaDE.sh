#!/bin/bash
#
#
#SBATCH -J DE
#SBATCH -A iPlant-Collabs
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -t 12:00:00
#SBATCH -p skx-normal
#SBATCH -o job_%j_%N.out

Rscript --vanilla --verbose ./lambdaDE.R > outputDE.Rout
