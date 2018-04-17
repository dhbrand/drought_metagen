#!/bin/bash
#SBATCH --mem=72g
#SBATCH -N 1
#SBATCH -c 16
#SBATCH -t 263:00:00

MICROBIOME=/nas/longleaf/home/roachjm/MICROBIOME
source $MICROBIOME/microbiome.source

WGS=/nas/longleaf/home/roachjm/MICROBIOME/scripts
WGS_DATA=/nas/depts/001/roachjm/WGS/pathoscopedb/

echo $1
python $WGS/pathoscope2/pathoscope2.py MAP -numThreads 16 -U FILTER2/$1.fastq -targetRefFiles $WGS_DATA/JMR_ti.fa -outDir PATHOSCOPE2 -outAlign $1.sam -expTag $1 -indexDir $WGS_DATA

rm PATHOSCOPE2/$1-JMR_ti_?.sam

python $WGS/pathoscope2/pathoscope2.py ID -alignFile PATHOSCOPE2/$1.sam -fileType sam -outDir PATHOSCOPE2 -expTag $1

rm PATHOSCOPE2/$1.sam

awk -f $WGS/pathoscope-scripts/addTaxInfo2.awk $WGS_DATA/taxa_names.txt PATHOSCOPE2/$1-sam-report.tsv > PATHOSCOPE2/$1-sam-report-FINAL.tsv
