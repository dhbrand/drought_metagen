#!/bin/bash
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=32g
#SBATCH -t 263:00:00

MICROBIOME=/nas/longleaf/home/roachjm/MICROBIOME
source $MICROBIOME/microbiome.source
module add trim_galore

ZM_REF=ZM/zm

echo $1

fastqc -t 8 -f fastq FASTQs/$1_R1.fastq -o QC0
fastqc -t 8 -f fastq FASTQs/$1_R2.fastq -o QC0

trim_galore --paired --illumina FASTQs/$1_R1.fastq FASTQs/$1_R2.fastq -o FILTER0

fastqc -t 8 -f fastq FILTER0/$1_R1_val_1.fq -o QC1
fastqc -t 8 -f fastq FILTER0/$1_R2_val_2.fq -o QC1

bowtie2 -p 8 --un-conc FILTER1/$1.NONREF.fastq -1 FILTER0/$1_R1_val_1.fq -2 FILTER0/$1_R2_val_2.fq -S /dev/null -x $ZM_REF

vsearch --threads 8 --fastq_mergepairs FILTER1/$1.NONREF.1.fastq --reverse FILTER1/$1.NONREF.2.fastq --fastqout JOINs/$1_JOIN.fastq --fastq_minovlen 60 --fastq_maxee 5.0 --fastqout_notmerged_fwd JOINs/$1_JOIN_R1.fastq --fastqout_notmerged_rev JOINs/$1_JOIN_R2.fastq

cat JOINs/$1_JOIN_R1.fastq JOINs/$1_JOIN_R2.fastq JOINs/$1_JOIN.fastq > COMBINE/$1.fastq

bowtie2 -p 8 --un FILTER2/$1.fastq -U COMBINE/$1.fastq -S /dev/null -x $ZM_REF

sbatch -o $1.pathoscope2.%A.out scripts/pathoscope2.sh $1

humann2 --threads 8 --input FILTER2/$1.fastq --output HUMANN2/$1 --nucleotide-database /nas/depts/001/roachjm/WGS/chocophlan/chocophlan --protein-database /nas/depts/001/roachjm/WGS/diamond_db/uniref --metaphlan /nas/longleaf/home/roachjm/MICROBIOME/scripts/metaphlan2 --metaphlan-options "--mpa_pkl /nas/depts/001/roachjm/WGS/metaphlan2db/db_v20/mpa_v20_m200.pkl --bowtie2db /nas/depts/001/roachjm/WGS/metaphlan2db/db_v20/mpa_v20_m200"

humann2_renorm_table --input HUMANN2/$1/$1_genefamilies.tsv --output HUMANN2/$1/$1_genefamilies-cpm.tsv --units cpm --update-snames

cp HUMANN2/$1/$1_humann2_temp/$1_metaphlan_bugs_list.tsv HUMANN2/$1/$1_metaphlan2.tsv
cp HUMANN2/$1/$1_humann2_temp/$1_metaphlan_bowtie2.txt HUMANN2/$1/$1_metaphlan_bowtie2.txt
#rm -rf HUMANN2/$1/$1_humann2_temp
