#!/bin/bash
#PBS -N sqanti
#PBS -q submit
#PBS -l nodes=1:ppn=2,mem=20gb,walltime=24:00:00
#PBS -M dong.x@wehi.edu.au
#PBS -m abe
#PBS -j oe

module load anaconda3
source activate
conda activate SQANTI
module load R
module load perl
module load gmap-gsnap

cd /wehisan/home/allstaff/d/dong.x/analysis/smchd1/long/scFLT

# python /wehisan/home/allstaff/d/dong.x/Programs/SQANTI/sqanti_qc.py \


python /wehisan/home/allstaff/d/dong.x/Programs/SQANTI/sqanti_qc.py \
  run1/isoform_annotated.gff3 \
  /wehisan/home/allstaff/d/dong.x/annotation/Mouse/gencode.vM24.annotation.gtf \
  /wehisan/home/allstaff/d/dong.x/annotation/Mouse/mm10/mm10.fa -g
