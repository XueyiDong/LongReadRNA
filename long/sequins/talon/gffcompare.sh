#!/bin/bash
#PBS -N gffcompare
#PBS -q submit
#PBS -l nodes=1:ppn=4,mem=60gb,walltime=1000:00:00
#PBS -M dong.x@wehi.edu.au
#PBS -m abe
#PBS -j oe

cd /wehisan/home/allstaff/d/dong.x/analysis/smchd1/long/sequins/talon

module load anaconda3
#conda activate flair

UNIXHOME=/wehisan/home/allstaff/d/dong.x

$UNIXHOME/Programs/gffcompare -r /stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/Mike_seqin/annotations/rnasequin_annotation_2.4.gtf -R -o /wehisan/home/allstaff/d/dong.x/analysis/smchd1/long/sequins/talon/compare/talon /wehisan/home/allstaff/d/dong.x/analysis/smchd1/long/sequins/talon/sequins_talon_talon.gtf
