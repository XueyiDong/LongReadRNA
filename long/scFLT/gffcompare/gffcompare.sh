#!/bin/bash
#PBS -N gffcompare
#PBS -q submit
#PBS -l nodes=1:ppn=4,mem=60gb,walltime=1000:00:00
#PBS -M dong.x@wehi.edu.au
#PBS -m abe
#PBS -j oe

cd /wehisan/home/allstaff/d/dong.x/analysis/smchd1/long/scFLT

module load anaconda3
#conda activate flair

UNIXHOME=/wehisan/home/allstaff/d/dong.x

$UNIXHOME/Programs/gffcompare -r $UNIXHOME/annotation/Mouse/gencode.vM23.annotation.gff3 -R -o $UNIXHOME/analysis/smchd1/long/scFLT/gffcompare/scFLT $UNIXHOME/analysis/smchd1/long/scFLT/run1/isoform_annotated.gff3
