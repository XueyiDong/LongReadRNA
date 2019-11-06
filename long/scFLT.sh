#!/bin/bash
#PBS -N scFLT
#PBS -q submit
#PBS -l nodes=1:ppn=16,mem=100gb,walltime=1000:00:00
#PBS -M dong.x@wehi.edu.au
#PBS -m abe
#PBS -j oe
#PBS -o scFLT.o

cd /wehisan/home/allstaff/d/dong.x/analysis/smchd1/long/scFLT

module load samtools
module load python

sftwr=/stornext/General/data/user_managed/grpu_mritchie_1/LuyiTian/git/scFLT
anno=/wehisan/home/allstaff/d/dong.x/annotation/Mouse
input=/wehisan/general/academic/seq_data/quentin/Nanopore/Smchd1-NSC-cDNA/Xueyi/Smchd1_NSC_cDNA_scBarcode.fq.gz
out=/wehisan/home/allstaff/d/dong.x/analysis/smchd1/long/scFLT

python $sftwr/sc_long_pipeline.py \
 --gff3 $anno/gencode.vM22.annotation.gff3\
 --infq $input\
 --outdir $out \
 --genomefa $anno/mm10/mm10.fa\
 --config_file $out/config_bulk_nanopore.json\
 --minimap2 /wehisan/home/allstaff/d/dong.x/Programs/minimap2