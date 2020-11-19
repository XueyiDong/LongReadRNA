#!/bin/bash
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=20G
#SBATCH --mail-user=dong.x@wehi.edu.au
#SBATCH --mail-type=BEGIN,END,FAIL

module load subread
module load samtools

folder=/stornext/General/data/user_managed/grpu_mritchie_1/long_RNA_benchmark/Illumina_April19
idx=/wehisan/home/allstaff/d/dong.x/annotation/sequins/rnasequin.index

for sample in 1-H1975-1_S1 2-H1975-2_S2 6-HCC-1_S6 7-HCC-2_S7
do subread-align -T 8 -t 0 -i $idx -r $folder/data/trimmed/$sample\_R1_val_1.fq.gz -R $folder/data/trimmed/$sample\_R2_val_2.fq.gz |samtools sort > $folder/results/subread_sequins/$sample.sorted.bam
done