#!/bin/bash
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=20G
#SBATCH --mail-user=dong.x@wehi.edu.au
#SBATCH --mail-type=BEGIN,END,FAIL

annodir=/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/Mike_seqin/20200228_YPRDP_2xsequin_mixAB/flames
trfa=transcript_assembly.fa

/wehisan/home/allstaff/d/dong.x/Programs/salmon/bin/salmon index -t $annodir/$trfa -i $annodir/transcript_assembly.index -p 8

datadir=/stornext/General/data/user_managed/grpu_mritchie_1/long_RNA_benchmark/Illumina_April19
mkdir -p $datadir/results/salmon

for sample in 1-H1975-1_S1 2-H1975-2_S2 6-HCC-1_S6 7-HCC-2_S7
do mkdir -p $datadir/results/salmon/$sample
/wehisan/home/allstaff/d/dong.x/Programs/salmon/bin/salmon quant -i $annodir/transcript_assembly.index -l A -1 $datadir/data/trimmed/$sample\_R1_val_1.fq.gz -2 $datadir/data/trimmed/$sample\_R2_val_2.fq.gz --validateMappings -o $datadir/results/salmon/$sample
done