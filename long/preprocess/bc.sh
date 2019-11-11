#!/bin/bash
#PBS -N basecalling
#PBS -q submit_2xp100
#PBS -l nodes=1:ppn=48,mem=16gb,walltime=72:00:00
#PBS -M dong.x@wehi.edu.au
#PBS -m abe
#PBS -j oe
#PBS -o basecalling.o
cd /wehisan/general/academic/seq_data/quentin/Nanopore/Smchd1-NSC-cDNA
module load guppy-gpu
guppy_basecaller -i data/multi_reads -s Xueyi/basecalled -r -c dna_r9.4.1_450bps_hac_prom.cfg  -x "auto"
