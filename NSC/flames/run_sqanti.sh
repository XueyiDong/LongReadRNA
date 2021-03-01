#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=20G
#SBATCH --mail-user=dong.x@wehi.edu.au
#SBATCH --mail-type=BEGIN,END,FAIL

module load anaconda3
source activate
conda activate SQANTI
module load R
module load perl
module load gmap-gsnap

cd /wehisan/home/allstaff/d/dong.x/analysis/2020/smchd1/NSC/flames

python /wehisan/home/allstaff/d/dong.x/Programs/SQANTI/sqanti_qc.py \
  results/isoform_annotated.filtered.gff3 \
  /wehisan/home/allstaff/d/dong.x/annotation/Mouse/gencode.vM24.annotation.gtf \
  /wehisan/home/allstaff/d/dong.x/annotation/Mouse/mm10/mm10.fa -g
