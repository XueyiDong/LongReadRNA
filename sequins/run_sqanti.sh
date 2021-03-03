#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
#SBATCH --mail-user=dong.x@wehi.edu.au
#SBATCH --mail-type=BEGIN,END,FAIL

module load anaconda3
source activate
conda activate SQANTI
module load R
module load perl
module load gmap-gsnap

cd /wehisan/home/allstaff/d/dong.x/analysis/2020/smchd1/sequins

python /wehisan/home/allstaff/d/dong.x/Programs/SQANTI/sqanti_qc.py \
  /stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/Mike_seqin/20200228_YPRDP_2xsequin_mixAB/flames/isoform_annotated.filtered.gff3 \
  /wehisan/home/allstaff/d/dong.x/annotation/sequins/rnasequin_annotation_2.4.gtf \
  /wehisan/home/allstaff/d/dong.x/annotation/sequins/rnasequin_decoychr_2.4.fa -g

python /wehisan/home/allstaff/d/dong.x/Programs/SQANTI/sqanti_qc.py \
  ./talon/sequins_talon_talon.gtf \
  /wehisan/home/allstaff/d/dong.x/annotation/sequins/rnasequin_annotation_2.4.gtf \
  /wehisan/home/allstaff/d/dong.x/annotation/sequins/rnasequin_decoychr_2.4.fa -g

cd ./flair

python /wehisan/home/allstaff/d/dong.x/Programs/SQANTI/sqanti_qc.py \
  flair.collapse.isoforms.gtf \
  /wehisan/home/allstaff/d/dong.x/annotation/sequins/rnasequin_annotation_2.4.gtf \
  /wehisan/home/allstaff/d/dong.x/annotation/sequins/rnasequin_decoychr_2.4.fa -g

