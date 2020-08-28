#!/bin/bash
#SBATCH --partition=submit
#SBATCH --ntasks=8
#SBATCH --mem-per-cpu=100GB
#SBATCH --time=120:00:00
#SBATCH --mail-user=dong.x@wehi.edu.au
#SBATCH --mail-type=BEGIN,END,FAIL

cd /stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/Mike_seqin/20200228_YPRDP_2xsequin_mixAB

module load samtools
module load anaconda3

source activate
conda activate FLAMES

FLAMES=/wehisan/home/allstaff/d/dong.x/Programs/FLAMES
anno=/wehisan/home/allstaff/d/dong.x/annotation/sequins
input=/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/Mike_seqin/20200228_YPRDP_2xsequin_mixAB/sequins_rebasecall/pass/merged
out=/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/Mike_seqin/20200228_YPRDP_2xsequin_mixAB/flames
config=/wehisan/home/allstaff/d/dong.x/analysis/smchd1/long/scFLT/config_bulk_nanopore.json

python $sftwr/sc_long_pipeline.py \
 --gff3 $anno/rnasequin_annotation_2.4.gtf\
 --infq $input\
 --outdir $out\
 --genomefa $anno/rnasequin_decoychr_2.4.fa\
 --minimap2 /wehisan/home/allstaff/d/dong.x/Programs/minimap2
