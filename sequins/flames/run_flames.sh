#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=100G
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
config=/wehisan/home/allstaff/d/dong.x/analysis/2020/smchd1/sequins/flames/config_sequin.json

# python $FLAMES/python/bulk_long_pipeline.py \
#  -a $anno/rnasequin_annotation_2.4.gtf\
#  -i $input\
#  -o $out\
#  -f $anno/rnasequin_decoychr_2.4.fa\
#  -c $config\
#  -m /wehisan/home/allstaff/d/dong.x/Programs/minimap2

python $FLAMES/python/sc_long_pipeline.py \
 -a $anno/rnasequin_annotation_2.4.gtf\
 -i $out/merged.fastq.gz\
 -o $out\
 -f $anno/rnasequin_decoychr_2.4.fa\
 -c $config\
 -m /wehisan/home/allstaff/d/dong.x/Programs/minimap2