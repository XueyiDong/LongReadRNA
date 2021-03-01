#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=100G
#SBATCH --mail-user=dong.x@wehi.edu.au
#SBATCH --mail-type=BEGIN,END,FAIL

cd /wehisan/home/allstaff/d/dong.x/analysis/2020/smchd1/NSC/flames

module load samtools
module load anaconda3

source activate
conda activate FLAMES

FLAMES=/wehisan/home/allstaff/d/dong.x/Programs/FLAMES
anno=/wehisan/home/allstaff/d/dong.x/annotation/Mouse
input=/wehisan/general/academic/seq_data/quentin/Nanopore/Smchd1-NSC-cDNA/data/rebasecall/pass/merged/used
config=/wehisan/home/allstaff/d/dong.x/analysis/2020/smchd1/NSC/flames/config_bulk_nanopore.json
dir_minimap2=/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/script/minimap2-2.17_x64-linux

python $FLAMES/python/bulk_long_pipeline.py \
 -a $anno/gencode.vM23.annotation.gff3\
 -i $input\
 -o ./results \
 --genomefa $anno/mm10/mm10.fa\
 --config_file $config\
 --minimap2 $dir_minimap2
