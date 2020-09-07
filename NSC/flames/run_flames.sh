#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=16
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
# out=/wehisan/home/allstaff/d/dong.x/analysis/smchd1/long/scFLT
# out=/wehisan/general/old-prkfs2/disk503/GP_Transfer/Smchd1/long/scFLT
config=/wehisan/home/allstaff/d/dong.x/analysis/2020/smchd1/NSC/flames/config_bulk_nanopore.json

# python $FLAMES/python/bulk_long_pipeline.py \
#  -a $anno/gencode.vM23.annotation.gff3\
#  -i $input\
#  -o ./results \
#  --genomefa $anno/mm10/mm10.fa\
#  --config_file $config\
#  --minimap2 /wehisan/home/allstaff/d/dong.x/Programs/minimap2

python $FLAMES/python/sc_long_pipeline.py \
 -a $anno/gencode.vM23.annotation.gff3\
 -i ./results/merged.fastq.gz\
 -o ./results \
 --genomefa $anno/mm10/mm10.fa\
 --config_file $config\
 --minimap2 /wehisan/home/allstaff/d/dong.x/Programs/minimap2
