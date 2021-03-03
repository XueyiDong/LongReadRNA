#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --mail-user=dong.x@wehi.edu.au
#SBATCH --mail-type=BEGIN,END,FAIL


export PATH=$PATH:/wehisan/home/allstaff/d/dong.x/Programs/minimap2

module load samtools
module load bedtools
module load anaconda3
source activate
conda activate flair

cd /wehisan/home/allstaff/d/dong.x/analysis/2020/smchd1/sequins/flair
anno=/wehisan/home/allstaff/d/dong.x/annotation/sequins
reads=/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/Mike_seqin/20200228_YPRDP_2xsequin_mixAB/sequins_rebasecall/pass/merged
flairdir=/wehisan/home/allstaff/d/dong.x/Programs/flair

# align
python $flairdir/flair.py align -g $anno/rnasequin_decoychr_2.4.fa -r $reads/barcode01.fq.gz $reads/barcode02.fq.gz $reads/barcode03.fq.gz $reads/barcode04.fq.gz -t 8 -v1.3

# # correct 
python $flairdir/flair.py correct -q flair.aligned.bed -f $anno/rnasequin_annotation_2.4.gtf -g $anno/rnasequin_decoychr_2.4.fa -t 8

# # collapse
python $flairdir/flair.py collapse -g $anno/rnasequin_decoychr_2.4.fa -f $anno/rnasequin_annotation_2.4.gtf -r $reads/barcode01.fq.gz $reads/barcode02.fq.gz $reads/barcode03.fq.gz $reads/barcode04.fq.gz -q flair_all_corrected.bed -t 8 --temp_dir ./tmp

# # quantify
python $flairdir/flair.py quantify -r reads_manifest.tsv -i flair.collapse.isoforms.fa -t 8 --temp_dir ./tmp
