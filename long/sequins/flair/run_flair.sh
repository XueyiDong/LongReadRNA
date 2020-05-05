#!/bin/bash
#PBS -N flair_collapse
#PBS -q submit
#PBS -l nodes=1:ppn=8,mem=80gb,walltime=120:00:00
#PBS -M dong.x@wehi.edu.au
#PBS -m abe
#PBS -j oe


export PATH=$PATH:/wehisan/home/allstaff/d/dong.x/Programs/minimap2

module load samtools
module load bedtools
module load anaconda3
source activate
conda activate flair

cd /wehisan/home/allstaff/d/dong.x/analysis/smchd1/long/sequins/flair
anno=/wehisan/home/allstaff/d/dong.x/annotation/sequins
reads=/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/Mike_seqin/20200228_YPRDP_2xsequin_mixAB/FLTSA_output/merged_fq.fastq.gz
flairdir=/wehisan/home/allstaff/d/dong.x/Programs/flair

# align
# python $flairdir/flair.py align -g $anno/rnasequin_decoychr_2.4.fa -r $reads -t 4

# # correct 
# python $flairdir/flair.py correct -q flair.aligned.bed -f $anno/rnasequin_annotation_2.4.gtf -g $anno/rnasequin_decoychr_2.4.fa -t 4

# # collapse
# python $flairdir/flair.py collapse -g $anno/rnasequin_decoychr_2.4.fa -r $reads -q flair_all_corrected.bed -t 8 --temp_dir ./tmp

# # quantify
python $flairdir/flair.py quantify -r reads_manifest.tsv -i flair.collapse.isoforms.fa -t 8 --temp_dir ./tmp
