#!/bin/bash
#PBS -N talon
#PBS -q submit
#PBS -l nodes=1:ppn=4,mem=80gb,walltime=120:00:00
#PBS -M dong.x@wehi.edu.au
#PBS -m abe
#PBS -j oe

# TALON_WF=/stornext/General/data/user_managed/grpu_mritchie_1/HasaruK/talon_workflow
SEQUINS_DIR=/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/Mike_seqin
DATA_DIR=$SEQUINS_DIR/20200228_YPRDP_2xsequin_mixAB/fastq_pass/merged
OUT_DIR=/wehisan/home/allstaff/d/dong.x/analysis/smchd1/long/sequins/talon
# # map
# $TALON_WF/minimap2/minimap2 -ax splice --MD -c $SEQUINS_DIR/annotations/rnasequin_decoychr_2.4.fa $DATA_DIR/barcode01.fastq > $OUT_DIR/barcode01.sam
# $TALON_WF/minimap2/minimap2 -ax splice --MD -c $SEQUINS_DIR/annotations/rnasequin_decoychr_2.4.fa $DATA_DIR/barcode02.fastq > $OUT_DIR/barcode02.sam
# $TALON_WF/minimap2/minimap2 -ax splice --MD -c $SEQUINS_DIR/annotations/rnasequin_decoychr_2.4.fa $DATA_DIR/barcode03.fastq > $OUT_DIR/barcode03.sam
# $TALON_WF/minimap2/minimap2 -ax splice --MD -c $SEQUINS_DIR/annotations/rnasequin_decoychr_2.4.fa $DATA_DIR/barcode04.fastq > $OUT_DIR/barcode04.sam
# # convert to bam
# samtools view -S -b $OUT_DIR/barcode01.sam > $OUT_DIR/barcode01.bam
# samtools view -S -b $OUT_DIR/barcode02.sam > $OUT_DIR/barcode02.bam
# samtools view -S -b $OUT_DIR/barcode03.sam > $OUT_DIR/barcode03.bam
# samtools view -S -b $OUT_DIR/barcode04.sam > $OUT_DIR/barcode04.bam
# # index bam
# samtools index $OUT_DIR/barcode01.bam
# samtools index $OUT_DIR/barcode02.bam
# samtools index $OUT_DIR/barcode03.bam
# samtools index $OUT_DIR/barcode04.bam
# # run TranscriptClean
# module load samtools
# module load bedtools
# module load anaconda3
# source activate
# conda activate TranscriptClean2
# module load R
# python $TALON_WF/TranscriptClean/TranscriptClean.py --sam $OUT_DIR/barcode01.sam --genome $SEQUINS_DIR/annotations/rnasequin_decoychr_2.4.fa --outprefix $OUT_DIR/barcode01
# Rscript $TALON_WF/TranscriptClean/generate_report.R $OUT_DIR/barcode01
# python $TALON_WF/TranscriptClean/TranscriptClean.py --sam $OUT_DIR/barcode02.sam --genome $SEQUINS_DIR/annotations/rnasequin_decoychr_2.4.fa --outprefix $OUT_DIR/barcode02
# Rscript $TALON_WF/TranscriptClean/generate_report.R $OUT_DIR/barcode02
# python $TALON_WF/TranscriptClean/TranscriptClean.py --sam $OUT_DIR/barcode03.sam --genome $SEQUINS_DIR/annotations/rnasequin_decoychr_2.4.fa --outprefix $OUT_DIR/barcode03
# Rscript $TALON_WF/TranscriptClean/generate_report.R $OUT_DIR/barcode03
# python $TALON_WF/TranscriptClean/TranscriptClean.py --sam $OUT_DIR/barcode04.sam --genome $SEQUINS_DIR/annotations/rnasequin_decoychr_2.4.fa --outprefix $OUT_DIR/barcode04
# Rscript $TALON_WF/TranscriptClean/generate_report.R $OUT_DIR/barcode04
# # run TALON
# module unload python
cd /wehisan/home/allstaff/d/dong.x/analysis/smchd1/long/sequins/talon

module load anaconda3
source activate
conda activate TALON


# initialize the database
talon_initialize_database --f $SEQUINS_DIR/annotations/rnasequin_annotation_2.4.gtf --g sequins2.4 --a 2.4 --o sequins_talon

# # internal priming check
# talon_label_reads --f sam/barcode01_clean.sam --g $SEQUINS_DIR/annotations/rnasequin_decoychr_2.4.fa --t 4 --deleteTmp --o labeled/barcode01_clean
# talon_label_reads --f sam/barcode02_clean.sam --g $SEQUINS_DIR/annotations/rnasequin_decoychr_2.4.fa --t 4 --deleteTmp --o labeled/barcode02_clean
# talon_label_reads --f sam/barcode03_clean.sam --g $SEQUINS_DIR/annotations/rnasequin_decoychr_2.4.fa --t 4 --deleteTmp --o labeled/barcode03_clean
# talon_label_reads --f sam/barcode04_clean.sam --g $SEQUINS_DIR/annotations/rnasequin_decoychr_2.4.fa --t 4 --deleteTmp --o labeled/barcode04_clean

# create config file
echo "barcode01,mixA,ONT,labeled/barcode01_clean_labeled.sam
barcode02,mixA,ONT,labeled/barcode02_clean_labeled.sam
barcode03,mixB,ONT,labeled/barcode03_clean_labeled.sam
barcode04,mixB,ONT,labeled/barcode04_clean_labeled.sam" > sequins_config.csv

# run talon
talon --f sequins_config.csv --db sequins_talon.db --build sequins2.4 --o sequins_talon

# summarize how many of each transcript were found
talon_summarize --db sequins_talon.db --o sequins_talon

# filter transcripts
talon_filter_transcripts \
       --db sequins_talon.db \
       -a 2.4 \
       --maxFracA 0.5 \
       --minCount 5 \
       --o filtered_transcripts.csv

# create abundance file using whitelist
talon_abundance \
       --db sequins_talon.db \
       --whitelist filtered_transcripts.csv \
       -a 2.4 \
       --build sequins2.4 \
       --o sequins_talon


# create GTF file from database
talon_create_GTF --db sequins_talon.db --annot 2.4 --whitelist filtered_transcripts.csv --build sequins2.4 --o sequins_talon


