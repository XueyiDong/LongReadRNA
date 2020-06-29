# mapping rep1 using minimap2
cd /wehisan/general/academic/seq_data/quentin/Nanopore/Smchd1-NSC-cDNA/Xueyi
mkdir -p results results/minimap2_mm10
cd results/minimap2_mm10
module load samtools/1.7

for sample in BC{07..12} none
do ~/Programs/minimap2/minimap2 -ax splice -uf -k14 --junc-bed  ~/annotation/Mouse/gencode.junction.bed ~/annotation/Mouse/mm10/mm10.fa ../../demultiplex/$sample.fastq | samtools view -b | samtools sort > $sample.sorted.bam
samtools index $sample.sorted.bam $sample.sorted.bai
done
