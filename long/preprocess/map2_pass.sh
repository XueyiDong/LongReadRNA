# mapping rep1 pass reads only using minimap2
cd /wehisan/general/academic/seq_data/quentin/Nanopore/Smchd1-NSC-cDNA_rep2
mkdir -p results results/minimap2_mm10_pass
cd results/minimap2_mm10_pass
module load samtools/1.7

for sample in BC{13..18}
do ~/Programs/minimap2/minimap2 -ax splice -uf -k14 --junc-bed  ~/annotation/Mouse/gencode.junction.bed ~/annotation/Mouse/mm10/mm10.fa ../../data/demultiplex_new/filtered/$sample.fastq.gz | samtools view -b | samtools sort > $sample.sorted.bam
samtools index $sample.sorted.bam $sample.sorted.bai
done
