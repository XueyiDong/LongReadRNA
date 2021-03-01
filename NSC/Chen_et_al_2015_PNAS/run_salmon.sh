#!/bin/bash
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=20G
#SBATCH --mail-user=dong.x@wehi.edu.au
#SBATCH --mail-type=BEGIN,END,FAIL

annodir=/wehisan/home/allstaff/d/dong.x/analysis/2020/smchd1/NSC/flames/results
trfa=transcript_assembly.fa
salmondir=/wehisan/home/allstaff/d/dong.x/Programs/salmon/bin/salmon
datadir=/wehisan/general/academic/seq_data/kelan/AGRF_CAGRF6583

$salmondir index -t $annodir/$trfa -i $annodir/transcript_assembly.index -p 8

mkdir -p $datadir/salmon

for SAMPLE in 113_C26VPACXX_ATCACG 114_C26VPACXX_CGATGT 11_C26VPACXX_GATCAG 12_C26VPACXX_TAGCTT 13_C26VPACXX_GGCTAC 14_C26VPACXX_CTTGTA
do mkdir -p $datadir/salmon/$SAMPLE
$salmondir -i $annodir/transcript_assembly.index -l A -1 $datadir/NPC_RNA_$SAMPLE\_L004_R1.fastq.gz -2 $datadir/NPC_RNA_$SAMPLE\_L004_R2.fastq.gz --validateMappings -o $datadir/salmon/$SAMPLE
done