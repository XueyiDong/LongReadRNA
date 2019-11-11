cd /wehisan/general/academic/seq_data/quentin/Nanopore/Smchd1-NSC-cDNA/Xueyi/demultiplex

sed -e '1~2 s/@/@ACGT_ABCD\#/g' BC07.fastq > BC07_scBarcode.fq 
sed -e '1~2 s/@/@AGTC_ABCD\#/g' BC08.fastq > BC08_scBarcode.fq 
sed -e '1~2 s/@/@ATCG_ABCD\#/g' BC09.fastq > BC09_scBarcode.fq 
sed -e '1~2 s/@/@CATG_ABCD\#/g' BC10.fastq > BC10_scBarcode.fq 
sed -e '1~2 s/@/@CGAT_ABCD\#/g' BC11.fastq > BC11_scBarcode.fq 
sed -e '1~2 s/@/@CTGA_ABCD\#/g' BC12.fastq > BC12_scBarcode.fq 
sed -e '1~2 s/@/@AAAA_ABCD\#/g' none.fastq > none_scBarcode.fq 

cd /wehisan/general/academic/seq_data/quentin/Nanopore/Smchd1-NSC-cDNA_rep2/data/demultiplex_new

sed -e '1~2 s/@/@GCTA_ABCD\#/g' BC13.fastq > BC13_scBarcode.fq 
sed -e '1~2 s/@/@GACT_ABCD\#/g' BC14.fastq > BC14_scBarcode.fq 
sed -e '1~2 s/@/@GTAC_ABCD\#/g' BC15.fastq > BC15_scBarcode.fq 
sed -e '1~2 s/@/@TCAG_ABCD\#/g' BC16.fastq > BC16_scBarcode.fq 
sed -e '1~2 s/@/@TGCA_ABCD\#/g' BC17.fastq > BC17_scBarcode.fq 
sed -e '1~2 s/@/@TAGC_ABCD\#/g' BC18.fastq > BC18_scBarcode.fq 
sed -e '1~2 s/@/@TTTT_ABCD\#/g' none.fastq > none_scBarcode.fq 

cd /wehisan/general/academic/seq_data/quentin/Nanopore/Smchd1-NSC-cDNA/Xueyi

dir1=/wehisan/general/academic/seq_data/quentin/Nanopore/Smchd1-NSC-cDNA/Xueyi/demultiplex
dir2=/wehisan/general/academic/seq_data/quentin/Nanopore/Smchd1-NSC-cDNA_rep2/data/demultiplex_new

cat $dir1/*_scBarcode.fq > scBarcode.fq
cat $dir2/*_scBarcode.fq >> scBarcode.fq
