# The long and the short of it: unlocking nanopore long-read RNA sequencing data with short-read differential expression analysis tools

This repository contains the code used to perform the analysis and generate the figures in the paper with the same title (currently available on [bioRxiv](https://www.biorxiv.org/content/10.1101/2020.06.28.176727v1) ). Some intermidiate results of time- and resource-consuming tools are also provided.

The specific version of *FLAMES* used in this project can be accessed [here](https://github.com/XueyiDong/FLAMES).

* [`NSC`](NSC) contains the code used to analyze our NSC data.
  * [`Chen_et_al_2015_PNAS`](NSC/Chen_et_al_2015_PNAS) contains code to analyze the NSC short-read data. Other files under `NSC` are for the ONT long-read data.
  * [`preprocess`](NSC/preprocess) contains the code to map the reads to the reference genome and do gene-level counting.
  * [`QC`](NSC/QC) contains the code to calculate QC metrices and make QC plots for the NSC long-read data.
  * [`flames`](NSC/flames) contains the script and config file to run [*FLAMES*](https://github.com/LuyiTian/FLAMES) and *SQANTI*, and the code to plot the number of detected isoforms.
* [`sequins`](sequins) contains the code used to analyze our sequin data.
  * [`illumina`](sequins/illumina) contains the code used to analyze the sequins short-read data. Other files under `sequins` are for the ONT long-read data.
  * [`annotations`](sequins/annotations) contains the annotations for sequins, including decoy chromosome sequence, transcript sequences, GTF annotation file and gene/transcript length and abundance.

Our RNA-seq data can be accessed from Gene Expression Omnibus (GEO) under accession numbers GSE151984 (sequin ONT long-read data), GSE151841 (NSC ONT long-read data), GSE164598 (sequin Illumina short-read data) and GSE65747 (NSC Illumina short-read data from a previous study).
