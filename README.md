# The long and the short of it: unlocking nanopore long-read RNA sequencing data with short-read differential expression analysis tools

This repository contains the code used to perform the analysis and generate the figures in the paper with the same name (currently available on [bioRxiv](https://www.biorxiv.org/content/10.1101/2020.06.28.176727v1) ). 

Our RNA-seq data can be accessed from Gene Expression Omnibus (GEO) under accession numbers GSE151984 (sequin ONT long-read data), GSE151841 (NSC ONT long-read data), GSE164598 (sequin Illumina short-read data) and GSE65747 (NSC Illumina short-read data from a previous study).



* `NSC` :This folder contains the code used to analyze our NSC datasets.
  * `Chen_et_al_2015_PNAS` contains scripts to analyze the NSC short-read dataset. 
  * `preprocess` contains the scripts to basecall, trim and demultiplex the raw fast5 data, map the reads to the reference genome and count them.