tracking <- read.table("comp.tracking", stringsAsFactors = FALSE)
library(ggplot2)
pdf("compare_category_sequins.pdf", height = 4)
ggplot(tracking, aes(x = V4, fill=V4)) +
  geom_bar() +
  geom_text(stat="count", aes(label=..count..), vjust=-0.2) +
  theme_bw()
dev.off()

# barplot by  read number but not number of transcripts

count <-read.csv("/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/Mike_seqin/20200228_YPRDP_2xsequin_mixAB/FLTSA_output/transcript_count.csv")

info <- strsplit2(tracking$V5, "|", fixed=TRUE)
tracking$transcript <- sub("transcript:", "", info[,2])
count$all <- rowSums(count[,3:ncol(count)])
count$category <- tracking$V4[match(count$transcript_id, tracking$transcript)]
head(count)
# known <- grepl("^ENSMUST", count$transcript_id)
# count$category[known] <- "="
count[(is.na(count$category)),]
count[(is.na(count$category)), "category"] = "="
table(count$category)
count_cat <- aggregate(count$all, by = list(count$category), FUN=sum)
colnames(count_cat) <- c("category", "count")
pdf("compare_category_count_sequins.pdf", height = 4)
ggplot(count_cat, aes(x = category, y=count, fill=category)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label=count), vjust = -0.1, size = 3) +
  theme_bw()
dev.off()

# calc length bias.
# tx length: know tx from anno, novel from GTF...but in this case let's calc from GTF
library(GenomicFeatures)
anno_FLTSA <- makeTxDbFromGFF("/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/Mike_seqin/20200228_YPRDP_2xsequin_mixAB/FLTSA_output/isoform_annotated.gff3")
tl <- transcriptLengths(anno_FLTSA)
count$tx_len <- tl$tx_len[match(count$transcript_id, tl$tx_name)]
anno <- read.delim("/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/Mike_seqin/annotations/rnasequin_isoforms_2.4.tsv", sep = "\t", stringsAsFactors = FALSE)
count$tx_len[is.na(count$tx_len)] <- anno$LENGTH[match(count$transcript_id[is.na(count$tx_len)], anno$NAME)]

library(tidyr)
count_tidy <- gather(count, "sample", "count", 3:6)
pdf("lengthVsCountSequins.pdf")
ggplot(count, aes(x=tx_len, y=all)) +
  geom_point(aes(colour=category)) +
  scale_x_continuous(trans="log2") +
  scale_y_continuous(trans="log2") +
  theme_bw() +
  labs(x="transcript length", y="") 
dev.off()
cor(count$tx_len, count$all)
# -0.07727279

sessionInfo()

# R version 3.6.1 (2019-07-05)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: CentOS Linux 7 (Core)
# 
# Matrix products: default
# BLAS:   /stornext/System/data/apps/R/R-3.6.1/lib64/R/lib/libRblas.so
# LAPACK: /stornext/System/data/apps/R/R-3.6.1/lib64/R/lib/libRlapack.so
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8   
# [6] LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] grid      parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] tidyr_1.0.0                              data.table_1.12.8                        AnnotationHub_2.16.1                    
# [4] BiocFileCache_1.8.0                      dbplyr_1.4.2                             Gviz_1.28.3                             
# [7] stageR_1.6.0                             SummarizedExperiment_1.14.1              DelayedArray_0.10.0                     
# [10] BiocParallel_1.20.1                      matrixStats_0.55.0                       ggplot2_3.2.1                           
# [13] RColorBrewer_1.1-2                       Mus.musculus_1.3.1                       TxDb.Mmusculus.UCSC.mm10.knownGene_3.4.7
# [16] org.Mm.eg.db_3.8.2                       GO.db_3.8.2                              OrganismDbi_1.26.0                      
# [19] GenomicFeatures_1.36.4                   GenomicRanges_1.36.1                     GenomeInfoDb_1.20.0                     
# [22] AnnotationDbi_1.46.1                     IRanges_2.18.3                           S4Vectors_0.22.1                        
# [25] Biobase_2.44.0                           BiocGenerics_0.30.0                      DRIMSeq_1.12.0                          
# [28] edgeR_3.26.8                             limma_3.40.6                            
# 
# loaded via a namespace (and not attached):
#   [1] colorspace_1.4-1              ellipsis_0.3.0                biovizBase_1.32.0             htmlTable_1.13.3             
# [5] XVector_0.24.0                base64enc_0.1-3               dichromat_2.0-0               rstudioapi_0.10              
# [9] farver_2.0.3                  bit64_0.9-7                   interactiveDisplayBase_1.22.0 fansi_0.4.1                  
# [13] splines_3.6.1                 knitr_1.26                    zeallot_0.1.0                 Formula_1.2-3                
# [17] Rsamtools_2.0.3               cluster_2.1.0                 png_0.1-7                     graph_1.62.0                 
# [21] shiny_1.4.0                   BiocManager_1.30.10           compiler_3.6.1                httr_1.4.1                   
# [25] backports_1.1.5               fastmap_1.0.1                 assertthat_0.2.1              Matrix_1.2-18                
# [29] lazyeval_0.2.2                cli_2.0.1                     later_1.0.0                   acepack_1.4.1                
# [33] htmltools_0.4.0               prettyunits_1.1.0             tools_3.6.1                   gtable_0.3.0                 
# [37] glue_1.3.1                    GenomeInfoDbData_1.2.1        reshape2_1.4.3                dplyr_0.8.3                  
# [41] rappdirs_0.3.1                Rcpp_1.0.3                    vctrs_0.2.1                   Biostrings_2.52.0            
# [45] rtracklayer_1.44.4            xfun_0.11                     stringr_1.4.0                 mime_0.8                     
# [49] lifecycle_0.1.0               ensembldb_2.8.1               XML_3.99-0.3                  zlibbioc_1.30.0              
# [53] scales_1.1.0                  BSgenome_1.52.0               VariantAnnotation_1.30.1      promises_1.1.0               
# [57] hms_0.5.3                     ProtGenerics_1.16.0           RBGL_1.60.0                   AnnotationFilter_1.8.0       
# [61] yaml_2.2.0                    curl_4.3                      memoise_1.1.0                 gridExtra_2.3                
# [65] biomaRt_2.40.5                rpart_4.1-15                  latticeExtra_0.6-29           stringi_1.4.5                
# [69] RSQLite_2.2.0                 checkmate_1.9.4               rlang_0.4.4                   pkgconfig_2.0.3              
# [73] bitops_1.0-6                  lattice_0.20-38               purrr_0.3.3                   labeling_0.3                 
# [77] GenomicAlignments_1.20.1      htmlwidgets_1.5.1             bit_1.1-15.1                  tidyselect_0.2.5             
# [81] plyr_1.8.5                    magrittr_1.5                  R6_2.4.1                      Hmisc_4.3-1                  
# [85] DBI_1.1.0                     pillar_1.4.3                  foreign_0.8-75                withr_2.1.2                  
# [89] survival_3.1-8                RCurl_1.98-1.1                nnet_7.3-12                   tibble_2.1.3                 
# [93] crayon_1.3.4                  jpeg_0.1-8.1                  progress_1.2.2                locfit_1.5-9.1               
# [97] blob_1.2.1                    digest_0.6.23                 xtable_1.8-4                  httpuv_1.5.2                 
# [101] munsell_0.5.0   
  