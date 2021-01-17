# This file is from Soneson et al., 2019, Nat Comm
readBam <- function(bamfile) {
  bam <- readGAlignments(bamfile, use.names = TRUE,
                         param = ScanBamParam(tag = c("NM"),
                                              what = c("qname","flag", "rname", 
                                                       "pos", "mapq")))
  ops <- GenomicAlignments::CIGAR_OPS
  wdths <- GenomicAlignments::explodeCigarOpLengths(cigar(bam), ops = ops)
  keep.ops <- GenomicAlignments::explodeCigarOps(cigar(bam), ops = ops)
  explodedcigars <- IRanges::CharacterList(relist(paste0(unlist(wdths), 
                                                         unlist(keep.ops)), wdths))
  for (opts in setdiff(GenomicAlignments::CIGAR_OPS, "=")) {
    mcols(bam)[[paste0("nbr", opts)]] <- 
      sapply(explodedcigars, function(cg) sum(as.numeric(gsub(paste0(opts, "$"), "", cg)), na.rm = TRUE))
  }
  mcols(bam)$readLength <- rowSums(as.matrix(mcols(bam)[, c("nbrS", "nbrH", "nbrM", "nbrI")]))
  bam
}

makeReadDf <- function(bam) {
  tmp <- data.frame(bam %>% setNames(NULL), stringsAsFactors = FALSE) %>%
    dplyr::rename(read = qname,
                  nbrJunctions = njunc) %>%
    dplyr::select(-cigar) %>%
    dplyr::mutate(alignedLength = nbrM + nbrI) ## equivalent to readLength-nbrS-nbrH
  
  tmp2 <- as.data.frame(table(names(subset(bam, flag %in% c(0, 16)))))
  if (nrow(tmp2) == 0) tmp2 <- data.frame(Var1 = tmp$read[1], Freq = 0)
  tmp <- tmp %>% 
    dplyr::left_join(tmp2 %>% dplyr::rename(read = Var1, nbrPrimaryAlignments = Freq))
  
  tmp3 <- as.data.frame(table(names(subset(bam, flag %in% c(256, 272)))))
  if (nrow(tmp3) == 0) tmp3 <- data.frame(Var1 = tmp$read[1], Freq = 0)
  tmp <- tmp %>% 
    dplyr::left_join(tmp3 %>% dplyr::rename(read = Var1, nbrSecondaryAlignments = Freq))
  
  tmp4 <- as.data.frame(table(names(subset(bam, flag %in% c(2048, 2064)))))
  if (nrow(tmp4) == 0) tmp4 <- data.frame(Var1 = tmp$read[1], Freq = 0)
  tmp <- tmp %>% 
    dplyr::left_join(tmp4 %>% dplyr::rename(read = Var1, nbrSupplementaryAlignments = Freq))
  
  tmp %>% dplyr::mutate(nbrSecondaryAlignments = replace(nbrSecondaryAlignments, 
                                                         is.na(nbrSecondaryAlignments), 0),
                        nbrSupplementaryAlignments = replace(nbrSupplementaryAlignments, 
                                                             is.na(nbrSupplementaryAlignments), 0))
}

makeSummaryList <- function(bam) {
  list(nAlignments = length(bam), 
       nAlignedReads = length(unique(names(bam))),
       nPrimaryAlignments = length(subset(bam, flag %in% c(0, 16))), 
       nReadsWithPrimaryAlignments = length(unique(names(subset(bam, flag %in% c(0, 16))))),
       nSecondaryAlignments = length(subset(bam, flag %in% c(256, 272))),
       nReadsWithSecondaryAlignments = length(unique(names(subset(bam, flag %in% c(256, 272))))),
       nSupplementaryAlignments = length(subset(bam, flag %in% c(2048, 2064))),
       nReadsWithSupplementaryAlignments = length(unique(names(subset(bam, flag %in% c(2048, 2064)))))
  )
}