library(data.table)
# File <-
#   "/wehisan/general/academic/seq_data/quentin/Nanopore/Smchd1-NSC-cDNA/data/rebasecall/sequencing_summary.txt.gz"
# Table <- fread(File, header = TRUE, sep = "\t")
# Read_Id <- as.character(Table$read_id)
# Read_length <- as.numeric(Table$sequence_length_template)
# Qscore <- as.numeric(Table$mean_qscore_template)
# Pass_filtering <- Table$pass_filtering
# Barcode <- as.character(Table$barcode_arrangement)
# Table <- cbind(Read_Id, Read_length, Qscore, Pass_filtering, Barcode)
# saveRDS(Table, "summaryInfo.RDS")

File2 <-
  "/wehisan/general/academic/seq_data/quentin/Nanopore/Smchd1-NSC-cDNA_rep2/data/rebasecall/sequencing_summary.txt.gz"
Table <- fread(File2, header = TRUE, sep = "\t")
Read_Id <- as.character(Table$read_id)
Read_length <- as.numeric(Table$sequence_length_template)
Qscore <- as.numeric(Table$mean_qscore_template)
Pass_filtering <- Table$pass_filtering
Barcode <- as.character(Table$barcode_arrangement)
Table <- cbind(Read_Id, Read_length, Qscore, Pass_filtering, Barcode)
saveRDS(Table, "summaryInfo2.RDS")
  
  
  