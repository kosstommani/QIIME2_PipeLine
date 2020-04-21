# Taeyune Kim - 04/20/2020
# From MicrobeAndMe v2.1
# Predict the microbial composition data to the database enterotype model.
# Alternatively, using three genera that Bacteriodes, Prevotella and Ruminococcus to get more clearly seperated cluster than using whole genera.
# For use in R.
#
# Usage: enterotype.R input_biom (output_csv)
#
# Taxonomy assignment must be performed using the 16S_ribosomal_RNA_20200305 database (META Analysis Team, Macrogen Inc.)
#

args <- commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file)\n", call.=FALSE)
} else if (length(args)==1) {
  args[2] <- "enterotype.csv" #default output file is enterotype.csv
}

biomIn <- args[1]
biomIn <- strsplit(biomIn,split=".biom")[[1]]

#BIOM without taxa
library(biomformat)
Genera_biom <- read_biom(paste0('./MAM_L6/', biomIn, '_L6.biom'))
Genera <- as.data.frame(as.matrix(biom_data(Genera_biom)))
Genera.rAB <- Genera

Genera.rAB_BPR <- Genera.rAB[
  c("Bacteria;__Bacteroidetes;__Bacteroidia;__Bacteroidales;__Bacteroidaceae;__Bacteroides",
    "Bacteria;__Bacteroidetes;__Bacteroidia;__Bacteroidales;__Prevotellaceae;__Prevotella",
    "Bacteria;__Firmicutes;__Clostridia;__Clostridiales;__Ruminococcaceae;__Ruminococcus",
  ),]

rownames(Genera.rAB_BPR) <- c("Bacteroides", "Prevotella", "Ruminococcus")
MAMv2.1_model <- readRDS(file='/garnet/Analysis/BI/AmpliconMetaGenome_MAM/test/MAMv2.1_ET_Automation/MAMv2.1_Enterotype_model.rds')
library(randomForest)
write.table(predict(MAMv2.1_model, t(Genera.rAB_BPR)), file = args[2], sep = ',', col.names = F, quote = F)
