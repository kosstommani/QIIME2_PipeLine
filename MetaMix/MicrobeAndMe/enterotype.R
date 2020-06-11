# Taeyune Kim - 04/20/2020
# From MicrobeAndMe v2.1
# Predict the microbial composition data to the database enterotype model.
# Alternatively, using three genera that Bacteriodes, Prevotella and Ruminococcus to get more clearly seperated cluster than using whole genera.
# For use in R.
#
# Usage: enterotype.R input_biom output_csv
#
# Taxonomy assignment must be performed using the 16S_ribosomal_RNA_20200305 database (META Analysis Team, Macrogen Inc.)
#

args <- commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file)\n", call.=FALSE)
}

biomIn <- args[1]
output <- args[2]

data <- read.table(biomIn, header=TRUE, row.names=1)
MAMv2.1_model <- readRDS(file='/garnet/Tools/Amplicon_MetaGenome/PipeLine_Dev/MetaMix/MicrobeAndMe/MAMv2.1_Enterotype_model.rds')
library(randomForest)
write.table(predict(MAMv2.1_model, t(data)), file = output, sep = ',', col.names = F, quote = F)