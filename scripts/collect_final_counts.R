# save counts
setwd("Martello_Polo_TGFb/analysis_oscab")

source("scripts/resources.R")

sce <- readRDS("results/correct-batches_experiment_1_2_untreated_0/cluster-cells/sce.rds")
counts <- counts(sce)
counts <- cbind.data.frame("gene_name"=rownames(counts),counts)
colnames(counts) <- gsub("_Martello_Polo","",colnames(counts))
colnames(counts) <- gsub("Day","DAY",colnames(counts))
write.table(
  counts,
  file = gzfile("results/correct-batches_experiment_1_2_untreated_0/cluster-cells/raw_counts.gz"),
  row.names = F,
  col.names = T,
  sep = "\t",
  quote = F
)