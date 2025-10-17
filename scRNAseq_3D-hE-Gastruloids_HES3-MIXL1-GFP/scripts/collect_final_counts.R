# save counts
setwd("/home/epigen/Martello_3D_human_epiblast/Martello_Polo_MIXL1/run240820_run241003/analysis_oscab")

source("scripts/resources.R")

sce <- readRDS("results/cluster-cells/sce.rds")
counts <- counts(sce)
counts <- cbind.data.frame("gene_name"=rownames(counts),counts)

write.table(
  counts,
  file = gzfile("results/cluster-cells/raw_counts.gz"),
  row.names = F,
  col.names = T,
  sep = "\t",
  quote = F
)