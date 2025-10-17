setwd("/home/epigen/Martello_3D_human_epiblast/Martello_Polo_MIXL1/run241003/analysis_oscab")

path_mapping_stats = "data/mapping_stats.txt"
path_biotypes_stats = "data/biotypes_stats.txt"

mapping_stats <- read.delim(path_mapping_stats, stringsAsFactors = F)
biotypes_stats  <- read.delim(path_biotypes_stats, stringsAsFactors = F)

biotypes_stats <- biotypes_stats[,c("Sample",colnames(biotypes_stats)[!colnames(biotypes_stats)%in%colnames(mapping_stats)])]
metadata_full <- merge(mapping_stats,biotypes_stats,by="Sample")

write.table(
  metadata_full,
  "data/metadata_full.txt",
  row.names = F,
  col.names = T,
  sep = "\t",
  quote = F
)
