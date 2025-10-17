setwd("/home/epigen/Martello_3D_human_epiblast/Martello_Polo_TGFb/analysis_oscab")
source("scripts/resources.R")

path_cluster_markers_l = c(
  "results/correct-batches_experiment_1_2_untreated_0/cluster-cells/markers-one-vs-rest/de.genes.by.group.pvalue_0.01.xlsx",
  "results/correct-batches_experiment_1_2_untreated_0/cluster-cells/markers-one-vs-rest/de.genes.by.group.pvalue_adj_0.05.xlsx"
  )
path_sce = "results/correct-batches_experiment_1_2_untreated_0/cluster-cells/sce.rds"

sce <- readRDS(path_sce)

for(path_cluster_markers in path_cluster_markers_l) {
  cluster_nm <- openxlsx::getSheetNames(path_cluster_markers)
  cluster_markers <- lapply(
    cluster_nm, 
    function(i) openxlsx::read.xlsx(path_cluster_markers, sheet = i)
  )
  names(cluster_markers) <- colData(sce)$ClusterLabels[match(cluster_nm,colData(sce)$Cluster)]
  cluster_markers <- cluster_markers[order(names(cluster_markers))]
  
  to_save <- cluster_markers
  names(to_save) <- gsub("\\:","\\-",names(to_save))
  saveXLSresEdgeR2(
    res = to_save,
    outfile = gsub(".pvalue_",".cluster_labels.pvalue_",path_cluster_markers)
  )
  rm(to_save)
  
}
