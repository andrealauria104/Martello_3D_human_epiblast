# sce <- readRDS("results/correct-batches_experiment_1_2_untreated_human_embryo/cluster-cells/sce.rds")
sce <- readRDS("results/correct-batches_experiment_1_2_untreated/cluster-cells/sce.rds")
epi_markers <- openxlsx::read.xlsx(
  "data/xiang_et_al_human_embryo/41586_2019_1875_MOESM7_ESM.Genes_EPI_Subtypes.xlsx",
  sheet = "Formatted"
)

epi_markers <- lapply(
  split(epi_markers,epi_markers$subtype),
  "[[",
  "gene_name"
)
epi_expr_mat <- assay(
  sce,
  "logcounts"
)
epi_rankings <- AUCell::AUCell_buildRankings(
  epi_expr_mat
)

tmp<-msigdbr::msigdbr(subcat="CP:KEGG")
tmp<-tmp[tmp$gs_name=="KEGG_APOPTOSIS",]
markers<-lapply(split(tmp,tmp$gs_name),function(x)unique(as.data.frame(x)[,"gene_symbol"]))
# epi_markers <- markers
epi_markers_auc <- AUCell::AUCell_calcAUC(
  geneSets = epi_markers,
  rankings = epi_rankings,aucMaxRank = 8000
)

to_plot <- reshape2::melt(epi_markers_auc@assays@data$AUC)
to_plot <- merge(
  to_plot,
  as.data.frame(colData(sce)),
  by.x="cells",
  by.y="Sample"
  )
colnames(to_plot)[2] <- "gene_set"
ggplot(
  to_plot,
  aes(
    x=Cluster,
    y=value
  )
) +
  geom_boxplot() +
  facet_grid(~gene_set) +
  ggpubr::stat_compare_means(
    comparisons = list(
      c("DAY0","DAY2"),
      c("DAY0","DAY4"),
      c("DAY0","DAY6")
      )
  )
