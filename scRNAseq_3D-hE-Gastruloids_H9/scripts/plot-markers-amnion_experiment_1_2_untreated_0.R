# Visualize results
setwd("/home/epigen/Martello_3D_human_epiblast/Martello_Polo_TGFb/analysis_oscab")
# 0. Resources ----
source("scripts/resources.R")

# paths
path_sce = "results/correct-batches_experiment_1_2_untreated_0/cluster-cells/sce.rds"
path_extended_data = "results/correct-batches_experiment_1_2_untreated_0/cluster-cells/extended_data.rds"
path_amnion_markers = "data/Amnion_markers_human.xlsx"
path_results = "results/correct-batches_experiment_1_2_untreated_0/visualize-results"

if(!dir.exists(path_results)) dir.create(path_results, recursive = T)

s33d = 1991
# 1. Plot amnion markers ----
sce <- readRDS(path_sce)
extended_data <- readRDS(path_extended_data)
palette_features <- extended_data$palette_features

amnion_markers <- openxlsx::read.xlsx(path_amnion_markers, sheet = 1)

# plot
expr_mat <- as.data.frame(logcounts(sce))
expr_mat <- expr_mat[rownames(expr_mat) %in% amnion_markers$gene_name,]
expr_mat <- cbind("gene_name"=rownames(expr_mat),expr_mat)
expr_mat <- reshape2::melt(
  expr_mat,
  id.var="gene_name",
  variable.name = "Sample",
  value.name = "expression"
  )
expr_mat <- as.data.frame(merge(expr_mat,colData(sce),by="Sample"))
expr_mat <- merge(
  expr_mat,
  amnion_markers,
  by = "gene_name"
)

# by cluster
markers_by_clust_plt <- ggplot(
  expr_mat,
  aes(
    x=ClusterNum,
    y=expression,
    col=ClusterNum
    )
) + facet_wrap(type~gene_name,scales = "free_x",ncol=6) +
  # geom_violin(scale = "width") +
  geom_boxplot(outlier.shape = NA,show.legend = F) +
  ggrastr::rasterise(geom_jitter(size=0.4,width = 0.3, show.legend = F),dpi=300) +
  scale_color_manual(values = palette_features$ClusterNum) +
  scale_y_continuous(labels = scales::number_format(accuracy = 1)
                     # ,limits = c(0,10)
                     ) +
  theme_bw() + my_theme 

outfile <- paste0(path_results,"/amnion_markers_by_clust_plt.pdf")
message(" -- saving to: ", outfile)
pdf(file = outfile, paper = "a4r", width = unit(8,'cm'), height = unit(5.5,'cm'))
print(markers_by_clust_plt)
dev.off()

# by time
markers_by_time_plt <- ggplot(
  expr_mat,
  aes(
    x=Time,
    y=expression,
    col=Time
  )
) + facet_wrap(type~gene_name,scales = "free_x",ncol=6) +
  # geom_violin(scale = "width") +
  geom_boxplot(outlier.shape = NA,show.legend = F) +
  ggrastr::rasterise(geom_jitter(size=0.4,width = 0.3, show.legend = F),dpi=300) +
  scale_color_manual(values = palette_features$Time) +
  scale_y_continuous(labels = scales::number_format(accuracy = 1)
                     # ,limits = c(0,10)
  ) +
  theme_bw() + my_theme + theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))

outfile <- paste0(path_results,"/amnion_markers_by_time_plt.pdf")
message(" -- saving to: ", outfile)
pdf(file = outfile, paper = "a4r", width = unit(8.5,'cm'), height = unit(5.8,'cm'))
print(markers_by_time_plt)
dev.off()
