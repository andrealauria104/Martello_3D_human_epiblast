#
# Pseudotime analysis: Monocle3
#
# 0. Resources ----
source("scripts/resources.R")
suppressPackageStartupMessages(library(monocle3))
suppressPackageStartupMessages(library(SeuratWrappers))
suppressPackageStartupMessages(library(Seurat))

# paths
path_sce = "results/cluster-cells/sce.rds"
path_extended_data = "results/cluster-cells/extended_data.rds"
path_results = "results/cluster-cells/pseudotime-monocle3"

if(!dir.exists(path_results)) dir.create(path_results, recursive = T)

# 0.1 Global params ----
metadata_features <- c("Time","ClusterNum") # used for grouping in plots, character vector

# clustering ---
cluster_feature_name = "ClusterNum" # Cluster
use_previous_clusters = T
reduction = "umap"
cluster_cells_resolution = NULL # NULL monocle3 default
cluster_cells_k = 5 # 20 monocle3 default

# visualization ---
p_cell_size = .65

# others ---
save_object = FALSE # add pseudotime data to input object and save
save_cds_object = TRUE # save cds object
s33d = 1991
# 1. Preprocess data ----
sce <- readRDS(path_sce)

extended_data <- readRDS(path_extended_data)
pal_default <- extended_data$palette_features
if(is.null(metadata_features)) metadata_features <- names(pal_default)

# Convert object to cell_data_set ---
cds <- SeuratWrappers::as.cell_data_set(
  Seurat::as.Seurat(sce)
  )
rowData(cds) <- rowData(sce)
rowData(cds)$gene_short_name <- rowData(cds)$gene_name

# 2. Pseudotime analysis ----
cds <- monocle3::cluster_cells(
  cds = cds,
  reduction_method = "UMAP",
  resolution=cluster_cells_resolution,
  k = cluster_cells_k
  )

if(use_previous_clusters) {
  # WORKAROUND to use existing clustering results ---
  to_clusters <- colData(cds)[[cluster_feature_name]]
  names(to_clusters) <- colnames(cds)
  cds@clusters[[toupper(reduction)]]$clusters <- to_clusters
  # ---
}

cds <- monocle3::learn_graph(
  cds,
  use_partition = F,
  close_loop = F,
  learn_graph_control = list(minimal_branch_len=12)
  
  )
cds <- monocle3::order_cells(cds, root_pr_nodes = "Y_1") # Y_1 #Y_44

# Retreive pseudotime ---
colData(cds)$Pseudotime <- monocle3::pseudotime(cds)
colData(cds)$Pseudotime_scaled <- (colData(cds)$Pseudotime-min(colData(cds)$Pseudotime))/max(colData(cds)$Pseudotime)-min(colData(cds)$Pseudotime)

# 3. Visualization ----
pumap_trajectory_principal_points <- monocle3::plot_cells(
  cds,
  color_cells_by = "pseudotime",
  cell_size = p_cell_size,
  label_branch_points=T,
  label_leaves = T,
  label_roots = F,
  label_cell_groups = F,
  trajectory_graph_segment_size = .5,
  label_principal_points = T,
  rasterize = T
  ) + 
  theme_bw() + 
  my_theme

outfile <- paste0(path_results,"/dimplot.",reduction,".trajectory_principal_points.pdf")
pdf(file = outfile, paper = "a4r",w=unit(4,'cm'),h=unit(6,"cm"),useDingbats = F)
print(pumap_trajectory_principal_points)
dev.off()

pumap_trajectory_pseudotime <- monocle3::plot_cells(
  cds,
  color_cells_by = "pseudotime",
  cell_size = p_cell_size,
  label_branch_points=F,
  label_leaves = F,
  label_roots = F,
  label_cell_groups = F,
  trajectory_graph_segment_size = .5,
  label_principal_points = F,
  rasterize = T
  ) + 
  theme_bw() + 
  my_theme +
  guides(col = guide_colourbar(barwidth = 0.6))


pumap_trajectory_l <- list()
for(feature in metadata_features) {
  
  pumap_trajectory_l[[feature]] <- monocle3::plot_cells(
    cds,
    color_cells_by = feature,
    cell_size = p_cell_size,
    label_branch_points=F,
    label_leaves = F,
    label_roots = F,
    label_cell_groups = F,
    trajectory_graph_segment_size = .5,
    rasterize = T,
    show_trajectory_graph = F
    ) + 
    theme_bw() + 
    my_theme +
    guides(col = guide_legend(ncol=1, override.aes = list(size=3))) +
    scale_color_manual(values = pal_default[[feature]]) 
}

pumap_trajectory_l <- c(list("Pseudotime"=pumap_trajectory_pseudotime),pumap_trajectory_l)
pumap_trajectory <- patchwork::wrap_plots(pumap_trajectory_l, nrow = 1)
pumap_trajectory <- pumap_trajectory + patchwork::plot_annotation(
  title = "Pseudotemporal ordering",
  subtitle = paste0("(n. cells = ",dim(cds)[2],")"),
  theme = theme(
    plot.title = element_text(size=8,hjust = 0.5),
    plot.subtitle = element_text(size=8,hjust = .5)
  )
)

outfile <- paste0(path_results,"/dimplot.",reduction,".trajectory.pdf")
pdf(file = outfile, paper = "a4r",w=unit(2.6*(length(metadata_features)+1),'cm'),h=unit(6,"cm"))
print(pumap_trajectory)
dev.off()

# 4. Save object ----
if(save_object) {
  sce@meta.data <- as.data.frame(colData(cds))
  outfile <- paste0(path_results,"/",basename(path_sce))
  message(" -- saving sce object to: ", outfile)
  saveRDS(sce, file = outfile)
}

if(save_cds_object) {
  outfile <- paste0(path_results,"/cds.rds")
  message(" -- saving cds object to: ", outfile)
  saveRDS(cds, file = outfile)
}
