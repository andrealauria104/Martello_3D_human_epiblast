#
# Pseudotime analysis: TSCAN
#
# 0. Resources ----
source("scripts/resources.R")
suppressPackageStartupMessages(library(TSCAN))

# paths
path_sce = "results/correct-batches_experiment_1_2_untreated/cluster-cells/sce.rds"
path_extended_data = "results/correct-batches_experiment_1_2_untreated/cluster-cells/extended_data.rds"
path_results = "results/correct-batches_experiment_1_2_untreated/cluster-cells/pseudotime-tscan"

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
s33d = 1991
# 1. Preprocess data ----
sce <- readRDS(path_sce)

extended_data <- readRDS(path_extended_data)
palette_features <- extended_data$palette_features

# 2. Pseudotime analysis ----
by.cluster <- aggregateAcrossCells(sce, ids=sce[[cluster_feature_name]], statistics="mean")
centroids <- reducedDim(by.cluster, "corrected")

# Retreive pseudotime ---
mst <- TSCAN::createClusterMST(centroids, clusters=NULL)
set.seed(s33d)
mst <- TSCAN::createClusterMST(
  reducedDim(sce, "corrected")[,1:50],
  clusters = sce$ClusterNum,
  with.mnn = T,
  mnn.k = 10
  )
plot(mst)
line.data <- TSCAN::reportEdges(
  by.cluster,
  mst=mst,
  clusters=NULL,
  use.dimred="UMAP"
  )

map.tscan <- TSCAN::mapCellsToEdges(
  sce,
  mst=mst,
  use.dimred="corrected_sel",
  clusters=NULL
  )
tscan.pseudo <- TSCAN::orderCells(map.tscan, mst,start = 0)

for(i in 1:ncol(tscan.pseudo)) {
  path_nm <- paste0("path_",i)
  sce[[path_nm]] <- FALSE
  sce[[path_nm]][-which(is.na(tscan.pseudo[,i]))] <- TRUE
}

# 3. Visualization ----
umap_plt_l <- list()
for(feature in metadata_features) {
  
  umap_plt_l[[feature]] <- plotsceReducedDim(
    sce = sce,
    dimred = "UMAP",
    col_by = feature,
    point_size = 0.3,
    pal = palette_features[[feature]]
  )
  
  if(feature=="Cluster") {
    umap_plt_l[[feature]] <- Seurat::LabelClusters(
      plot = umap_plt_l[[feature]],
      id = feature,
      repel = F,
      col="black"
    )
  }
  umap_plt_l[[feature]] <- umap_plt_l[[feature]] +
    geom_line(line.data,mapping = aes(x=dim1, y=dim2, group=edge),inherit.aes = F)
}

umap_plt <- patchwork::wrap_plots(umap_plt_l, nrow = 1)
umap_plt <- umap_plt + patchwork::plot_annotation(
  title = "Pseudotemporal ordering (TSCAN)",
  subtitle = paste0("(n. cells = ",dim(sce)[2],")"),
  theme = theme(
    plot.title = element_text(size=8,hjust = 0.5),
    plot.subtitle = element_text(size=8,hjust = .5)
  )
)

outfile <- paste0(path_results,"/pumap_",paste0(metadata_features,collapse = "_"),".pdf")
pdf(file = outfile, paper = "a4r",w=unit(2.8*length(metadata_features),'cm'),h=unit(4.5,"cm"))
print(umap_plt)
dev.off()

# 4. Save object ----
# if(save_cds_object) {
#   outfile <- paste0(path_results,"/cds.rds")
#   message(" -- saving cds object to: ", outfile)
#   saveRDS(cds, file = outfile)
# }
