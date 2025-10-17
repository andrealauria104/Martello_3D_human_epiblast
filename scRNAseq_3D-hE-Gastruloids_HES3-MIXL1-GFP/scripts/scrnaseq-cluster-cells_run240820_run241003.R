# Analyze scRNAseq - Bioconductor workflow
# Cluster cells
# 0. Resources ----
source("scripts/resources.R")

# paths
path_sce = "results//sce.rds"
path_results = "results/cluster-cells"

if(!dir.exists(path_results)) dir.create(path_results, recursive = T)

s33d = 1991
# 0.1 Global params ----
assay.type <- "logcounts"

batch_var  <- NULL #"batch"
design_var <- "Time"
shape_by   <- NULL
blocking_method <- "design" # block | design
  
metadata_features <- c("batch","Time","Cluster")

k_snn_graph <- 10

fdrTh        <- 0.05
fcTh         <- 0
test.type    <- "wilcox" # wilcox, t
markers_pval <- c("all","any","some")
markers_ntop <- 50
gene_labels <- NULL

save_additional_files <- FALSE
markers_heatmap       <- TRUE
markers_tsne          <- FALSE
# 1. Load and prepare data -----
# Load data ---
sce <- readRDS(path_sce)

# 2. Clustering ----
dimred_method <- "PCA.elbow"

# Graph-based clustering ---
if(is.null(k_snn_graph)) {
  if(floor(sqrt(ncol(sce)))%%2==0) {
    # keep it odd
    k_snn_graph <- floor(sqrt(ncol(sce)))+1  
  } else {
    k_snn_graph <- floor(sqrt(ncol(sce)))
  }
}
print(paste0("k: ",k_snn_graph))
g <- scran::buildSNNGraph(sce, k=k_snn_graph, use.dimred = dimred_method, type="rank")
clust <- igraph::cluster_louvain(g)$membership

sce$Cluster <- factor(clust)
table(sce$Cluster)

cluster <- colData(sce)[,"Cluster"]
names(cluster) <- colnames(sce)
# Palettes for metadata ---
metadata_features <- get_comma_arglist(metadata_features)
palette_features  <- get_palette_features(metadata_features)

for(feature in metadata_features) {
  names(palette_features[[feature]]) <- unique(colData(sce)[,feature])[order(unique(colData(sce)[,feature]))]
}

# palette_features$Cluster <- RColorBrewer::brewer.pal(9,"Set3")[-c(2,5)][1:length(levels(sce$Cluster))]
# names(palette_features$Cluster) <- levels(sce$Cluster)
palette_features$Cluster <- c(
  RColorBrewer::brewer.pal(9,"Reds")[4],
  RColorBrewer::brewer.pal(9,"Greens")[c(7,3,5)],
  RColorBrewer::brewer.pal(9,"Reds")[2],
  RColorBrewer::brewer.pal(9,"Reds")[6]
)
names(palette_features$Cluster) <- c("1","2","3","4","5","6")

palette_features$ClusterLab <- palette_features$Cluster
# names(palette_features$ClusterLab) <- c("A:PSCs","B:Transition1",
#                                         "C:Transition2","D:Post.EpiLike",
#                                         "E:Mes","F:Spont.Diff.2D")

pal_polo <- RColorBrewer::brewer.pal(9,"Blues")[c(3,5,8)]
names(pal_polo) <- c("Day0","Day2","Day5")
palette_features[["Time"]] <- pal_polo

# Visualization ---
# tSNE ---
tsne_plt_l <- list()
for(feature in metadata_features) {
  tsne_plt_l[[feature]] <- plotsceReducedDim(
    sce = sce,
    dimred = "TSNE",
    col_by = feature,
    shape_by = shape_by,
    point_size = 0.3,
    pal = palette_features[[feature]]
    )
  
  if(feature=="Cluster") {
    tsne_plt_l[[feature]] <- Seurat::LabelClusters(
      plot = tsne_plt_l[[feature]],
      id = feature,
      repel = F,
      col="black"
    )
  }
  
  outfile <- paste0(path_results,"/ptsne_",feature,".pdf")
  pdf(file = outfile, paper = "a4r",w=unit(3,'cm'),h=unit(6,"cm"))
  print(tsne_plt_l[[feature]])
  dev.off()
}
tsne_plt <- patchwork::wrap_plots(tsne_plt_l, nrow = 1)
tsne_plt <- tsne_plt + patchwork::plot_annotation(
  title = "Clustering",
  subtitle = paste0("(n. cells = ",dim(sce)[2],")"),
  theme = theme(
    plot.title = element_text(size=8,hjust = 0.5),
    plot.subtitle = element_text(size=8,hjust = .5)
    )
  )
  
outfile <- paste0(path_results,"/ptsne_",paste0(metadata_features,collapse = "_"),".pdf")
pdf(file = outfile, paper = "a4r",w=unit(2.8*length(metadata_features),'cm'),h=unit(4.5,"cm"))
print(tsne_plt)
dev.off()

# UMAP ---
umap_plt_l <- list()
for(feature in metadata_features) {
  
  umap_plt_l[[feature]] <- plotsceReducedDim(
    sce = sce,
    dimred = "UMAP",
    col_by = feature,
    shape_by = shape_by,
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
  outfile <- paste0(path_results,"/pumap_",feature,".pdf")
  pdf(file = outfile, paper = "a4r",w=unit(3,'cm'),h=unit(6,"cm"))
  print(umap_plt_l[[feature]])
  dev.off()
}

umap_plt <- patchwork::wrap_plots(umap_plt_l, nrow = 1)
umap_plt <- umap_plt + patchwork::plot_annotation(
  title = "Clustering",
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

# Assessing cluster separation ---
# Cluster modularity from graph
ratio <- scran::clusterModularity(g, clust, as.ratio=TRUE)

hm_ratio <- pheatmap::pheatmap(log10(ratio+1), cluster_cols=FALSE, cluster_rows=FALSE,
                               col=rev(heat.colors(100)))

outfile <- paste0(path_results,"/cluster_sep_modularity_heatmap.pdf")
pdf(file = outfile, paper = "a4",w=unit(5,'cm'),h=unit(5,"cm"),useDingbats = F)
hm_ratio
dev.off()
# Bootstrapping
# ass.prob <- bluster::bootstrapStability(sce, FUN=function(x) {
#   g <- scran::buildSNNGraph(x, k=k_snn_graph,use.dimred=dimred_method)
#   igraph::cluster_louvain(g)$membership
# }, clusters=sce$Cluster)
# 
# hm_ass.prob <- pheatmap::pheatmap(ass.prob, cluster_cols=FALSE, cluster_rows=FALSE,
#                                   col=colorRampPalette(c("white", "blue"))(100))
# 
# outfile <- paste0(path_results,"/cluster_sep_bootstrapping_heatmap.pdf")
# pdf(file = outfile, paper = "a4",w=unit(5,'cm'),h=unit(5,"cm"),useDingbats = F)
# hm_ass.prob
# dev.off()

# Cluster composition ---
cluster_composition_feature <- design_var
cluster_composition <- list()
for(i in levels(sce$Cluster)) {
  cluster_composition[[i]] <- table(sce[[cluster_composition_feature]][sce$Cluster==i])
}
pie_cluster_composition <- list()
for(i in names(cluster_composition)) {
  pie_cluster_composition[[i]] <- plot_ggpie(
    as.data.frame(cluster_composition[[i]]),
    variable = cluster_composition_feature,
    pal = palette_features[[cluster_composition_feature]][names(cluster_composition[[i]])],
    title = paste0("Cluster ",i),
    nudge_x_val = 0.5
    )
}

pie_cluster_composition_arranged <- ggpubr::ggarrange(
  plotlist = pie_cluster_composition,
  legend = "left",
  common.legend = F)

outfile <- paste0(path_results,"/pie_cluster_composition_arranged.pdf")
pdf(file = outfile, paper = "a4r",w=unit(8,'cm'),h=unit(8,"cm"),useDingbats = F)
print(pie_cluster_composition_arranged)
dev.off()

# Barplot ---
cluster_composition_df <- reshape2::melt(cluster_composition)
colnames(cluster_composition_df) <- c(design_var,"value","Cluster")

bar_cluster_composition <- ggplot(cluster_composition_df, aes_string(x="Cluster",y="value",fill=design_var)) + 
  geom_bar(position="fill", stat="identity",width=1) + theme_bw() + my_theme +
  scale_fill_manual(values =  palette_features[[cluster_composition_feature]]) + theme(legend.key.size = unit(4,'mm')) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + ylab("Percent of cells") +
  guides(fill=guide_legend(ncol=1))

outfile <- paste0(path_results,"/bar_cluster_composition.pdf")
pdf(file = outfile,paper = "a4r",w = unit(2.6,'cm'),h=unit(3.2,'cm'))
print(bar_cluster_composition)
dev.off()

# 3. Find markers ----
# 4. Save clustered data ----
# Assign cluster labels
palette_features$ClusterLabels <- palette_features$Cluster
names(palette_features$ClusterLabels) <- c(
  "B:Post.EpiLike 1",
  "F:PS-like 3",
  "D:PS-like 1",
  "E:PS-like 2",
  "A:PSCs",
  "C:Post.EpiLike 2"
)

sce$ClusterLabels <- sce$Cluster
levels(sce$ClusterLabels) <- c(
  "B:Post.EpiLike 1",
  "F:PS-like 3",
  "D:PS-like 1",
  "E:PS-like 2",
  "A:PSCs",
  "C:Post.EpiLike 2"
)
sce$ClusterLabels <- factor(
  sce$ClusterLabels,
  levels = sort(levels(sce$ClusterLabels))
  )

palette_features$ClusterNum <- palette_features$ClusterLabels
names(palette_features$ClusterNum) <- gsub("\\:.*","",names(palette_features$ClusterLabels))
sce$ClusterNum <- sce$ClusterLabels
levels(sce$ClusterNum) <- gsub("\\:.*","",levels(sce$ClusterNum))

# Dataset ---
outfile <- paste0(path_results,"/sce.rds")
message(" -- saving: ",outfile)
saveRDS(sce, file = outfile)

# Additional data ---
to_save_extended <- list("palette_features"=palette_features)
outfile <- paste0(path_results,"/extended_data.rds")
message(" -- saving: ",outfile)
saveRDS(to_save_extended, file = outfile)

cells_by_cluster <- colData(sce)[,c("Sample","Time","Cluster","ClusterLabels","batch")]
write.table(
  cells_by_cluster,
  file = gzfile(paste0(path_results,"/cells_by_cluster.gz")),
  row.names = F,
  col.names = T,
  sep = "\t",
  quote = F
)
