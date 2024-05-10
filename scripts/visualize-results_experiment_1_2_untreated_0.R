# Visualize results
setwd("/home/epigen/Martello_3D_human_epiblast/Martello_Polo_TGFb/analysis_oscab")
# 0. Resources ----
source("scripts/resources.R")

# paths
path_sce = "results/correct-batches_experiment_1_2_untreated_0/cluster-cells/sce.rds"
path_cds = "results/correct-batches_experiment_1_2_untreated_0/cluster-cells/pseudotime-monocle3/cds.rds"
path_extended_data = "results/correct-batches_experiment_1_2_untreated_0/cluster-cells/extended_data.rds"
path_cluster_markers = "results/correct-batches_experiment_1_2_untreated_0/cluster-cells/markers-one-vs-rest/de.genes.by.group.cluster_labels.pvalue_0.01.xlsx"
path_results = "results/correct-batches_experiment_1_2_untreated_0/visualize-results"

if(!dir.exists(path_results)) dir.create(path_results, recursive = T)

s33d = 1991

# 1. Clustering ----
sce <- readRDS(path_sce)
extended_data <- readRDS(path_extended_data)
palette_features <- extended_data$palette_features

metadata_features = c("batch","Time","ClusterNum")
shape_by = NULL
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
  
  if(feature=="ClusterNum") {
    tsne_plt_l[[feature]] <- Seurat::LabelClusters(
      plot = tsne_plt_l[[feature]],
      id = feature,
      repel = F,
      col="black"
    )
  }
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
  
  if(feature=="ClusterNum") {
    umap_plt_l[[feature]] <- Seurat::LabelClusters(
      plot = umap_plt_l[[feature]],
      id = feature,
      repel = F,
      col="black"
    )
  }
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

# Cluster composition ---
cluster_composition_feature = "Time"
cluster_composition <- list()
for(i in levels(sce$ClusterNum)) {
  cluster_composition[[i]] <- table(sce[[cluster_composition_feature]][sce$ClusterNum==i])
}
# Barplot ---
cluster_composition_df <- reshape2::melt(cluster_composition)
colnames(cluster_composition_df) <- c("Time","value","ClusterNum")

bar_cluster_composition <- ggplot(cluster_composition_df, aes(x=ClusterNum,y=value,fill=Time)) + 
  geom_bar(position="fill", stat="identity",width=1) + theme_bw() + my_theme +
  scale_fill_manual(values =  palette_features[[cluster_composition_feature]]) + theme(legend.key.size = unit(4,'mm')) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + ylab("Percent of cells") +
  guides(fill=guide_legend(ncol=1))

outfile <- paste0(path_results,"/bar_cluster_composition.pdf")
pdf(file = outfile,paper = "a4r",w = unit(2.5,'cm'),h=unit(3.3,'cm'))
print(bar_cluster_composition)
dev.off()

# 2. Cluster markers ----
cluster_nm <- openxlsx::getSheetNames(path_cluster_markers)
cluster_markers <- lapply(
  cluster_nm, 
  function(i) openxlsx::read.xlsx(path_cluster_markers, sheet = i)
  )
names(cluster_markers) <- cluster_nm
# names(cluster_markers) <- colData(sce)$ClusterLabels[match(cluster_nm,colData(sce)$Cluster)]
# cluster_markers <- cluster_markers[order(names(cluster_markers))]

# to_save <- cluster_markers
# names(to_save) <- gsub("\\:","\\-",names(to_save))
# saveXLSresEdgeR2(
#   res = to_save,
#   outfile = gsub(".pvalue_adj",".cluster_labels.pvalue_adj",path_cluster_markers)
#   )
# rm(to_save)

gene_sets_l <- list(
  "Pluripotency" = c("NANOG", "ZNF398", "POU5F1", "PRDM14", "SOX2", "ZIC2", "OTX2"),
  "Neural" = c("SIX3", "PAX3", "PAX6", "ENC1", "GBX2", "NNAT", "NESTIN", "SOX1", "CHD1"),
  "Epithelial" = c("EPCAM", "CDH1", "ESRP1"),
  "Polarity"=  c("CLDN6","CLDN7","CLDN10","PODXL", "PARD3", "EZR", "OCLN"),
  "Mesenchymal" = c("S100A4", "VIM", "CDH2", "ZEB2"),
  "Primitive streak" = c("GATA6", "MIXL1", "VIM", "EOMES", "CDH2", "CER1", "SOX4", "LIN28B"),
  "TGF-Beta targets" = c("SKIL", "LEFTY1", "LEFTY2", "SMAD7"),
  "ICM" = c("SMAD5", "ESRRB", "GRHL2", "DHRS3", "TRERF1"),
  "preEpi"= c("KLF4", "KLF17", "TFCP2L1", "DPPA2", "DPPA4"),
  "Hypoblast" = c("PDGFRA", "SOX17", "COL4A1", "NID2", "HNF4A")
)

markers_gene_sets_l <- lapply(
  gene_sets_l,
  function(gene_set) {
    markers_gene_sets <- lapply(names(cluster_markers), function(x) {
      res <- gene_set[gene_set%in%cluster_markers[[x]]$gene_name]
      data.frame(
        "gene_name"= res,
        "marker_of_cluster" = rep(x,length(res))
      )
    })
    markers_gene_sets <- do.call(rbind,markers_gene_sets)
  })
markers_gene_sets_l <- do.call(rbind,markers_gene_sets_l)
markers_gene_sets_l$type <- gsub("\\..*","",rownames(markers_gene_sets_l))
rownames(markers_gene_sets_l) <- NULL

markers_gene_sets_l <- markers_gene_sets_l[order(markers_gene_sets_l$marker_of_cluster),]

# selection ---
gene_sel <- c("SOX2",
              "GATA6","MIXL1","EOMES",
              "ZNF398","POU5F1",
              "CDH1","PODXL",
              "CLDN7")
names(gene_sel) <- markers_gene_sets_l$type[match(gene_sel,markers_gene_sets_l$gene_name)]
markers_gene_sets_sel <- markers_gene_sets_l[markers_gene_sets_l$gene_name %in% gene_sel,]
markers_gene_sets_sel$type <- factor(markers_gene_sets_sel$type,
                                     levels = unique(markers_gene_sets_sel$type))
# plot
expr_mat <- as.data.frame(logcounts(sce))
expr_mat <- expr_mat[rownames(expr_mat) %in% gene_sel,]
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
  markers_gene_sets_sel,
  by = "gene_name"
)

expr_mat$gene_name <- factor(expr_mat$gene_name,
                             levels = unique(expr_mat$gene_name)[c(7,8,9,1,2,6,3,4,5)])

markers_plt <- ggplot(
  expr_mat,
  aes(
    x=ClusterNum,
    y=expression,
    col=ClusterNum
    )
) + facet_wrap(type~gene_name,scales = "free_x",nrow=1) +
  # geom_violin(scale = "width") +
  geom_boxplot(outlier.shape = NA,show.legend = F) +
  ggrastr::rasterise(geom_jitter(size=0.4,width = 0.3, show.legend = F),dpi=300) +
  scale_color_manual(values = palette_features$ClusterNum) +
  scale_y_continuous(labels = scales::number_format(accuracy = 1)
                     # ,limits = c(0,10)
                     ) +
  theme_bw() + my_theme 

outfile <- paste0(path_results,"/markers_plt.pdf")
message(" -- saving to: ", outfile)
pdf(file = outfile, paper = "a4r", width = unit(8,'cm'), height = unit(5.5,'cm'))
print(markers_plt)
dev.off()
