# Analyze scRNAseq - Bioconductor workflow
# Integrate scRNAseq datasets
# 0. Resources ----
source("scripts/resources.R")

# paths ---
paths_sce <- c("results/experiment_1/sce.rds","results/experiment_2_untreated/sce.rds")
names(paths_sce) <- c("experiment_1","experiment_2")

path_results  <- "results/correct-batches_experiment_1_2_untreated"
if(!dir.exists(path_results)) dir.create(path_results, recursive = T)

# Global vars ---
metadata_features <- c("batch","Time","Plate")
n_top_hvgs = 5000
perplexity = 20
ndim_corrected = 50
s33d = 1991

# 1. Read datasets ----
cells_by_cluster <- read.delim("results/correct-batches_experiment_1_2_untreated_0/cluster-cells/cells_by_cluster.gz")
sample_list <- cells_by_cluster$Sample[cells_by_cluster$ClusterLabels=="F:Spont.Diff.2D"]
sce_datasets <- lapply(paths_sce, readRDS)
gene_universe <- Reduce(intersect,lapply(sce_datasets,function(x) rownames(x)))
sce_datasets  <- lapply(sce_datasets, function(x) x[gene_universe,])
sce_datasets  <- lapply(sce_datasets, function(x) x[,!colnames(x)%in%sample_list])
# sce_datasets <- lapply(
#   sce_datasets,
#   function(sce) {
#     clust <- scran::quickCluster(sce)
#     # rowData(sce)$use_for_norm <- !(grepl("rRNA_",rowData(sce)$gene_type) | grepl("^MT-",rowData(sce)$gene_name))
#     # sce <- scran::computeSumFactors(sce, subset.row=rowData(sce)$use_for_norm)
#     sce <- scran::computeSumFactors(sce, clusters=clust, min.mean=0.1)
#     sce <- scater::logNormCounts(sce)
#   }
# )
# Palettes for metadata ---
metadata_features <- get_comma_arglist(metadata_features)
palette_features  <- get_palette_features(metadata_features)

pal_polo <- RColorBrewer::brewer.pal(9,"Blues")[c(3,5,7,9)]
names(pal_polo) <- c("DAY0","DAY2","DAY4","DAY6")
palette_features[["Time"]] <- pal_polo

# 2. Variance modelling ----
dec_datasets <- lapply(sce_datasets, modelGeneVar)
combined.dec <- do.call(combineVar, dec_datasets)
chosen.hvgs  <- getTopHVGs(combined.dec, n=n_top_hvgs)

# 3. Combine datasets -----
sce <- batchelor::correctExperiments(sce_datasets, PARAM=NoCorrectParam())
sce <- batchelor::multiBatchNorm(
  sce,
  batch = sce$batch,norm.args=list(use_altexps=FALSE)
  )

# 4. Explore combined dataset - Dimensionality reduction ----
# PCA ---
sce <- scater::runPCA(
  sce,
  ncomponents=100,
  subset_row=chosen.hvgs,
  exprs_values="logcounts",
  BSPARAM=BiocSingular::ExactParam()
  )

# Number of PCs - Elbow method ---
percent.var <- attr(reducedDim(sce), "percentVar")
chosen.elbow <- PCAtools::findElbowPoint(percent.var)
reducedDim(sce, "PCA.elbow") <- reducedDim(sce)[,1:chosen.elbow]
reducedDimNames(sce)

# Number of PCs - Population structure: clustered PCs ---
# pcs <- reducedDim(sce, "PCA")
# choices <- getClusteredPCs(pcs)
# val <- metadata(choices)$chosen
# reducedDim(sce, "PCA.clust") <- pcs[,1:val]

pca_plt_uncorrected_l <- list()
for(feature in metadata_features) {
  pca_plt_uncorrected_l[[feature]] <- plotsceReducedDim(
    sce = sce,
    dimred = "PCA",
    col_by = feature,
    point_size = 0.3,
    pal = palette_features[[feature]]
    )
}

pca_plt_uncorrected <- patchwork::wrap_plots(pca_plt_uncorrected_l, nrow = 1)

outfile <- paste0(path_results,"/pca_uncorrected_",paste0(metadata_features,collapse = "_"),".pdf")
pdf(file = outfile, paper = "a4r",w=unit(3*length(metadata_features),'cm'),h=unit(6,"cm"))
print(pca_plt_uncorrected)
dev.off()

dimred_method = "PCA.elbow"
# t-SNE ---
set.seed(s33d)
sce <- scater::runTSNE(
  sce,
  dimred=dimred_method,
  perplexity=perplexity,
  name = "TSNE_uncorrected"
  )

tsne_plt_uncorrected_l <- list()
for(feature in metadata_features) {
  tsne_plt_uncorrected_l[[feature]] <- plotsceReducedDim(
    sce = sce,
    dimred = "TSNE_uncorrected",
    col_by = feature,
    point_size = 0.3,
    pal = palette_features[[feature]]
  ) + xlab("tSNE1" ) + ylab("tSNE2")
}
tsne_plt_uncorrected <- patchwork::wrap_plots(tsne_plt_uncorrected_l, nrow = 1)
tsne_plt_uncorrected <- tsne_plt_uncorrected + patchwork::plot_annotation(
  title = "Uncorrected",
  theme = theme(plot.title = element_text(size=8,hjust = 0.5))
)

outfile <- paste0(path_results,"/tsne_uncorrected_",paste0(metadata_features,collapse = "_"),".pdf")
pdf(file = outfile, paper = "a4r",w=unit(3*length(metadata_features),'cm'),h=unit(6,"cm"))
print(tsne_plt_uncorrected)
dev.off()

# UMAP ---
set.seed(s33d)
sce <- scater::runUMAP(
  sce,
  dimred=dimred_method,
  name = "UMAP_uncorrected",
  )

umap_plt_uncorrected_l <- list()
for(feature in metadata_features) {
  umap_plt_uncorrected_l[[feature]] <- plotsceReducedDim(
    sce = sce,
    dimred = "UMAP_uncorrected",
    col_by = feature,
    point_size = 0.3,
    pal = palette_features[[feature]]
  ) + xlab("UMAP1" ) + ylab("UMAP2")
}

umap_plt_uncorrected <- patchwork::wrap_plots(umap_plt_uncorrected_l, nrow = 1)
umap_plt_uncorrected <- umap_plt_uncorrected + patchwork::plot_annotation(
  title = "Uncorrected",
  theme = theme(plot.title = element_text(size=8,hjust = 0.5))
)

outfile <- paste0(path_results,"/umap_uncorrected_",paste0(metadata_features,collapse = "_"),".pdf")
pdf(file = outfile, paper = "a4r",w=unit(3*length(metadata_features),'cm'),h=unit(6,"cm"))
print(umap_plt_uncorrected)
dev.off()

# 5. Batch correction ----
# fastMNN
set.seed(s33d)
f.out <- batchelor::fastMNN(
  sce,
  batch=sce$batch,
  subset.row=chosen.hvgs,
  correct.all=TRUE
)
reducedDim(sce, "corrected")  <- reducedDim(f.out, "corrected")
reducedDim(sce, "corrected_sel")  <- reducedDim(f.out, "corrected")[,1:ndim_corrected]
assay(sce, "reconstructed", withDimnames = F) <- as.matrix(assay(f.out, "reconstructed"))

# Dimensionality reduction 
set.seed(s33d)
sce <- scater::runTSNE(
  sce,
  dimred="corrected_sel",
  perplexity=perplexity
  )

tsne_plt_l <- list()
for(feature in metadata_features) {
  tsne_plt_l[[feature]] <- plotsceReducedDim(
    sce = sce,
    dimred = "TSNE",
    col_by = feature,
    point_size = 0.3,
    pal = palette_features[[feature]]
  )
}
tsne_plt <- patchwork::wrap_plots(tsne_plt_l, nrow = 1)
tsne_plt <- tsne_plt + patchwork::plot_annotation(
  title = "Corrected [MNN]",
  subtitle = paste0("(n. cells = ",dim(sce)[2],")"),
  theme = theme(
    plot.title = element_text(size=8,hjust = 0.5),
    plot.subtitle = element_text(size=8,hjust = .5)
  )
)

outfile <- paste0(path_results,"/tsne_corrected_",paste0(metadata_features,collapse = "_"),".pdf")
pdf(file = outfile, paper = "a4r",w=unit(3*length(metadata_features),'cm'),h=unit(6.3,"cm"))
print(tsne_plt)
dev.off()

tsne_plt_correction <- ggpubr::ggarrange(
  tsne_plt_uncorrected,
  tsne_plt,
  align = "v",
  nrow = 2
)

outfile <- paste0(path_results,"/tsne_correction.pdf")
pdf(file = outfile, paper = "a4r",w=unit(6*length(metadata_features),'cm'),h=unit(12,"cm"))
print(tsne_plt_correction)
dev.off()

# UMAP ---
set.seed(s33d)
sce <- scater::runUMAP(
  sce,
  dimred="corrected",
  metric = "euclidean",
  min_dist=0.05,
  n_neighbors=15,
  spread=1.5
)

umap_plt_l <- list()
for(feature in metadata_features) {
  umap_plt_l[[feature]] <- plotsceReducedDim(
    sce = sce,
    dimred = "UMAP",
    col_by = feature,
    point_size = 0.3,
    pal = palette_features[[feature]]
  )
}

umap_plt <- patchwork::wrap_plots(umap_plt_l, nrow = 1)
umap_plt <- umap_plt + patchwork::plot_annotation(
  title = "Corrected [MNN]",
  subtitle = paste0("(n. cells = ",dim(sce)[2],")"),
  theme = theme(
    plot.title = element_text(size=8,hjust = 0.5),
    plot.subtitle = element_text(size=8,hjust = .5)
    )
)

outfile <- paste0(path_results,"/umap_corrected_",paste0(metadata_features,collapse = "_"),".pdf")
pdf(file = outfile, paper = "a4r",w=unit(3*length(metadata_features),'cm'),h=unit(6.3,"cm"))
print(umap_plt)
dev.off()

umap_plt_correction <- ggpubr::ggarrange(
  umap_plt_uncorrected,
  umap_plt,
  align = "v",
  nrow = 2
)

outfile <- paste0(path_results,"/umap_correction.pdf")
pdf(file = outfile, paper = "a4r",w=unit(6*length(metadata_features),'cm'),h=unit(12,"cm"))
print(umap_plt_correction)
dev.off()

# 6. Save corrected data ----
outfile <- paste0(path_results,"/sce.rds")
message(" -- saving: ",outfile)
saveRDS(sce, file = outfile)
