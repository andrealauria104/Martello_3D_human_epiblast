

cds <- readRDS("results/correct-batches_experiment_1_2_untreated/cluster-cells/pseudotime-monocle3/cds.rds")
extended_data <- readRDS("results/correct-batches_experiment_1_2_untreated/cluster-cells/extended_data.rds")

gene_sel <- "POU5F1"
expr_mat <- as.data.frame(logcounts(cds))
cells <- colData(cds)$Sample[colData(cds)$Branch_2 != "other"]
expr_mat <- expr_mat[rownames(expr_mat) %in% gene_sel,colnames(expr_mat) %in% cells]
expr_mat <- cbind("gene_name"=rownames(expr_mat),expr_mat)
expr_mat <- reshape2::melt(
  expr_mat,
  id.var="gene_name",
  variable.name = "Sample",
  value.name = "expression"
)
expr_mat <- as.data.frame(merge(expr_mat,colData(cds),by="Sample"))

ggplot(
  expr_mat,
  aes(
    x=Pseudotime_scaled,
    y=expression,
    col=Cluster
  )
) +
  geom_point(size=.4) +
  scale_color_manual(
    values = extended_data$palette_features$Cluster
      ) +
  theme_bw() + my_theme


genes <- setdiff(rownames(pr_graph_test_res_1_sig),rownames(pr_graph_test_res_2_sig))
genes <- c("MIXL1","EOMES")
c_nm <- colData(cds)$Sample[colData(cds)$Branch_1=="PSC_to_PS-like"]
pt.matrix <- as.matrix(cds@assays@data$logcounts[match(genes,rowData(cds)$gene_short_name),c_nm])
c_ord_idx <- order(subset(colData(cds),Branch_1=="PSC_to_PS-like")$Pseudotime)
pt.matrix <- pt.matrix[,c_ord_idx]

pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- genes


library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

ht_1 <- Heatmap(
  pt.matrix,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "RdBu"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 6),
  row_title_rot                = 0,
  cluster_rows                 = T,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE,
  use_raster = T,
  show_column_dend = F,
  show_row_dend = F,
  # left_annotation = rowAnnotation(mark=anno),
  # top_annotation = ha_column,
  column_title = "PSC_to_PS-like",
  column_title_gp = gpar(fontsize=8),
  clustering_distance_rows = "euclidean",
  clustering_method_rows = "ward.D2",
  border = T)
draw(ht_1)
