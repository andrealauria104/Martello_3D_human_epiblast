# Visualize results
setwd("/home/epigen/Martello_3D_human_epiblast/Martello_Polo_MIXL1/run240820_run241003/analysis_oscab")
# 0. Resources ----
source("scripts/resources.R")

# paths
path_sce = "results/cluster-cells/sce.rds"
path_extended_data = "results/cluster-cells/extended_data.rds"
path_results = "results//visualize-results-post"
path_de_markers = "results/cluster-cells/markers-one-vs-rest/de.genes.by.group.rds"

if(!dir.exists(path_results)) dir.create(path_results, recursive = T)

s33d = 1991

metadata_features <- c("Time","Plate","batch")

# 1. Markers ----
sce <- readRDS(path_sce)
extended_data <- readRDS(path_extended_data)
palette_features <- extended_data$palette_features

gene_sets_l <- list(
  "Pluripotency" = c("NANOG", "ZNF398", "POU5F1", "PRDM14", "SOX2", "ZIC2", "OTX2"),
  # "Neural" = c("SIX3", "PAX3", "PAX6", "ENC1", "GBX2", "NNAT", "NESTIN", "SOX1", "CHD1"),
  "Epithelial" = c("EPCAM", "CDH1", "ESRP1","ESRP2","CDH3","CLDN3"),
  "Polarity"=  c("CLDN6","CLDN7","CLDN10","PODXL", "PARD3", "EZR", "OCLN"),
  "Mesenchymal" = c("S100A4", "VIM", "CDH2", "ZEB2","ZEB1","TWIST1"),
  "Primitive streak" = c("GATA6", "MIXL1", "VIM", "EOMES", "CDH2", "CER1", "SOX4", "LIN28B","TBXT","SOX17")
  # "TGF-Beta targets" = c("SKIL", "LEFTY1", "LEFTY2", "SMAD7"),
  # "ICM" = c("SMAD5", "ESRRB", "GRHL2", "DHRS3", "TRERF1"),
  # "preEpi"= c("KLF4", "KLF17", "TFCP2L1", "DPPA2", "DPPA4"),
  # "Hypoblast" = c("PDGFRA", "SOX17", "COL4A1", "NID2", "HNF4A")
)

# markers_gene_sets_l <- lapply(
#   gene_sets_l,
#   function(gene_set) {
#     markers_gene_sets <- lapply(names(cluster_markers), function(x) {
#       res <- gene_set[gene_set%in%cluster_markers[[x]]$gene_name]
#       data.frame(
#         "gene_name"= res,
#         "marker_of_cluster" = rep(x,length(res))
#       )
#     })
#     markers_gene_sets <- do.call(rbind,markers_gene_sets)
#   })
markers_gene_sets_l <- reshape2::melt(gene_sets_l)
colnames(markers_gene_sets_l) <- c("gene_name","type")
# markers_gene_sets_l$type <- gsub("\\..*","",rownames(markers_gene_sets_l))
rownames(markers_gene_sets_l) <- NULL

# markers_gene_sets_l <- markers_gene_sets_l[order(markers_gene_sets_l$marker_of_cluster),]

# selection ---
gene_sel <- c("SOX2","PRDM14","OCLN","ESRP1","EPCAM","CLDN3",
              "GATA6","MIXL1","EOMES","TBXT","TWIST1","SOX17",
              "ZNF398","POU5F1","NANOG","ZEB1","ZEB2","CER1",
              "CDH1","PODXL","VIM","CDH2","CLDN7")
names(gene_sel) <- markers_gene_sets_l$type[match(gene_sel,markers_gene_sets_l$gene_name)]
markers_gene_sets_sel <- markers_gene_sets_l[markers_gene_sets_l$gene_name %in% gene_sel,]
markers_gene_sets_sel$type <- factor(markers_gene_sets_sel$type,
                                     levels = unique(markers_gene_sets_sel$type))
# plot
expr_mat <- as.data.frame(logcounts(sce))
# expr_mat <- expr_mat[rownames(expr_mat) %in% gene_sel,]
expr_mat <- expr_mat[rownames(expr_mat) %in% unique(markers_gene_sets_l$gene_name),]
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
  markers_gene_sets_l,
  by = "gene_name"
)

# expr_mat$gene_name <- factor(expr_mat$gene_name,
#                              levels = unique(expr_mat$gene_name)[c(7,8,9,1,2,6,3,4,5)])

markers_plt <- ggplot(
  expr_mat,
  aes(
    x=Time,
    y=expression,
    col=Time
    )
) + facet_wrap(type~gene_name,scales = "free_x",nrow=3) +
  # geom_violin(scale = "width") +
  geom_boxplot(outlier.shape = NA,show.legend = F) +
  ggrastr::rasterise(geom_jitter(size=0.25,width = 0.3, show.legend = F),dpi=300) +
  scale_color_manual(values = palette_features$Time) +
  scale_y_continuous(labels = scales::number_format(accuracy = 1)
                     # ,limits = c(0,10)
                     ) +
  theme_bw() + my_theme

outfile <- paste0(path_results,"/markers_plt.pdf")
message(" -- saving to: ", outfile)
pdf(file = outfile, paper = "a4r", width = unit(14,'cm'), height = unit(14,'cm'))
print(markers_plt)
dev.off()

# 2. Dimensionality reduction ----
# UMAP
pumap <- list()
for(feature in metadata_features) {
  
  pumap[[feature]] <- plotsceReducedDim(
    sce = sce,
    dimred = "UMAP",
    col_by = feature,
    point_size = 0.3,
    pal = palette_features[[feature]]
  )
}
pumap_list1 <- ggpubr::ggarrange(plotlist = pumap[2:3], nrow=1, align = "hv")
pumap_list1 <- ggpubr::ggarrange(plotlist = list(pumap$Time,pumap_list1), nrow=1, widths = c(1,1.5))

outfile <- paste0(path_results,"/pumap_",paste0(metadata_features,collapse = "_"),".pdf")
pdf(file = outfile, paper = "a4r",w=unit(9,"cm"),h=unit(6,"cm"))
print(pumap_list1)
dev.off()

feature <- 'ClusterLabels'
outfile <- paste0(path_results,"/pumap_",feature,".pdf")
pdf(file = outfile, paper = "a4r",w=unit(4,"cm"),h=unit(4,"cm"))
print(
  plotsceReducedDim(
    sce = sce,
    dimred = "UMAP",
    col_by = feature,
    point_size = 0.3,
    pal = palette_features[[feature]]
  )
)
dev.off()
# 3. Markers by clusters ----
markers_clust_plt <- ggplot(
  expr_mat,
  aes(
    x=ClusterNum,
    y=expression,
    col=ClusterNum
  )
) + facet_wrap(type~gene_name,scales = "free_x",nrow=5) +
  # geom_violin(scale = "width") +
  geom_boxplot(outlier.shape = NA,show.legend = F) +
  ggrastr::rasterise(geom_jitter(size=0.25,width = 0.3, show.legend = F),dpi=300) +
  scale_color_manual(values = palette_features$ClusterNum) +
  scale_y_continuous(labels = scales::number_format(accuracy = 1)
                     # ,limits = c(0,10)
  ) +
  theme_bw() + my_theme

outfile <- paste0(path_results,"/markers_clust_plt.pdf")
message(" -- saving to: ", outfile)
pdf(file = outfile, paper = "a4", width = unit(7,'cm'), height = unit(15,'cm'))
print(markers_clust_plt)
dev.off()

# 4. Markers summary ----
de_markers <- readRDS(path_de_markers)
de_markers_sel <- lapply(de_markers, function(x) {
 x$gene_name <- rownames(x)
 x <- x[,c("gene_name","pvalue","fdr")]
 merge(x, markers_gene_sets_sel, by="gene_name", all=F)
})
de_markers_sel <- do.call(rbind, de_markers_sel)
de_markers_sel$clust <- gsub("\\.\\d+$","",rownames(de_markers_sel))
rownames(de_markers_sel) <- NULL
de_markers_sel$score <- -log10(de_markers_sel$fdr)
de_markers_sel_mat <- reshape2::dcast(de_markers_sel, clust~gene_name,
                                      value.var = "score", fun.aggregate = mean)
rownames(de_markers_sel_mat) <- de_markers_sel_mat[,1]
de_markers_sel_mat <- de_markers_sel_mat[,-1]
de_markers_sel_mat <- as.matrix(de_markers_sel_mat)

cl_map <-  unique(colData(sce)[,c("Cluster","ClusterNum")])
rownames(de_markers_sel_mat) <- cl_map$ClusterNum[match(rownames(de_markers_sel_mat),cl_map$Cluster)]
de_markers_sel_mat <- de_markers_sel_mat[order(rownames(de_markers_sel_mat)),]
de_markers_sel_mat[is.nan(de_markers_sel_mat)] <- 0
pal_hm <- circlize::colorRamp2(
  c(0,5,15,20),
  c('white',RColorBrewer::brewer.pal(3,"Reds"))
)
fc_hm <- ComplexHeatmap::Heatmap(
  de_markers_sel_mat,
  col = pal_hm,
  show_row_dend = F,
  show_column_dend = F,
  cluster_rows = F,
  border = T,
  column_split = markers_gene_sets_sel$type[match(colnames(de_markers_sel_mat),markers_gene_sets_sel$gene_name)],
  cluster_column_slices = F,
  row_names_gp = gpar(fontsize=8),
  column_names_rot = 45,
  column_names_gp = gpar(fontsize=8),
  row_title_gp = gpar(fontsize=8),
  column_title_gp = gpar(fontsize=8),
  row_names_side = 'left',
  heatmap_legend_param = list(title='-log10(FDR)',
                              title_gp=gpar(fontsize=8,face='plain'),
                              labels_gp=gpar(fontsize=8,face='plain'),
                              legend_direction='horizontal')
)

outfile  <- paste0(path_results, "/pv_hm.pdf")
message('-- saving to: ',outfile)
pdf(file = outfile, paper = "a4r", w = unit(6.2,'cm'), h = unit(1.9,'cm'))
draw(fc_hm, heatmap_legend_side='bottom')
dev.off()

# 5. Markers by clusters - second selection ----
markers_sel2 <- c(
  "NANOG","POU5F1",
  "CDH1","ESRP1","CLDN7",
  "VIM","MIXL1","EOMES","GATA6","SOX17"
)
expr_mat2 <-expr_mat[expr_mat$gene_name%in%markers_sel2,]
expr_mat2$type <- factor(expr_mat2$type, levels = c("Pluripotency","Epithelial",
                                                    "Polarity","Mesenchymal","Primitive streak"))
# expr_mat2 <- expr_mat2[-which(expr_mat2$type=="Primitive streak" & expr_mat2$gene_name=="CDH2"),]
expr_mat2 <- expr_mat2[-which(expr_mat2$type=="Primitive streak" & expr_mat2$gene_name=="VIM"),]
markers_clust2_plt <- ggplot(
  expr_mat2,
  aes(
    x=ClusterNum,
    y=expression,
    col=ClusterNum
  )
) + facet_wrap(type~gene_name,scales = "free_x",nrow=2) +
  # geom_violin(scale = "width") +
  ylab("Expression") +
  geom_boxplot(outlier.shape = NA,show.legend = F, fatten=2) +
  ggrastr::rasterise(geom_jitter(size=0.18,width = 0.1, show.legend = F, alpha=.6),dpi=300) +
  scale_color_manual(values = palette_features$ClusterNum) +
  scale_y_continuous(labels = scales::number_format(accuracy = 1)
                     # ,limits = c(0,10)
  ) +
  theme_bw() + my_theme

outfile <- paste0(path_results,"/markers_clust_second_sel_plt.pdf")
message(" -- saving to: ", outfile)
pdf(file = outfile, paper = "a4", width = unit(4,'cm'), height = unit(6,'cm'))
print(markers_clust2_plt)
dev.off()
