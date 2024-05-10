# Visualize results
setwd("/home/epigen/Martello_3D_human_epiblast/Martello_Polo_TGFb/analysis_oscab")
# 0. Resources ----
source("scripts/resources.R")

# paths
path_sce = "results/correct-batches_experiment_1_2_untreated/cluster-cells/sce.rds"
path_cds = "results/correct-batches_experiment_1_2_untreated/cluster-cells/pseudotime-monocle3/cds.rds"
path_extended_data = "results/correct-batches_experiment_1_2_untreated_0/cluster-cells/extended_data.rds"
path_cluster_markers = "results/correct-batches_experiment_1_2_untreated_0/cluster-cells/markers-one-vs-rest/de.genes.by.group.cluster_labels.pvalue_0.01.xlsx"
path_results = "results/correct-batches_experiment_1_2_untreated/visualize-results"

if(!dir.exists(path_results)) dir.create(path_results, recursive = T)

s33d = 1991

# 1. Cluster markers ----
cluster_nm <- openxlsx::getSheetNames(path_cluster_markers)
cluster_markers <- lapply(
  cluster_nm, 
  function(i) openxlsx::read.xlsx(path_cluster_markers, sheet = i)
)
names(cluster_markers) <- cluster_nm

gene_sets_l <- list(
  "Pluripotency" = c("KLF4","DPPA4","NANOG", "ZNF398", "POU5F1", "PRDM14", "SOX2", "ZIC2", "OTX2"),
  "Neural" = c("SIX3", "PAX3", "PAX6", "ENC1", "GBX2", "NNAT", "NESTIN", "SOX1", "CHD1"),
  "Epithelial" = c("EPCAM", "CDH1", "ESRP1"),
  "Polarity"=  c("CLDN6","CLDN7","CLDN10", "PODXL", "PARD3", "EZR", "OCLN"),
  "Mesenchymal" = c("S100A4", "VIM", "CDH2", "ZEB2"),
  "Primitive streak" = c("GATA6", "MIXL1", "VIM", "EOMES", "CDH2", "CER1", "SOX4", "LIN28B"),
  "TGF-Beta targets" = c("SKIL", "LEFTY1", "LEFTY2", "SMAD7"),
  "ICM" = c("SMAD5", "ESRRB", "GRHL2", "DHRS3", "TRERF1"),
  "preEpi"= c(
    # "KLF4", 
              "KLF17", "TFCP2L1", "DPPA2"
              # , "DPPA4"
              ),
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
levels(markers_gene_sets_sel$type) <- gsub("^Epithelial$","Epithelial\\/Polarity",levels(markers_gene_sets_sel$type))
levels(markers_gene_sets_sel$type) <- gsub("^Polarity$","Epithelial\\/Polarity",levels(markers_gene_sets_sel$type))
levels(markers_gene_sets_sel$type)
# 2. Pseudotime ----
cds <- readRDS(path_cds)
extended_data <- readRDS(path_extended_data)
palette_features <- extended_data$palette_features

expr_mat <- as.data.frame(logcounts(cds))
expr_mat <- expr_mat[rownames(expr_mat) %in% markers_gene_sets_sel$gene_name,]
expr_mat <- as.data.frame(t(scale(t(expr_mat))))
expr_mat <- cbind("gene_name"=rownames(expr_mat),expr_mat)
expr_mat <- reshape2::melt(
  expr_mat,
  id.var="gene_name",
  variable.name = "Sample",
  value.name = "expression"
)
expr_mat <- as.data.frame(merge(expr_mat,colData(cds),by="Sample"))
expr_mat <- merge(
  expr_mat,
  unique(markers_gene_sets_sel[,c(1,3)]),
  by = "gene_name"
)

expr_mat_1 <- expr_mat[expr_mat$ClusterNum%in%c("A","B","C","D"),]
expr_mat_1$Path <- "PSCs to PostEpi-like"
expr_mat_2 <- expr_mat[expr_mat$ClusterNum%in%c("A","B","C","E"),]
expr_mat_2$Path <- "PSCs to PS-like"

expr_mat_l <- rbind(
  expr_mat_1,
  expr_mat_2
)

expr_mat_l$Path <- factor(expr_mat_l$Path, levels = unique(expr_mat_l$Path))
types <- levels(expr_mat_l$type)
paths <- c("PSCs to PostEpi-like","PSCs to PostEpi-like","PSCs to PS-like")
names(paths) <- types
markers_pseudo_plt <- list()
for(i in types) {
  markers_pseudo_plt[[i]] <- ggplot(
    expr_mat_l[expr_mat_l$type==i & expr_mat_l$Path == paths[i],],
    aes(
      x=Pseudotime_scaled,
      y=expression,
      col=ClusterNum
    )
  ) + 
    # facet_grid(Path~gene_name,scales = "free") + 
    facet_grid(~gene_name,scales = "free") + 
    ggrastr::rasterise(geom_point(size=0.4,show.legend = F),dpi=600) +
    scale_color_manual(values = extended_data$palette_features$ClusterNum) +
    geom_smooth(aes(col=NULL),se = T,fullrange=T,show.legend = F) +
    theme_bw() + my_theme + 
    scale_x_continuous(breaks = seq(0,1,.5), labels = scales::number_format(accuracy = 1)) +
    ggtitle(i) + ylab("Relative expression") +
    xlab("Pseudotime") 
}

markers_pseudo_plt_arranged <- ggpubr::ggarrange(
  plotlist = markers_pseudo_plt,
  nrow=1,
  align = "hv")

outfile <- paste0(path_results,"/markers_pseudo_plt.pdf")
message(" -- saving to: ", outfile)
pdf(file = outfile, paper = "a4r", width = unit(8.5,'cm'), height = unit(7,'cm'))
print(markers_pseudo_plt_arranged)
dev.off()

# 3. Clusters in pseudotime ----
to_plot <- as.data.frame(colData(cds))
to_plot <- plyr::ddply(
  to_plot,
  .(ClusterLabels),
  mutate,
  mean_pseudotime_by_clust=mean(Pseudotime)
  )

to_plot$ClusterLabels <- factor(
  to_plot$ClusterLabels,
  levels = unique(to_plot$ClusterLabels[order(to_plot$mean_pseudotime_by_clust)])
  )

clusters_in_pseudotime_plt <- ggplot(to_plot, aes(x = Pseudotime_scaled, y = ClusterLabels, colour = ClusterLabels)) +
  ggrastr::rasterize(geom_jitter(height =0.3, size=.5),dpi=300) +
  theme_bw() + my_theme +
  theme(aspect.ratio = NULL)+
  scale_color_manual(values = palette_features$ClusterLabels)  +
  xlab("Pseudotime [scaled]") + ylab("Cell type") +
  ggtitle("Cells ordered by pseudotime") +
  guides(col=guide_legend(title="Cell type",override.aes = list(size=2)))

outfile <- paste0(path_results,"/clusters_in_pseudotime_plt.pdf")
pdf(file = outfile, paper = "a4r",w=unit(3.8,'cm'),h=unit(1.6,"cm"))
print(clusters_in_pseudotime_plt)
dev.off()
