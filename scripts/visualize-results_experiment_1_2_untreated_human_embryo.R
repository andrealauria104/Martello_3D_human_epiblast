# Visualize results
setwd("/home/epigen/Martello_3D_human_epiblast/Martello_Polo_TGFb/analysis_oscab")
# 0. Resources ----
source("scripts/resources.R")

# paths
path_sce = "results/correct-batches_experiment_1_2_untreated_human_embryo/cluster-cells/sce.rds"
path_results = "results/correct-batches_experiment_1_2_untreated_human_embryo/visualize-results"
path_cluster_markers = "results/correct-batches_experiment_1_2_untreated_human_embryo/cluster-cells/markers-one-vs-rest/de.genes.by.group.pvalue_adj_0.05.xlsx"

if(!dir.exists(path_results)) dir.create(path_results, recursive = T)

s33d = 1991

# 1. Integration ----
sce <- readRDS(path_sce)

sce$data <- "epiblast_model"
sce$data[sce$batch=="human_embryo"] <- "embryo"
sce$time_epiblast_model <- "Embryo"
sce$time_epiblast_model[sce$data=="epiblast_model"] <- as.character(sce$Time[sce$data=="epiblast_model"])
sce$time_embryo <- "Epiblast model"
sce$time_embryo[sce$data=="embryo"] <- as.character(sce$Time[sce$data=="embryo"])

to.sort <- as.character(sce$time_embryo)
to.sort <- factor(to.sort, levels = unique(to.sort))
sample.ord <- sce$Sample[order(to.sort)]
# palette ---
palette_features <- list()

pal_polo <- c(
  RColorBrewer::brewer.pal(9,"Blues")[c(3,5,7,9)],
  "darkgrey"
  )
names(pal_polo) <- c("DAY0","DAY2","DAY4","DAY6","Embryo")

palette_features[["time_epiblast_model"]]<-pal_polo
sce$time_epiblast_model <- factor(sce$time_epiblast_model, levels = names(pal_polo))

pal_embryo <- c(
  RColorBrewer::brewer.pal(9,"Purples")[c(4,6,9)],
  RColorBrewer::brewer.pal(8,"Dark2")[c(1,4)],
  RColorBrewer::brewer.pal(9,"Oranges")[c(3,6,8)],
  "darkgrey"
)
names(pal_embryo) <- c(
  "Pre-EPIs", "Post-EPIs","PSA-EPI",
  "ICM","Hypoblast","CTBs","STBs","EVTs",
  "Epiblast model"
)
palette_features[["time_embryo"]]<-pal_embryo
sce$time_embryo <- factor(sce$time_embryo, levels = names(pal_embryo))
palette_features[["Cluster"]] <- c(
  RColorBrewer::brewer.pal(7,"Set2")[-c(3,5)],
  RColorBrewer::brewer.pal(8,"Set3")[-c(1:5)],
  RColorBrewer::brewer.pal(8,"Set1")[-c(1:2)]
  )
# t-SNE ---
tsne_plt_l <- list()
for(feature in c("Cluster","time_embryo","time_epiblast_model")) {
  p_title <- gsub("time_","",feature)
  p_title <- gsub("\\_"," ",p_title)
  p_title <- gsub("^(\\w)","\\U\\1", p_title, perl=TRUE)
  if(feature!="Cluster") {
    n_cells <- table(sce$data==gsub("time_","",feature))["TRUE"]
    p_subtitle <- paste0("(n. cells = ",n_cells,")")
  } else {
    p_subtitle <- NULL
  }
  
  tsne_plt_l[[feature]] <- plotsceReducedDim(
    sce = sce,
    dimred = "TSNE",
    sample_ord = sample.ord,
    col_by = feature,
    point_size = 0.4,
    pal = palette_features[[feature]]
  ) +
    guides(col = guide_legend(title = NULL, override.aes = list(size=3))) +
    ggtitle(p_title, subtitle = p_subtitle) +
    theme(plot.subtitle = element_text(size=7,hjust = .5))
  
  if(feature=="Cluster") {
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
  title = "Transcriptome integration",
  theme = theme(
    plot.title = element_text(size=8,hjust = 0.5)
  )
)

outfile <- paste0(path_results,"/tsne_integration.pdf")
pdf(file = outfile, paper = "a4r",w=unit(8.5,'cm'),h=unit(5,"cm"))
print(tsne_plt)
dev.off()

# 2. Integration colored by cell types ----
sce_untreated_0 <- readRDS("results/correct-batches_experiment_1_2_untreated_0/cluster-cells/sce.rds")
extended_data_0 <- readRDS("results/correct-batches_experiment_1_2_untreated_0/cluster-cells/extended_data.rds")

sce$time_epiblast_model_cell_type <- sce$time_epiblast_model
sce$time_epiblast_model_cell_type <- as.character(sce_untreated_0$ClusterNum[match(sce$Sample, sce_untreated_0$Sample)])
sce$time_epiblast_model_cell_type[is.na(sce$time_epiblast_model_cell_type)] <- "Embryo"
sce$time_epiblast_model_cell_type <- factor(sce$time_epiblast_model_cell_type)
palette_features[["time_epiblast_model_cell_type"]] <- c(extended_data_0$palette_features$ClusterNum,
                                                         palette_features$time_epiblast_model[5])

feature = "time_epiblast_model_cell_type"
p_title <- "3D Epiblast model"
n_cells <- table(sce$data)["epiblast_model"]
p_subtitle <- paste0("(n. cells = ",n_cells,")")

tsne_plt_l[[feature]] <- plotsceReducedDim(
  sce = sce,
  dimred = "TSNE",
  sample_ord = sample.ord,
  col_by = feature,
  point_size = 0.4,
  pal = palette_features[[feature]]
) +
  guides(col = guide_legend(title = NULL, override.aes = list(size=3))) +
  ggtitle(p_title, subtitle = p_subtitle) +
  theme(plot.subtitle = element_text(size=7,hjust = .5))


tsne_cell_types_plt <- patchwork::wrap_plots(tsne_plt_l[c(1,2,4)], nrow = 1)
tsne_cell_types_plt <- tsne_cell_types_plt + patchwork::plot_annotation(
  title = "Transcriptome integration",
  theme = theme(
    plot.title = element_text(size=8,hjust = 0.5)
  )
)

outfile <- paste0(path_results,"/tsne_integration_cell_types.pdf")
pdf(file = outfile, paper = "a4r",w=unit(8.5,'cm'),h=unit(5,"cm"))
print(tsne_cell_types_plt)
dev.off()

# 3. Cluster composition ----
shared_clusters_idx <- rowSums(table(sce$Cluster,gsub("_1|_2","_polo",sce$batch))>0)==2
cluster_type <- list(
  "Shared clusters"=names(which(shared_clusters_idx)),
  "Non shared clusters"=names(which(!shared_clusters_idx))
)

# Shared clusters ---
cluster_composition_shared_l <- lapply(
  c("time_embryo","time_epiblast_model","time_epiblast_model_cell_type"),
  function(design_var) {
    cluster_composition_feature <- design_var
    cluster_composition <- list()
    for(i in cluster_type[["Shared clusters"]]) {
      cluster_composition[[as.character(i)]] <- table(sce[[cluster_composition_feature]][sce$Cluster==i])
    }
    
    cluster_composition_df <- reshape2::melt(cluster_composition)
    colnames(cluster_composition_df) <- c("Time","value","Cluster")
    p_title <- gsub("time_","",design_var)
    p_title <- gsub("\\_"," ",p_title)
    p_title <- gsub("^(\\w)","\\U\\1", p_title, perl=TRUE)
    cluster_composition_df$time_type <- p_title
    cluster_composition_df$Cluster <- factor(
      cluster_composition_df$Cluster,
      levels = stringr::str_sort(unique(cluster_composition_df$Cluster),numeric = T)
    )
    return(cluster_composition_df)
    
  }
)
cluster_composition_shared <- do.call(rbind,cluster_composition_shared_l)
cluster_composition_shared$Time <- factor(
  cluster_composition_shared$Time,
  levels = levels(cluster_composition_shared$Time)[c(1:9,14,10:13,15:19)]
)
bar_cluster_composition_shared <- ggplot(cluster_composition_shared, aes(x=Cluster,y=value,fill=Time)) + 
  geom_bar(position="fill", stat="identity",width=1) + theme_bw() + my_theme +
  scale_fill_manual(values =  c(palette_features[[1]],palette_features[[2]],palette_features[[4]])) + 
  theme(
    legend.key.height = unit(1,'mm'),
    legend.key.width  = unit(3,'mm'),
    aspect.ratio = NULL
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + ylab("Percent of cells") +
  guides(fill=guide_legend(title = NULL, ncol=2,nrow=10)) +
  facet_grid(~time_type) +
  ggtitle("Shared clusters")

# Non shared clusters
cluster_composition_nonshared_l <- lapply(
  c("time_embryo"),
  function(design_var) {
    cluster_composition_feature <- design_var
    cluster_composition <- list()
    for(i in cluster_type[["Non shared clusters"]]) {
      cluster_composition[[as.character(i)]] <- table(sce[[cluster_composition_feature]][sce$Cluster==i])
    }
    
    cluster_composition_df <- reshape2::melt(cluster_composition)
    colnames(cluster_composition_df) <- c("Time","value","Cluster")
    p_title <- gsub("time_","",design_var)
    p_title <- gsub("\\_"," ",p_title)
    p_title <- gsub("^(\\w)","\\U\\1", p_title, perl=TRUE)
    cluster_composition_df$time_type <- p_title
    cluster_composition_df$Cluster <- factor(
      cluster_composition_df$Cluster,
      levels = stringr::str_sort(unique(cluster_composition_df$Cluster),numeric = T)
    )
    return(cluster_composition_df)
    
  }
)
cluster_composition_nonshared <- do.call(rbind,cluster_composition_nonshared_l)
bar_cluster_composition_nonshared <- ggplot(cluster_composition_nonshared, aes(x=Cluster,y=value,fill=Time)) + 
  geom_bar(position="fill", stat="identity",width=1) + theme_bw() + my_theme +
  scale_fill_manual(values =  c(palette_features[[1]],palette_features[[2]])) + 
  theme(
    legend.key.height = unit(1,'mm'),
    legend.key.width  = unit(3,'mm'),
    aspect.ratio = NULL
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + ylab("Percent of cells") +
  guides(fill=guide_legend(title = NULL, ncol=2,nrow=9)) +
  facet_grid(~time_type) +
  ggtitle("Non shared clusters")


bar_cluster_composition <- patchwork::wrap_plots(
  bar_cluster_composition_nonshared,
  bar_cluster_composition_shared,
  nrow = 1,
  widths =c(1,2)
)
bar_cluster_composition <- bar_cluster_composition + patchwork::plot_annotation(
  title = "Cluster composition",
  theme = theme(
    plot.title = element_text(size=8,hjust = 0.5)
  )
) + patchwork::plot_layout(widths = c(1, 3.2))

outfile <- paste0(path_results,"/bar_cluster_composition.pdf")
pdf(file = outfile,paper = "a4r",w = unit(8,'cm'),h=unit(2.3,'cm'))
print(bar_cluster_composition)
dev.off()

# 4. Markers expression ----
cluster_nm <- openxlsx::getSheetNames(path_cluster_markers)
cluster_markers <- lapply(
  cluster_nm, 
  function(i) openxlsx::read.xlsx(path_cluster_markers, sheet = i)
)
names(cluster_markers) <- cluster_nm

gene_sets_l <- list(
  "Pluripotency" = c("NANOG", "ZNF398", "POU5F1", "PRDM14", "SOX2", "ZIC2", "OTX2"),
  # "Neural" = c("SIX3", "PAX3", "PAX6", "ENC1", "GBX2", "NNAT", "NESTIN", "SOX1", "CHD1"),
  "Epithelial" = c("EPCAM", "CDH1", "ESRP1"),
  "Polarity"=  c("CLDN6","CLDN7","CLDN10", "PODXL", "PARD3", "EZR", "OCLN"),
  "Mesenchymal" = c("S100A4", "VIM", "CDH2", "ZEB2"),
  "Primitive streak" = c("GATA6", "MIXL1", "VIM", "EOMES", "CDH2", "CER1", "SOX4", "LIN28B"),
  # "TGF-Beta targets" = c("SKIL", "LEFTY1", "LEFTY2", "SMAD7"),
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
markers_gene_sets_l_unique <- lapply(
  split(markers_gene_sets_l,markers_gene_sets_l$type),
  function(x) unique(x$gene_name)
  )

# Main figure ---
selection_list <- list(
  "Pluripotency" = c(2,6),
  "Polarity" = c(2,5),
  "Primitive streak" = c(4,5)
)
markers_gene_sets_l_unique_sel <- markers_gene_sets_l_unique[names(selection_list)]
for(i in names(markers_gene_sets_l_unique_sel)) {
  markers_gene_sets_l_unique_sel[[i]] <- markers_gene_sets_l_unique_sel[[i]][selection_list[[i]]]
}
# tSNE plots ---
tsne_markers_plt_l <- list()
for(i in names(markers_gene_sets_l_unique_sel)) {
  markers.toplot <- prepareDimRedPltGenes(
    sce,
    genes = markers_gene_sets_l_unique_sel[[i]],
    assay.var = "logcounts",
    scale = T
  )
  tsne_markers_plt_l[[i]] <- plotsceReducedDim(
    sce = sce,
    marker = markers.toplot, 
    dimred = "TSNE",
    point_size = .2,
    point_alpha=.5,
    scale_col_low = "darkblue",
    scale_col_mid = "white",
    scale_col_high = "#FD8D3C",
    guide_colourbar_title = "Scaled expression \n\t[Z-score]"
    # guide_colourbar_title = "Expression"
    ) + blank_theme +
    ggtitle(i)
  
}

# tsne_markers_plt <- patchwork::wrap_plots(tsne_markers_plt_l, ncol = 1)
# tsne_markers_plt <- tsne_markers_plt + patchwork::plot_layout(guides = "collect")

tsne_markers_plt <- ggpubr::ggarrange(
  plotlist = tsne_markers_plt_l,
  ncol = 1,
  align = "v",
  common.legend=T,
  legend = "right"
  )

outfile <- paste0(path_results,"/tsne_markers.pdf")
pdf(file = outfile, paper = "a4",width = unit(8,'cm'), height = unit(5,'cm'))
print(tsne_markers_plt)
dev.off()

# Supplementary figure ---
selection_list <- list(
  "ICM" = c(3),
  "preEpi" = c(3,4),
  "Polarity" = c(4,6),
  "Pluripotency"=c(3,7),
  "Primitive streak" = c(6,7,3)
)
markers_gene_sets_l_unique_sel <- markers_gene_sets_l_unique[names(selection_list)]
for(i in names(markers_gene_sets_l_unique_sel)) {
  markers_gene_sets_l_unique_sel[[i]] <- markers_gene_sets_l_unique_sel[[i]][selection_list[[i]]]
}
markers_gene_sets_l_unique_sel <- c(
  list("ICM/preEpi/Pluripotency"=c(
    markers_gene_sets_l_unique_sel[[1]],
    markers_gene_sets_l_unique_sel[[2]],
    markers_gene_sets_l_unique_sel[[4]])),
  list("Polarity/Primitive streak"=c(markers_gene_sets_l_unique_sel[[3]],
       markers_gene_sets_l_unique_sel[[5]]))
  )

# tSNE plots ---
tsne_markers_plt_l <- list()
for(i in names(markers_gene_sets_l_unique_sel)) {
  markers.toplot <- prepareDimRedPltGenes(
    sce,
    genes = markers_gene_sets_l_unique_sel[[i]],
    assay.var = "logcounts",
    scale = T
  )
  
  tsne_markers_plt_l[[i]] <- plotsceReducedDim(
    sce = sce,
    marker = markers.toplot,
    dimred = "TSNE",
    point_size = .2,
    point_alpha=.5,
    scale_col_low = "darkblue",
    scale_col_mid = "white",
    scale_col_high = "#FD8D3C",
    guide_colourbar_title = "Scaled expression \n\t[Z-score]",
    facet_ncol = 5
    # guide_colourbar_title = "Expression"
  ) + blank_theme +
    ggtitle(i)
  
}

# tsne_markers_plt <- patchwork::wrap_plots(tsne_markers_plt_l, ncol = 1)
# tsne_markers_plt <- tsne_markers_plt + patchwork::plot_layout(guides = "collect")

tsne_markers_plt <- ggpubr::ggarrange(
  plotlist = tsne_markers_plt_l,
  nrow = 2,
  align = "v",
  common.legend=T,
  legend = "right"
)

outfile <- paste0(path_results,"/tsne_markers_suppl.pdf")
pdf(file = outfile, paper = "a4",width = unit(12,'cm'), height = unit(4.5,'cm'))
print(tsne_markers_plt)
dev.off()
