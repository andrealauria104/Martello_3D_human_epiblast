setwd("/home/epigen/Martello_3D_human_epiblast/public_datasets/Zhao_et_al_Nature_Methods_2025")
library(Seurat)
library(reshape2)
library(ggplot2)
library(plyr)

path_results <- "finalized_scran"

if(!dir.exists(path_results)) dir.create(path_results, recursive = T)
# Read seurat object ----
seurat_obj <- readRDS("psd.R3.6.em.seurat.ob.rds")

sce <- as.SingleCellExperiment(seurat_obj)
sce <- scran::computeSumFactors(sce)
sce <- scater::logNormCounts(sce)
# Define markers panel ----
gene_sets_l <- list(
  "Epithelial" = c("CDH1", "ESRP1", "EPCAM"),
  "Polarity"=  c("CLDN6","CLDN7", "OCLN", "PODXL"),
  "Mesenchymal" = c("CDH2", "ZEB1", "VIM")
)
gene_sets_l <- reshape2::melt(gene_sets_l)
colnames(gene_sets_l) <- c("gene_name", "marker_type")

# markers_expr <- Seurat::GetAssayData(seurat_obj, assay = "RNA", slot = "data")
markers_expr <- as.matrix(assay(sce, 'logcounts'))
markers_expr <- markers_expr[which(rownames(markers_expr) %in% unique(gene_sets_l$gene_name)),]
markers_expr <- reshape2::melt(markers_expr, varnames = c("gene_name", "cell_id"))
metadata <- seurat_obj[[]]
metadata$cell_id <- rownames(metadata)

markers_expr <- merge(markers_expr, metadata, by = "cell_id")
markers_expr <- merge(markers_expr, gene_sets_l, by = "gene_name")
# genes <- unique(unlist(gene_sets_l))
# expr_meta <- FetchData(seurat_obj, vars = c(genes, "cell_id","rename_EML", "sub_rename_EML"))

# Refine annotation ----
markers_expr <- markers_expr[!markers_expr$rename_EML %in% c("Ambiguous", "Unknown"), ]
ct_extraembryonic <- c("Hypoblast", "EPI.PrE.INT", "TE", "CTB", "STB", "EVT", "Amnion","ExE_Mes", "YSE")
markers_expr$annotation <- markers_expr$rename_EML
epi_idx <- which(markers_expr$annotation == "Epiblast")
exe_idx <- which(markers_expr$annotation %in% ct_extraembryonic)
markers_expr$annotation[epi_idx] <- markers_expr$sub_rename_EML[epi_idx]
markers_expr$annotation[exe_idx] <- "Extra-embryonic"
annot_ord <- c("Zygote", "2-4 cell", "8 cell", "Morula", "Prelineage",
               "ICM", "Early_EPI", "Late_EPI", "PriS",
               "Mesoderm", "AdvMes")

markers_expr <- markers_expr[markers_expr$annotation %in% annot_ord, ]

markers_expr$annotation <- factor(
  markers_expr$annotation,
  levels = annot_ord
)
markers_expr$annot_lab <- as.character(markers_expr$annotation)
markers_expr$annot_lab[markers_expr$annotation %in% c("Zygote", "2-4 cell", "8 cell", "Morula", "Prelineage", "ICM")] <- "prelineage"
markers_expr$annot_lab[markers_expr$annotation %in% c("Early_EPI","Late_EPI")] <- "epi"
markers_expr$annot_lab[markers_expr$annotation %in% c("PriS","Mesoderm", "AdvMes")] <- "mes"

markers_expr$marker_type <- factor(
  markers_expr$marker_type,
  levels = sort(unique(markers_expr$marker_type))[c(1,3,2)]
)

# Aggregate
markers_expr_agg <- plyr::ddply(
  markers_expr[grep("^ZNF398_", markers_expr$marker_type, invert = T),],
  .(gene_name, annotation, marker_type),
  summarize,
  median = median(value),
  perc_expr = sum(value>0)/length(value)
)
markers_expr_agg$gene_name <- factor(
  markers_expr_agg$gene_name,
  levels = sort(unique(as.character(markers_expr_agg$gene_name)), decreasing = T)
)
# Visualization bubble ----
markers_expr_bubble_plt <- ggplot(
  markers_expr_agg[markers_expr_agg$perc_expr>0,],
  aes(x=annotation, y=gene_name)
) +
  geom_point(shape=21, aes(fill=median, size=perc_expr)) +
  facet_grid(marker_type~., scales = 'free', space = 'free') +
  theme_bw() +
  theme(
    legend.position  = "right",
    line             = element_line(size = 0.5),
    legend.title     = element_text(color="black",size=7),
    legend.text      = element_text(color="black",size=7),
    axis.title       = element_blank(),
    axis.text        = element_text(color="black",size=7),
    axis.text.x      = element_text(angle = 45, hjust = 1, vjust = 1, size=7),
    axis.line        = element_line(size=0.25),
    panel.grid.minor = element_blank(),
    strip.text.x       = element_text(size=7),
    strip.text.y       = element_text(size=7),
    strip.background = element_blank(),
    plot.title       = element_text(size = 8, hjust = 0.5, face = "plain"),
    plot.background  = element_rect(size = 0.25),
    panel.background = element_rect(size = 0.25),
    panel.border     = element_rect(size=0.25),
    axis.ticks       = element_line(size = 0.25),
    legend.key.size = unit(4,'mm'),
    strip.background.y = element_rect(color = "black", size = .8, fill = NA)
  ) +
  scale_size(range = c(1, 3.8)) +
  guides(
    fill = guide_colorbar(title='Median (log) Expr.', barwidth = .5, barheight = 2.8),
    size = guide_legend(title = '% Expr. Cells')
  ) +
  scale_fill_gradientn(
    colours =  RColorBrewer::brewer.pal(9, 'Oranges'),
  )

outfile <- paste0(path_results, "/markers_expr_bubble_plt.pdf")
message("-- saving to: ", outfile)
pdf(file = outfile, paper = "a4", width = unit(3.75,'cm'), height = unit(2.2,'cm'))
print(markers_expr_bubble_plt)
dev.off()

# Visualization boxplot ----
markers_expr_box_plt <- list()
for(i in unique(markers_expr$marker_type)) {
  markers_expr_box_plt[[i]] <- ggplot(
    markers_expr[markers_expr$marker_type==i,],
    aes(x=annotation, y=value)
  ) +
    geom_boxplot(outlier.shape = NA, fill="white", show.legend = F) +
    ggrastr::rasterise(geom_jitter(width = .2,height = 0,size=0.3, alpha=0.4, col="#606060"),dpi=300) +
    ylab("Expression") +
    stat_summary(col="red", fun = "median", size=0.2) +
    facet_grid(marker_type~gene_name, scales = "free") +
    theme_bw() +
    theme(
      legend.position  = "right",
      line             = element_line(size = 0.5),
      legend.title     = element_text(color="black",size=7),
      legend.text      = element_text(color="black",size=7),
      axis.title.x     = element_blank(),
      axis.title.y     = element_text(color="black",size=7),
      axis.text        = element_text(color="black",size=7),
      axis.text.x      = element_text(angle = 45, hjust = 1, vjust = 1, size=7),
      axis.line        = element_line(size=0.25),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      strip.text       = element_text(size=7),
      strip.background = element_blank(),
      plot.title       = element_text(size = 8, hjust = 0.5, face = "plain"),
      plot.background  = element_rect(size = 0.25),
      panel.background = element_rect(size = 0.25),
      panel.border     = element_rect(size=0.25),
      axis.ticks       = element_line(size = 0.25),
      legend.key.size = unit(4,'mm')
    )
  outfile <- paste0(path_results, "/markers_expr_box_plt.",i,".pdf")
  message("-- saving to: ", outfile)
  pdf(file = outfile, paper = "a4", width = unit(6.5,'cm'), height = unit(2.6,'cm'))
  print(markers_expr_box_plt[[i]])
  dev.off()
}

# Visualization boxplot - selection ----
gene_sel <- c("EPCAM", "PODXL", "CDH2")
palette_features <- c("prelineage"="#9ECAE1", "epi"="#3F007D", "mes"="#41AB5D")

markers_expr_sel_box_plt <- ggplot(
  markers_expr[markers_expr$gene_name %in% gene_sel,],
  aes(x=annotation, y=value, col=annot_lab)
) +
  geom_boxplot(outlier.shape = NA, fill="white", show.legend = F) +
  geom_line(data=markers_expr_agg[markers_expr_agg$gene_name %in% gene_sel, ],
            aes(x=annotation, y=median, group=gene_name),
            inherit.aes = F,
            col="#3c78c6ff",
            lwd=.65) +
  geom_jitter(width = .15,height = 0,size=0.1, alpha=0.4, show.legend = F) +
  ylab("Expression") +
  stat_summary(col="#D94801", fun = "median", size=0.15) +
  scale_color_manual(values = palette_features) +
  facet_wrap(marker_type~gene_name, ncol = 1, scales = 'free') +
  theme_bw() +
  theme(
    legend.position  = "right",
    line             = element_line(size = 0.5),
    legend.title     = element_text(color="black",size=7),
    legend.text      = element_text(color="black",size=7),
    axis.title.x     = element_blank(),
    axis.title.y     = element_text(color="black",size=7),
    axis.text        = element_text(color="black",size=7),
    axis.text.x      = element_text(angle = 45, hjust = 1, vjust = 1, size=7),
    axis.line        = element_line(size=0.25),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    strip.text       = element_text(size=7),
    strip.background = element_blank(),
    plot.title       = element_text(size = 8, hjust = 0.5, face = "plain"),
    plot.background  = element_rect(size = 0.25),
    panel.background = element_rect(size = 0.25),
    panel.border     = element_rect(size=0.25),
    axis.ticks       = element_line(size = 0.25),
    legend.key.size = unit(4,'mm')
  )
# markers_expr_sel_box_plt <- ggrastr::rasterise(markers_expr_sel_box_plt, dpi=600)
outfile <- paste0(path_results, "/markers_expr_sel_box_plt.pdf")
message("-- saving to: ", outfile)
pdf(file = outfile, paper = "a4", width = unit(2,'cm'), height = unit(4.6,'cm'))
print(markers_expr_sel_box_plt)
dev.off()
