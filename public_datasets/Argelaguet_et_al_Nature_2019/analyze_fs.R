setwd("/home/epigen/Martello_3D_human_epiblast/public_datasets/Argelaguet_et_al_Nature_2019")
library(Seurat)
library(reshape2)
library(ggplot2)
library(plyr)
library(scater)
library(scran)

path_results <- "finalized_scran"

if(!dir.exists(path_results)) dir.create(path_results, recursive = T)
# Read seurat object ----
seurat_obj <- readRDS("GSE133725.seurat_obj.rds")

sce <- as.SingleCellExperiment(seurat_obj)
sce <- scran::computeSumFactors(sce)
sce <- scater::logNormCounts(sce)

# Define markers panel ----
gene_sets_l <- list(
  "Epithelial" = c("CDH1", "ESRP1", "EPCAM"),
  "Polarity"=  c("CLDN6","CLDN7", "OCLN", "PODXL"),
  "Mesenchymal" = c("CDH2", "ZEB2", "VIM")
)
gene_sets_l <- reshape2::melt(gene_sets_l)
colnames(gene_sets_l) <- c("human_gene_name", "marker_type")

# mouse homologs ----
gene_info <- read.delim("/home/reference_data/bioinfotree/task/gencode/dataset/mmusculus/M23/primary_assembly.annotation.ensg2gene_symbol2biotype.map.header_added")
gene_info$ens_id <- gsub("\\.\\d+$","", gene_info$gene_id)

gene_homologs <- read.delim("/home/reference_data/bioinfotree/task/gencode/dataset/mmusculus/M23/ensembl_mart.hsapiens.homolog.one2one_through_gene_name.header_added.gz")
gene_sets_l <- merge(
  gene_sets_l,
  gene_homologs[,c("ensembl_gene_id", "external_gene_name",
                   "homolog_ensembl_gene", "homolog_associated_gene_name")],
  by.x = "human_gene_name",
  by.y = "homolog_associated_gene_name",
  all.x = T
)
# colnames(gene_sets_l)[3] <- "gene_name"
gene_sets_l <- gene_sets_l[order(gene_sets_l$marker_type, gene_sets_l$human_gene_name),]

# Markers expr ----
# markers_expr <- as.matrix(Seurat::GetAssayData(seurat_obj, assay = "RNA", slot = "data"))
markers_expr <- as.matrix(assay(sce, 'logcounts'))
markers_expr <- markers_expr[which(rownames(markers_expr) %in% unique(gene_sets_l$ensembl_gene_id)),]
markers_expr <- reshape2::melt(markers_expr, varnames = c("ensembl_gene_id", "cell_id"))
metadata <- seurat_obj[[]]
metadata <- metadata[, c("orig.ident", "plate", "stage", "lineage10x", "lineage10x_2")]
metadata$cell_id <- rownames(metadata)
metadata <- metadata[-which(is.na(metadata$lineage10x_2)),]

markers_expr <- merge(markers_expr, metadata, by = "cell_id")
markers_expr <- merge(markers_expr, gene_sets_l, by = "ensembl_gene_id")
# genes <- unique(unlist(gene_sets_l))
# expr_meta <- FetchData(seurat_obj, vars = c(genes, "cell_id","rename_EML", "sub_rename_EML"))

# Refine annotation ----
markers_expr$annotation <- paste0(
  markers_expr$stage,
  "_",
  markers_expr$lineage10x_2
  )

stage_ord <- c(
  "E4.5_Epiblast",
  "E5.5_Epiblast",
  "E6.5_Epiblast", "E6.5_Primitive_Streak", "E6.5_Mesoderm"
)
markers_expr <- markers_expr[markers_expr$annotation %in% stage_ord, ]
markers_expr$annotation <- factor(
  markers_expr$annotation,
  levels = stage_ord
)

markers_expr$annot_lab <- "mes"
markers_expr$annot_lab[grepl("Epiblast", as.character(markers_expr$annotation))] <- "epi"

markers_expr$marker_type <- factor(
  markers_expr$marker_type,
  levels = sort(unique(markers_expr$marker_type))[c(1,3,2)]
)

levels(markers_expr$annotation) <- gsub("_Epiblast","_EPI",levels(markers_expr$annotation))
levels(markers_expr$annotation) <- gsub("_Primitive_Streak","_PriS",levels(markers_expr$annotation))
levels(markers_expr$annotation) <- gsub("_Mesoderm","_Mesod.",levels(markers_expr$annotation))

# Aggregate
markers_expr_agg <- plyr::ddply(
  markers_expr,
  .(external_gene_name, annotation, marker_type),
  summarize,
  median = median(value),
  perc_expr = sum(value>0)/length(value)
)

markers_expr_agg$external_gene_name <- factor(
  markers_expr_agg$external_gene_name,
  levels = sort(unique(as.character(markers_expr_agg$external_gene_name)), decreasing = T)
)

# Visualization bubble ----
markers_expr_bubble_plt <- ggplot(
  markers_expr_agg[markers_expr_agg$perc_expr>0,],
  aes(x=annotation, y=external_gene_name)
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
pdf(file = outfile, paper = "a4", width = unit(2.75,'cm'), height = unit(2.3,'cm'))
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
    # stat_summary(col="red", fun = "median", size=0.2) +
    facet_grid(marker_type~external_gene_name, scales = "free") +
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
gene_sel <- c("Epcam", "Podxl", "Cdh2")
palette_features <- c("prelineage"="#9ECAE1", "epi"="#3F007D", "mes"="#41AB5D")

markers_expr_sel_box_plt <- ggplot(
  markers_expr[markers_expr$external_gene_name %in% gene_sel,],
  aes(x=annotation, y=value, col=annot_lab)
) +
  geom_boxplot(outlier.shape = NA, fill="white", show.legend = F) +
  geom_line(data=markers_expr_agg[markers_expr_agg$external_gene_name %in% gene_sel, ],
            aes(x=annotation, y=median, group=external_gene_name),
            inherit.aes = F,
            col="#3c78c6ff",
            lwd=.65) +
  geom_jitter(width = .15,height = 0,size=0.1, alpha=0.4, show.legend = F) +
  ylab("Expression") +
  stat_summary(col="#D94801", fun = "median", size=0.15) +
  scale_color_manual(values = palette_features) +
  facet_wrap(marker_type~external_gene_name, ncol = 1, scales = "free") +
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

outfile <- paste0(path_results, "/markers_expr_sel_box_plt.pdf")
message("-- saving to: ", outfile)
pdf(file = outfile, paper = "a4", width = unit(1.3,'cm'), height = unit(4.8,'cm'))
print(markers_expr_sel_box_plt)
dev.off()
