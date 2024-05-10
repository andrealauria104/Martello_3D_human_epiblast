source("scripts/resources.R")

setwd("/home/epigen/Martello_3D_human_epiblast/Martello_Polo_TGFb/analysis_oscab")

path_results = "results/quality/correct-batches_experiment_1_2_untreated_0"
if(!dir.exists(path_results)) dir.create(path_results, recursive = T)

sce <- readRDS("results/correct-batches_experiment_1_2_untreated_0/cluster-cells/sce.rds")

umap_plt_l <- list()
pal <- list()
pal_qt <- RColorBrewer::brewer.pal(9,"Oranges")
# assigned rate ----
sce$assigned_rate_qt <- cut(
  sce$Assigned_rate_total,
  breaks = quantile(sce$Assigned_rate_total, seq(0,1,.2)),
  include.lowest = T)


pal[["assigned_rate_qt"]]<-pal_qt[1:length(levels(sce$assigned_rate_qt))]

umap_plt_l[["assigned_rate_qt"]] <- plotsceReducedDim(
  sce = sce,
  dimred = "UMAP",
  col_by = "assigned_rate_qt",
  pal = pal[["assigned_rate_qt"]],
  point_size = 0.4
)

# assigned count ----
sce$assigned_reads_qt <- cut(
  sce$assigned_reads,
  breaks = quantile(sce$assigned_reads, seq(0,1,.2)),
  include.lowest = T)

pal[["assigned_reads_qt"]]<-pal_qt[1:length(levels(sce$assigned_reads_qt))]

umap_plt_l[["assigned_reads_qt"]] <- plotsceReducedDim(
  sce = sce,
  dimred = "UMAP",
  col_by = "assigned_reads_qt",
  pal = pal[["assigned_reads_qt"]],
  point_size = 0.4
)

# mito ----
sce$perc_mitochondrial_qt <- cut(
  sce$perc_mitochondrial,
  breaks = quantile(sce$perc_mitochondrial, seq(0,1,.2)),
  include.lowest = T)

pal[["perc_mitochondrial_qt"]]<-pal_qt[1:length(levels(sce$perc_mitochondrial_qt))]

umap_plt_l[["perc_mitochondrial_qt"]] <- plotsceReducedDim(
  sce = sce,
  dimred = "UMAP",
  col_by = "perc_mitochondrial_qt",
  pal = pal[["perc_mitochondrial_qt"]],
  point_size = 0.4
)

# combine plots
umap_qt <- patchwork::wrap_plots(umap_plt_l, nrow = 1)

outfile <- paste0(path_results,"/umap_plt_l.pdf")
pdf(file = outfile, paper = "a4r",width = unit(10,'cm'), height = unit(6,'cm'))
print(umap_qt)
dev.off()
