setwd("/home/epigen/Martello_3D_human_epiblast/Martello_Polo_TGFb/analysis_oscab")
source("scripts/resources.R")
path_results <- "results/correct-batches_experiment_1_2_untreated_human_embryo/visualize-results"
sce <- readRDS("results/correct-batches_experiment_1_2_untreated_human_embryo/cluster-cells/sce.rds")

markers <- c(
  "POU5F1","ZNF398",
  "PARD3","PODXL",
  "GATA6","MIXL1","EOMES"
  )

markers.toplot <- prepareDimRedPltGenes(
  sce,
  genes = markers,
  assay.var = "logcounts",
  scale = T
  )

plotsceReducedDim(
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
) + blank_theme




tsne_markers_plt_l <- lapply(
  markers,
  function(genes) {
    plt <- plotTSNE(sce, col=genes,point_size=0.4)+
      theme_bw() + 
      my_theme +
      xlab("tSNE 1") +
      ylab("tSNE 2") +
      scale_color_distiller(
        'pr',
        palette='RdBu',
        breaks = c(1,2,4,6),
        limits = c(0,max(logcounts(sce)[genes,]))
        ) +
      guides(
        col=guide_colourbar(
          title=NULL, barwidth=.5,
          barheight = 2.2, label.theme = element_text(size=6))
        ) +
      ggtitle(genes)
    plt <- ggrastr::rasterise(plt, dpi=300)
    return(plt)
  }
)

tsne_markers_plt <- patchwork::wrap_plots(
  tsne_markers_plt_l,
  ncol = 2
)

outfile <- paste0(path_results,"/tsne_markers.pdf")
pdf(file = outfile, paper = "a4",width = unit(8,'cm'), height = unit(5,'cm'))
print(tsne_markers_plt)
dev.off()
