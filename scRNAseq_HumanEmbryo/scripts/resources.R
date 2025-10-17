# Libraries ----
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(batchelor))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(ggsci))
suppressPackageStartupMessages(library(ggrastr))
suppressPackageStartupMessages(library(openxlsx))

# Themes ----
# ggplot themes
my_theme <- ggplot2::theme(
  legend.position  = "right",
  line             = element_line(size = 0.5),
  legend.title     = element_text(color="black",size=8),
  legend.text      = element_text(color="black",size=8),
  axis.title       = element_text(color="black",size=8),
  axis.text        = element_text(color="black",size=8),
  axis.line        = element_line(size=0.25),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  strip.text       = element_text(size=8),
  strip.background = element_blank(),
  plot.title       = element_text(size = 8, hjust = 0.5, face = "plain"),
  plot.background  = element_rect(size = 0.25),
  panel.background = element_rect(size = 0.25),
  panel.border     = element_rect(size=0.25),
  axis.ticks       = element_line(size = 0.25),
  legend.key.size = unit(4,'mm'),
  aspect.ratio = 1
)

# Functions ----
# Plots
plot_qc_th <- function(sce, x_var, discard_var) 
{
  qc_th_plot_total_count <- scater::plotColData(sce, x=x_var, y="sum", 
                                        colour_by=discard_var) + scale_y_log10() + 
    ggtitle("Total count") + theme_bw() + my_theme + theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1))
  qc_th_plot_detected <- scater::plotColData(sce, x=x_var, y="detected", 
                                     colour_by=discard_var) + scale_y_log10() + 
    ggtitle("Detected features") + theme_bw() + my_theme + theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1))
  qc_th_plot_mito <- scater::plotColData(sce, x=x_var, y="subsets_Mito_percent", 
                                 colour_by=discard_var) + ggtitle("Mito percent") + 
    theme_bw() + my_theme + theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1))
  qc_th_plot_mito <- scater::plotColData(sce, x=x_var, y="subsets_Mito_percent", 
                                 colour_by=discard_var) + ggtitle("Mito percent") + 
    theme_bw() + my_theme + theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1))
  qc_th_plot_top <- scater::plotColData(sce, x=x_var, y="percent.top_50", 
                                colour_by=discard_var) + ggtitle("Top 50%") + 
    theme_bw() + my_theme + theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1))
  qc_th_plot_assigned_tot <- scater::plotColData(sce, x=x_var, y="Assigned_rate_total", 
                                         colour_by=discard_var) + ggtitle("Assigned rate - total") + 
    theme_bw() + my_theme + theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1))
  qc_th_plot_assigned_uniqmap <- scater::plotColData(sce, x=x_var, y="Assigned_rate_uniqmap", 
                                             colour_by=discard_var) + ggtitle("Assigned rate - uniquely mapped") + 
    theme_bw() + my_theme + theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1))
  #qc_plot_ercc <- scater::plotColData(sce, x="block", y="altexps_ERCC_percent", 
  #             colour_by=discard_var) + ggtitle("ERCC percent")
  qc_th_plot_list <- list(qc_th_plot_total_count,qc_th_plot_detected,qc_th_plot_mito,qc_th_plot_top,qc_th_plot_assigned_tot,qc_th_plot_assigned_uniqmap)
  qc_th_plots <- ggpubr::ggarrange(plotlist = qc_th_plot_list, common.legend = T, legend = "bottom", ncol=2,nrow=3, align = "hv")
  
  return(qc_th_plots)
}
plotsceReducedDim <- function(sce, dimred,
                              num_dim=2, dim_1 = 1, dim_2 = 2,
                              col_by=NULL, shape_by=NULL, point_size=1,
                              pal = NULL, marker=NULL, facet_ncol = 3 ,
                              scale_col_limits = c(-1.5,1.5),
                              legend_position = "right") 
{
  dimred_toplot <- as.data.frame(SingleCellExperiment::reducedDim(sce, type=dimred)[,1:num_dim])
  if(dimred=="TSNE") colnames(dimred_toplot) <- gsub("^V","tSNE",colnames(dimred_toplot))
  if(dimred=="UMAP") colnames(dimred_toplot) <- gsub("^V","UMAP",colnames(dimred_toplot))
  x_var <- colnames(dimred_toplot)[dim_1]
  y_var <- colnames(dimred_toplot)[dim_2]
  if(dimred=="PCA") {
    percent.var <- attr(reducedDim(sce,type=dimred), "percentVar")
    x_lab <- paste0(x_var," (",round(percent.var[dim_1],2),"%)")
    y_lab <- paste0(x_var," (",round(percent.var[dim_2],2),"%)")
  } else {
    x_lab <- x_var
    y_lab <- y_var
  }
  dimred_toplot$Sample <- rownames(dimred_toplot)
  dimred_toplot <- merge(dimred_toplot, as.data.frame(colData(sce)), by="Sample")
  rownames(dimred_toplot) <- dimred_toplot$Sample
  if(is.null(pal)) {
    if(!is.null(col_by)) {
      pal <- ggsci::pal_d3()(length(unique(dimred_toplot[,col_by])))
    } else {
      pal <- "black"
    }
  }
  
  if(is.null(marker)) {
    p0 <- ggplot(dimred_toplot, aes_string(x=x_var, y=y_var, col=col_by, shape=shape_by)) +
      ggrastr::rasterize(geom_point(size=point_size),dpi=300) +
      scale_color_manual(values=pal)+
      xlab(x_lab) + ylab(y_lab) +
      theme_bw() + my_theme +
      theme(legend.position = legend_position) +
      guides(col = guide_legend(ncol=1, override.aes = list(size=3)))
    
  } else {
    mnames <- unique(marker$name)
    if(length(mnames)>1) {
      tmp <- list()
      for(i in mnames) {
        tmp[[i]] <- dimred_toplot
        tmp_mark <- subset(marker, name==i)
        idx <- match(rownames(tmp[[i]]), tmp_mark$cell)
        tmp[[i]]$marker <- i
        tmp[[i]]$expression <- tmp_mark$expression[idx]
      }
      dimred_toplot <- do.call(rbind, tmp)
      rm(tmp_mark)
    } else {
      idx <- match(rownames(dimred_toplot), marker$cell)
      dimred_toplot$marker <- marker$name[idx]
      dimred_toplot$expression <- marker$expression[idx]
    }
    p0 <- ggplot(dimred_toplot, aes_string(x=x_var, y=y_var, col="expression")) +
      ggrastr::rasterize(geom_point(size=point_size),dpi=300) + 
      facet_wrap(~marker, ncol = facet_ncol) +
      scale_colour_gradient2(low="#2166AC", mid="#F7F7F7", high="#B2182B"
                             , midpoint = 0
                             , breaks = seq(-2,2,1)
                             , limits = scale_col_limits
                             , oob=squish) +
      guides(col = guide_colourbar(barwidth = 0.8)) +
      xlab(x_lab) + ylab(y_lab) + theme_bw() + my_theme + 
      theme(legend.position = legend_position, panel.grid.major = element_blank(),aspect.ratio = 1)
  }
  return(p0)
}
plotsceExpressionHeatmap <- function(sce, assay.type
                                     , scale=T
                                     , ridx=NULL
                                     , cidx=NULL
                                     , myLegend=NULL
                                     , myPalette=NULL
                                     , myZscale=NULL
                                     , text_size = 8
                                     , annotDF = NULL
                                     , annotCol =  NULL) 
{
  m <- as.matrix(assay(sce,assay.type))
  
  if(!is.null(ridx))  m <- m[ridx,]
  if(!is.null(cidx))  m <- m[,cidx]
  
  if(scale) {
    m_scaled <- t(apply(m, 1, scale))
    colnames(m_scaled) <- colnames(m)
  } else {
    m_scaled <- m
  }
  
  if(is.null(myPalette)) {
    myPalette <- c("blue","black","red")
  }
  
  if(is.null(myZscale)) {
    myZscale <- c(-2, 0, 2)
  }
  ramp <- circlize::colorRamp2(myZscale, myPalette)
  
  if(is.null(myLegend)) myLegend <- gsub("_"," ",assay.type)
  
  if(scale) {
    myLegend_title <- paste0("Z-score (",myLegend,")")
  } else {
    myLegend_title <- myLegend
  }
  if (!is.null(annotDF)) {
    if (!is.null(annotCol)) {
      ha_column <- ComplexHeatmap::HeatmapAnnotation(df  = annotDF, 
                                                     col = annotCol, 
                                                     annotation_legend_param = list(title_gp  = gpar(fontsize=text_size),
                                                                                    labels_gp = gpar(fontsize=text_size)))
    } else {
      ha_column <- ComplexHeatmap::HeatmapAnnotation(df  = annotDF, 
                                                     annotation_legend_param = list(title_gp  = gpar(fontsize=text_size),
                                                                                    labels_gp = gpar(fontsize=text_size)))
    }
  } else {
    ha_column <- new("HeatmapAnnotation")
  }
  hm <- ComplexHeatmap::Heatmap(m_scaled
                                , cluster_rows = F
                                , cluster_columns = F
                                , show_row_names = F
                                , show_column_names = F
                                , col = ramp
                                , heatmap_legend_param = list(title = myLegend_title,
                                                              title_gp = gpar(fontsize=text_size),
                                                              title_position = "topcenter",
                                                              legend_width  = unit(3, "cm"),
                                                              legend_height = unit(0.5, "mm"),
                                                              values_gp     = gpar(fontsize=text_size),
                                                              labels_gp     = gpar(fontsize=text_size),
                                                              legend_direction = "horizontal"
                                )
                                , top_annotation = ha_column
                                , use_raster = T
                                , raster_device = "tiff"
                                , raster_quality = 10)
  ComplexHeatmap::draw(hm, heatmap_legend_side = "bottom")
  return(hm)
}
# Utils
quantize_tech_feature <- function(sce, tech_quant_feature) 
{
  tech_quant_values <- quantile(colData(sce)[[tech_quant_feature]])
  colData(sce)[[paste0("quant_",tech_quant_feature)]] <- cut(colData(sce)[[tech_quant_feature]]
                                                             , breaks=c(-Inf, tech_quant_values[2], tech_quant_values[3], tech_quant_values[4], Inf)
                                                             , labels=c("q1","q2","q3","q4"))
  return(sce)
}
# helper functions 
get_comma_arglist <- function(argument)
{
  if(!is.null(argument) && grepl("\\,",argument)) {
    argument <- unlist(strsplit(argument,"\\,"))
  }
  return(argument)
}
get_palette_features <- function(metadata_features)
{
  pals <- ls('package:ggsci', pattern = 'pal')[1:length(metadata_features)]
  palette_features <- lapply(pals, function(p) {
    rcb_pal <- sample(rownames(subset(RColorBrewer::brewer.pal.info,category=="div")),1)
    c(get(p)()(9),RColorBrewer::brewer.pal(9,rcb_pal))
  })
  names(palette_features) <- metadata_features
  return(palette_features)
}

saveXLSresEdgeR2 <- function(res, outfile
                             , name = "edgeR"
                             , force = TRUE
                             , keepNA = TRUE)
{
  message("[+] Saving results to file: ", outfile)
  
  if(file.exists(outfile) & force==F) {
    warning(" -- appending to existing file.")
    wb <- openxlsx::loadWorkbook(outfile)
  } else {
    wb <- openxlsx::createWorkbook()
  }
  
  if(is.list(res) & !is.data.frame(res)) {
    for(i in names(res)) {
      openxlsx::addWorksheet(wb = wb, sheetName = i)
      openxlsx::writeData(wb = wb
                          , sheet = i
                          , x = res[[i]]
                          , keepNA = keepNA
                          , colNames = TRUE
                          , rowNames = FALSE)
    }
  } else {
    openxlsx::addWorksheet(wb = wb, sheetName = name)
    openxlsx::writeData(wb = wb
                        , sheet = name
                        , x = res
                        , keepNA = keepNA
                        , colNames = TRUE
                        , rowNames = FALSE)
  }
  
  openxlsx::saveWorkbook(wb = wb
                         , file = outfile
                         , overwrite = TRUE)
}
