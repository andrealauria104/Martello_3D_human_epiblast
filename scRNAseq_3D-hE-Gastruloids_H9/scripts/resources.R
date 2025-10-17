# Libraries ----
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(batchelor))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(ggsci))
suppressPackageStartupMessages(library(ggrastr))
suppressPackageStartupMessages(library(plyr))
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

blank_theme <- theme_minimal()+
  theme(
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=8, face="plain",hjust = 0.5),
    line = element_line(size=0.5),
    axis.text      = element_blank(),
    legend.title   = element_text(color="black",size=8),
    legend.text    = element_text(color="black",size=8),
    strip.text = element_text(color="black",size=8),
    axis.line = element_blank(),
    axis.title = element_text(hjust = 0,size=8),
    aspect.ratio = 1,
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
                              num_dim=2, dim_1=1, dim_2=2, sample_ord=NULL,
                              col_by=NULL, shape_by=NULL, point_size=1,
                              point_alpha=1, pal=NULL, marker=NULL,
                              facet_ncol=3, guide_colourbar_title = "Expression",
                              scale_col_low="#2166AC", scale_col_mid="#F7F7F7",
                              scale_col_high="#B2182B", scale_col_midpoint=0,
                              scale_col_breaks=seq(-2,2,1),
                              scale_col_limits=c(-1.5,1.5),
                              legend_position="right") {
  
  dimred_toplot <- as.data.frame(
    SingleCellExperiment::reducedDim(sce, type=dimred)[,1:num_dim]
    )
  
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
  
  if(!is.null(sample_ord)) {
    dimred_toplot <- dimred_toplot[match(sample_ord,dimred_toplot$Sample),]
    message(" --- head: ")
    print(head(dimred_toplot))
    message(" --- tail: ")
    print(tail(dimred_toplot))
  }
  
  if(is.null(pal)) {
    if(!is.null(col_by)) {
      pal <- ggsci::pal_d3()(length(unique(dimred_toplot[,col_by])))
    } else {
      pal <- "black"
    }
  }
  
  if(is.null(marker)) {
    p0 <- ggplot(dimred_toplot, aes_string(x=x_var, y=y_var, col=col_by, shape=shape_by)) +
      ggrastr::rasterize(geom_point(size=point_size,alpha=point_alpha),dpi=600) +
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
      ggrastr::rasterize(geom_point(size=point_size,alpha=point_alpha),dpi=600) + 
      facet_wrap(~marker, ncol = facet_ncol) +
      scale_colour_gradient2(
        low=scale_col_low, mid=scale_col_mid, high=scale_col_high,
        midpoint = scale_col_midpoint,
        breaks = scale_col_breaks,
        limits = scale_col_limits,
        oob=scales::squish) +
      guides(
        col=guide_colourbar(
          title=guide_colourbar_title, barwidth=.6,
          barheight = 3,title.hjust = 0)
      ) +
      xlab(x_lab) + ylab(y_lab) + theme_bw() + my_theme
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

plotClusterMarkersHeatmap <- function(res
                                      , markers_ntop = 50
                                      , myPalette = NULL
                                      , myScale = NULL
                                      , gene_labels = NULL
                                      , myLegend_title = "logFC"
                                      , fdrTh = 0.05
                                      , text_size = 8) 
{
  if("Top"%in%colnames(res[[1]])) {
    tmp <- lapply(res, function(x) getMarkerEffects(x[x$Top<=markers_ntop,])) 
    column_title <- paste0("Cluster markers \n[Any pairwise comparison, Top <=",markers_ntop,"]")
  } else {
    tmp <- lapply(res, function(x) getMarkerEffects(x[x$FDR<=fdrTh,]))
    column_title <- paste0("Cluster markers [FDR <=",fdrTh,"]")
  }
  tmp <- tmp[unlist(lapply(tmp,nrow))!=0]
  tmp2 <- lapply(names(tmp), function(x) {
    y <- cbind(tmp[[x]],x=rep(NA,nrow(tmp[[x]])))
    colnames(y) <- gsub("x",x,colnames(y))
    return(as.data.frame(y))
  })
  names(tmp2) <- names(tmp)
  tmp3 <- do.call(rbind.data.frame, tmp2)
  nm_split <- gsub("\\..*","",rownames(tmp3))
  names(nm_split) <- rownames(tmp3)
  
  if(is.null(myPalette)) myPalette <- rev(RColorBrewer::brewer.pal(5,"RdBu"))
  if(is.null(myScale)) myScale <- c(-4, -1, 0, 1, 4)
  
  ramp <- circlize::colorRamp2(myScale, myPalette)
  
  if(is.null(gene_labels)) {
    set.seed(s33d)
    # gene_labels <- unlist(lapply(tmp, function(x) sample(rownames(x),5)))
    gene_labels <- unlist(lapply(tmp, function(x) rownames(x)[1:5]))
  } 
  hm_list <- list()
  for(i in names(tmp2)) {
    
    idx <- na.omit(match(gene_labels, rownames(tmp2[[i]])))
    if(length(idx)>0) {
      if(any(is.na(idx))) idx <- idx[!is.na(idx)]
      ha <- rowAnnotation(foo = anno_mark(at = idx, labels = rownames(tmp2[[i]])[idx],labels_gp=gpar(fontsize=6),padding = 1,labels_rot=-45))
    } else {
      ha <- NULL
    }
    hm_list[[i]] <- Heatmap(as.matrix(tmp2[[i]][,order(colnames(tmp2[[i]]))])
                            , show_row_dend = F
                            , show_column_names = T
                            , show_row_names = F
                            , cluster_rows = T
                            , cluster_columns = F
                            , na_col = "white"
                            , col = ramp
                            , border = "black"
                            , column_names_rot = 0
                            , row_title = i
                            , column_title = column_title
                            , row_title_rot = 0
                            , row_title_gp = gpar(fontsize=text_size)
                            , column_title_gp = gpar(fontsize=text_size)
                            , column_names_gp = gpar(fontsize=text_size)
                            , row_names_gp = gpar(fontsize=text_size)
                            , show_heatmap_legend = ifelse(match(i,names(tmp2))==1,T,F)
                            , right_annotation = ha
                            , heatmap_legend_param = list(title = myLegend_title, 
                                                          title_gp = gpar(fontsize = text_size), title_position = "topcenter", 
                                                          legend_width = unit(2, "cm"), legend_height = unit(0.2, 
                                                                                                             "mm"), values_gp = gpar(fontsize = text_size), 
                                                          labels_gp = gpar(fontsize = text_size), legend_direction = "horizontal")
                            , use_raster = T
                            , raster_device = "tiff"
                            , raster_quality = 10)
  }
  
  hm_list_v <- Reduce(f = "%v%", hm_list)
  draw(hm_list_v, heatmap_legend_side="bottom")
  return(hm_list_v)
}

plot_ggpie <- function(data, variable, title=NULL, pal = NULL, precedence = NULL, nudge_x_val = 0, label_size = 2.5)
{
  if(is.vector(data)) {
    mp <- cbind.data.frame("var1" = names(data), reshape2::melt(data))
    mp$var1 <- as.factor(gsub("_"," ", mp$var1))
    if(is.null(precedence)) precedence <- as.character(mp$var1)
    mp$var1 <- factor(mp$var1, levels = precedence)
  } else {
    mp <- data
    colnames(mp) <- c("var1","value")
    
  } 
  
  n <- length(unique(mp$var1))
  mp$ypos <- cumsum(mp$value)[n] - cumsum(mp$value) + mp$value/2
  mp$perc <- round(mp$value/sum(mp$value), 4)
  
  if(is.null(pal)) pal <- RColorBrewer::brewer.pal(n,"Set2")
  
  pie <- ggplot(mp, aes(x = "", y = value, fill = var1)) + 
    geom_bar(width = 1, stat = "identity", col = "black", 
             lwd = 0.25) + 
    coord_polar("y", start = 0) + scale_fill_manual(values = pal) + ggtitle(title) +
    blank_theme + theme(axis.text.x = element_blank()
                        , plot.title = element_text(size = 8, face = "bold", hjust = 0.5, vjust = 0.7), text = element_text(size = 8)
                        , legend.key.size = unit(4,'mm')) + guides(fill=guide_legend(title=variable)) +
    geom_label(aes(label = scales::percent(perc), fill = var1, y = ypos), size = label_size, show.legend = F, nudge_x = nudge_x_val)
  return(pie)
}

prepareDimRedPltGenes <- function(sce, genes,
                                   cumulative_expression = F,
                                   assay.var = "logcounts",
                                   scale = T) {
  
  if(!cumulative_expression || length(genes)==1) {
    idx <- intersect(genes, rownames(sce))
    if(scale) {
      mat.genes <- t(scale(t(assay(sce, assay.var)[idx,, drop=F])))
    } else {
      mat.genes <- assay(sce, assay.var)[idx,, drop=F]
    }
    
    dimredplt.genes <- reshape2::melt(mat.genes,
                                     varnames = c("name","cell"),
                                     value.name = "expression")
  } else {
    if(!is.list(genes)) {
      message("[!] Converting genes to gene set list.")
      genes <- list(genes)
    }
    geneset.length <- unlist(lapply(genes,length))
    if(any(geneset.length<=1)) {
      message(
        paste0(
          "[!] Gene sets: ",
          paste0(which(geneset.length<=1),collapse = ","),
          " contain only 1 gene!"
          )
        )
    }
    geneset <- lapply(genes, function(x) {
      idx <- intersect(x, rownames(sce))
      colSums(assay(sce, assay.var)[idx,, drop=F])
    })
    geneset <- do.call(rbind, geneset)
    
    if(scale) {
      mat.geneset <- t(scale(t(geneset)))
    } else {
      mat.geneset <- geneset
    }
    
    dimredplt.genes <- reshape2::melt(
      mat.geneset,
      varnames = c("name","cell"),
      value.name = "expression"
      )
  }
  
  return(dimredplt.genes)
}
