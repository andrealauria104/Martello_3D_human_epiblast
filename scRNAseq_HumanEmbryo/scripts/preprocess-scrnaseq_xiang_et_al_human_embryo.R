# Analyze scRNAseq - Bioconductor workflow
# Pre-processing
# 0. Resources ----
library(RNAseqRtools)
library(scRNAseqRtools)
library(SingleCellExperiment)
library(scran)
library(scater)
library(pheatmap)
library(limma)
library(batchelor)

# paths ---
path_biotypes_stats <- "../Rmarkdown/Results/Xiang-HumanEmbryo/biotypes_stats.txt"
path_mapping_stats <- "../Rmarkdown/Results/Xiang-HumanEmbryo/mapping_stats.txt"
path_counts <-  "data/xiang_et_al_human_embryo/raw_counts.csv.gz"
path_ginfo <- "data/gene_info.gencode.hsapiens.32.basic.annotation.txt.gz"
# output dirs ---
FIGDIR  <- "results/xiang_et_al_human_embryo/preprocess"
RESDIR  <- "results/xiang_et_al_human_embryo/preprocess"

if(!dir.exists(FIGDIR)) dir.create(FIGDIR, recursive = T)
if(!dir.exists(RESDIR)) dir.create(RESDIR, recursive = T)

s33d <- 1991
# 0.1 Functions ----
# Plots
plot_qc_th <- function(sce, x_var, discard_var) 
{
  qc_th_plot_total_count <- plotColData(sce, x=x_var, y="sum", 
                                        colour_by=discard_var) + scale_y_log10() + 
    ggtitle("Total count") + theme_bw() + my_theme + theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1))
  qc_th_plot_detected <- plotColData(sce, x=x_var, y="detected", 
                                     colour_by=discard_var) + scale_y_log10() + 
    ggtitle("Detected features") + theme_bw() + my_theme + theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1))
  qc_th_plot_mito <- plotColData(sce, x=x_var, y="subsets_Mito_percent", 
                                 colour_by=discard_var) + ggtitle("Mito percent") + 
    theme_bw() + my_theme + theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1))
  qc_th_plot_mito <- plotColData(sce, x=x_var, y="subsets_Mito_percent", 
                                 colour_by=discard_var) + ggtitle("Mito percent") + 
    theme_bw() + my_theme + theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1))
  qc_th_plot_top <- plotColData(sce, x=x_var, y="percent.top_50", 
                                colour_by=discard_var) + ggtitle("Top 50%") + 
    theme_bw() + my_theme + theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1))
  qc_th_plot_assigned_tot <- plotColData(sce, x=x_var, y="Assigned_rate_total", 
                                         colour_by=discard_var) + ggtitle("Assigned rate - total") + 
    theme_bw() + my_theme + theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1))
  qc_th_plot_assigned_uniqmap <- plotColData(sce, x=x_var, y="Assigned_rate_uniqmap", 
                                             colour_by=discard_var) + ggtitle("Assigned rate - uniquely mapped") + 
    theme_bw() + my_theme + theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1))
  #qc_plot_ercc <- plotColData(sce, x="block", y="altexps_ERCC_percent", 
  #             colour_by=discard_var) + ggtitle("ERCC percent")
  qc_th_plot_list <- list(qc_th_plot_total_count,qc_th_plot_detected,qc_th_plot_mito,qc_th_plot_top,qc_th_plot_assigned_tot,qc_th_plot_assigned_uniqmap)
  qc_th_plots <- ggpubr::ggarrange(plotlist = qc_th_plot_list, common.legend = T, legend = "bottom", ncol=2,nrow=3, align = "hv")
  
  return(qc_th_plots)
}
plotsceReducedDim <- function(sce, dimred
                              , num_dim=2
                              , dim_1 = 1
                              , dim_2 = 2
                              , col_by=NULL
                              , shape_by=NULL
                              , point_size=1
                              , pal = NULL
                              , marker=NULL
                              , facet_ncol = 3 
                              , scale_col_limits = c(-1.5,1.5)
                              , legend_position = "right") 
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
      geom_point(size=point_size) + scale_color_manual(values=pal)+
      xlab(x_lab) + ylab(y_lab) + guides(col = guide_legend(ncol=1,override.aes=list(size=2.5))) + 
      theme_bw() + my_theme +
      theme(legend.position = legend_position, panel.grid.major = element_blank(),aspect.ratio = 1)
    
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
      geom_point(size=point_size) + facet_wrap(~marker, ncol = facet_ncol) +
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
  panel.grid       = element_line(size = 0.25),
  axis.ticks       = element_line(size = 0.25),
  legend.key.size = unit(4,'mm'),
  aspect.ratio = 1
)
# 0.2 Global params ----
metadata_features <- c("instrument_model","age.ch1","cell.type.ch1")
organism <- "hsapiens" # hsapiens | mmusculus

assay.type  <- "logcounts"
norm_method <- "libsize" # pool | pool_cluster | libsize
norm_filter <- NULL # rRNA_mito | NULL
npc_sel     <- 10L

batch_var  <- NULL
design_var <- "age.ch1"
batch_corr_method <- "fastMNN" # limma | regressBatches | rescaleBatches | fastMNN
batch_corr_use_design <- FALSE
shape_by <- NULL

explained_var_features <- c(metadata_features
                            , "Assigned_rate_total"
                            , "assigned_reads"
                            , "perc_mitochondrial"
                            , "ngenes")

additional_filters <- FALSE
save_additional_files <- FALSE
# 1. Load and prepare data -----
# Load data ---
mapping_stats <- read.delim(path_mapping_stats)
counts        <- read.delim(path_counts, stringsAsFactors = F, header = T,sep = ",", row.names = 1)
colnames(counts) <- gsub("^X","",colnames(counts))
colnames(counts) <- gsub("_S\\d+.*","",colnames(counts))
biotypes_stats  <- read.delim(path_biotypes_stats, stringsAsFactors = F, header = T)
biotypes_stats <- biotypes_stats[,c("Sample",colnames(biotypes_stats)[!colnames(biotypes_stats)%in%colnames(mapping_stats)])]
metadata <- merge(mapping_stats,biotypes_stats,by="Sample")
rownames(metadata) <- metadata$Sample
metadata <- metadata[colnames(counts),]

gene_info  <- read.delim(path_ginfo, stringsAsFactors = F, header = T)
rownames(gene_info) <- gene_info$gene_name
gene_info <- gene_info[rownames(counts),]

# Create SCE object ---
sce <- SingleCellExperiment::SingleCellExperiment(
  assays = list(counts = as.matrix(counts)),
  colData = metadata,
  rowData=gene_info
  )

# Filter empty|minibulk|zero-read cells
# sce_full <- sce
if(sum(sce$assigned_reads==0)>0) sce <- sce[,-which(sce$assigned_reads==0)]
if(any(grepl("MBULK|EMPTY|MINIBULK",colnames(sce)))) sce <- sce[,-grep("MBULK|EMPTY|MINIBULK",colnames(sce))]

# Mark genes to exclude for normalization
if(!is.null(norm_filter) && norm_filter=="rRNA_mito") {
  rowData(sce)$use_for_norm <- !(grepl("rRNA_",rowData(sce)$gene_type) | grepl("^MT-",rowData(sce)$gene_name))  
}

# Palettes for metadata ---
metadata_features <- get_comma_arglist(metadata_features)
palette_features  <- get_palette_features(metadata_features)
# 2. QC for marking and/or filtering low quality cells -----
# 2.1 Level 1: Fixed thresholds ----
is.mito <- grepl("^mt-|^MT-",rownames(sce))
names(is.mito) <- rownames(sce)
qc_df <- scater::perCellQCMetrics(sce, subsets=list(Mito=is.mito),percent.top=50)
sce <- scater::addPerCellQC(sce, subsets=list(Mito=is.mito),percent.top=50)

qc.assigned <- sce$Assigned_rate_total < 20
qc.lib <- qc_df$sum < 1e5
qc.nexprs <- qc_df$detected < 2e3
# qc.spike <- qc_df$altexps_ERCC_percent > 10
qc.mito <- qc_df$subsets_Mito_percent > 25 #10
qc.mito[which(is.na(qc.mito))] <- FALSE
qc.top <- qc_df$percent.top_50 > 90 #10
qc.top[which(is.na(qc.top))] <- FALSE

# default: qc.lib | qc.nexprs | qc.spike | qc.mito | qc.top | qc.assigned
discard_fixed_th <- qc.lib | qc.nexprs | qc.mito | qc.top | qc.assigned
# custom
# discard_fixed_th <- qc.lib | qc.nexprs

fixed_th_df <- data.frame(discard_fixed_th_LibSize=qc.lib
                          , discard_fixed_th_NExprs=qc.nexprs
                          , discard_fixed_th_Mito_percent=qc.mito
                          , discard_fixed_th_percent_top_50=qc.top
                          , discard_fixed_th_assigned_rate=qc.assigned
                          , discard_fixed_th=discard_fixed_th)
rownames(fixed_th_df) <- rownames(qc_df)

discard_fixed_th_df_stats <- data.frame(LibSize=sum(qc.lib)
                                        , NExprs=sum(qc.nexprs)
                                        # ,SpikeProp=sum(qc.spike)
                                        , MitoProp=sum(qc.mito)
                                        , TopProp=sum(qc.top)
                                        , Total_cells=length(discard_fixed_th)
                                        , Total_discard=sum(discard_fixed_th)
                                        , Total_keep=sum(!discard_fixed_th))
colData(sce) <- cbind.DataFrame(colData(sce),fixed_th_df)

qc_fixed_th_plots <- list()
for(feature in metadata_features) {
  qc_fixed_th_plots[[feature]] <- plot_qc_th(sce
                                             , x_var = feature
                                             ,discard_var = "discard_fixed_th")
  outfile <- paste0(FIGDIR,"/QC_plots_fixed_th_",feature,".pdf")
  message(" -- saving figure: ", outfile)
  pdf(file = outfile, paper = "a4",w=unit(6,'cm'),h=unit(10,"cm"),useDingbats = F)
  print(qc_fixed_th_plots[[feature]])
  dev.off()
}

# 2.2 Level 2: Adaptive thresholds ----
# Identify cells that are outliers for the various QC metrics, 
# based on the median absolute deviation (MAD)
# Use as primary or additional filter (if shallow filter on fixed thresholds) to discard/mark cells
# from the median value of each metric across all cells
adaptive_th_df <- scater::quickPerCellQC(qc_df, percent_subsets = c("subsets_Mito_percent","percent.top_50"))
colnames(adaptive_th_df) <- gsub("discard","discard_adaptive_th",colnames(adaptive_th_df))
discard_adaptive_th_df_stats <- colSums(as.matrix(adaptive_th_df))
colData(sce) <- cbind.DataFrame(colData(sce),adaptive_th_df)

# qc_adaptive_th_plots_time_treatment <- plot_qc_th(sce, x_var = "cell.type.ch1",discard_var = "discard_adaptive_th")

# qc_adaptive_th_plots <- list()
# for(feature in metadata_features) {
#   qc_adaptive_th_plots[[feature]] <- plot_qc_th(sce
#                                              , x_var = feature
#                                              , discard_var = "discard_adaptive_th")  
#   outfile <- paste0(FIGDIR,"/QC_plots_MAD_adaptive_th_",qc_x_var,".pdf")
#   message(" -- saving figure: ", outfile)
#   pdf(file = outfile, paper = "a4",w=unit(6,'cm'),h=unit(10,"cm"),useDingbats = F)
#   print(qc_adaptive_th_plots[[feature]])
#   dev.off()
# }

# 2.3 Discard -----
# discard_var <- "discard_fixed_th"
# sce$discard <- sce[[discard_var]]
# unfiltered  <- sce # save unfiltered copy
# table(unfiltered$discard)
# sce <- sce[,!sce$discard]

# The hard way
# table(!(unfiltered$assigned_reads<1e5 | unfiltered$ngenes < 2000 | unfiltered$perc_mitochondrial > 25))
# 3. Normalization ----
print(norm_method)
if(norm_method=="pool_cluster") {
  norm_clusters <- scran::quickCluster(sce) # for larger data sets
  if(is.null(norm_filter)) {
    sce <- scran::computeSumFactors(sce, clusters=norm_clusters)
  } else {
    sce <- scran::computeSumFactors(sce, clusters=norm_clusters, subset.row=rowData(sce)$use_for_norm)
  }
} else if(norm_method=="pool")  {
  if(is.null(norm_filter)) {
    sce <- scran::computeSumFactors(sce)
  } else {
    sce <- scran::computeSumFactors(sce, subset.row=rowData(sce)$use_for_norm)
  }
} else if(norm_method=="libsize") {
  message(" -- using librarySizeFactors only")
}
sce <- scater::logNormCounts(sce)

size_factors_df <- data.frame("librarySizeFactors"=librarySizeFactors(sce)
                              ,"DeconvolutionFactors"=sizeFactors(sce))
p_size_factor <- list()
for(feature in metadata_features) {
  size_factors_df[[feature]] <- colData(sce)[,feature]
  p_size_factor[[feature]] <- ggplot(size_factors_df, aes_string(x = "librarySizeFactors"
                                                                 , y = "DeconvolutionFactors"
                                                                 , col = feature)) +
    geom_point(size=1) + theme_bw() + my_theme + 
    scale_x_log10() + scale_y_log10() + 
    scale_color_manual(values = palette_features[[feature]])
  
  outfile <- paste0(FIGDIR,"/size_factors_col_feature_",feature,".pdf")
  pdf(file = outfile, paper = "a4",w=unit(3,'cm'),h=unit(3,"cm"),useDingbats = F)
  print(p_size_factor[[feature]])
  dev.off()
}

# 4. Variance Modeling for feature selection ----
if(is.null(batch_var)) {
  if(is.null(norm_filter)) {
    dec <- scran::modelGeneVar(sce)  
  } else {
    dec <- scran::modelGeneVar(sce, subset.row=rowData(sce)$use_for_norm)
  }
  
  pdec <- ggplot(as.data.frame(dec),aes(x=mean,y=total)) + geom_point(size=0.5) +
    geom_line(aes(x=mean,y=tech),col="red") +
    theme_bw() + my_theme + xlab("Mean of log-expression")+
    ylab("Variance of log-expression")
  
  outfile <- paste0(FIGDIR,"/variance_modeling.pdf")
  pdf(file = outfile, paper = "a4",w=unit(3,'cm'),h=unit(3,"cm"),useDingbats = F)
  print(pdec)
  dev.off()
  
} else {
  if(is.null(norm_filter)) {
    dec <- scran::modelGeneVar(sce, block=sce[[batch_var]],subset.row=rowData(sce)$use_for_norm) # block for batch  
  } else {
    dec <- scran::modelGeneVar(sce, block=sce[[batch_var]]) # block for batch
  }
  pdec <- list()
  blocked.stats <- dec$per.block
  for (i in colnames(blocked.stats)) {
    current <- blocked.stats[[i]]
    curfit <- metadata(current)
    current$trend <- curfit$trend(curfit$mean)
    pdec[[i]] <- ggplot(as.data.frame(current),aes(x=mean,y=total)) + 
      ggrastr::rasterise(geom_point(size=0.5),dpi=300) +
      geom_line(aes(x=mean,y=tech),col="red") +
      theme_bw() + my_theme + xlab("Mean of log-expression")+
      ylab("Variance of log-expression") + ggtitle(paste0(batch_var," ",i))
  }
  if(length(pdec)>6) {
    arrange_ncol <- 6
    arrange_nrow <- ceiling(length(pdec)/6)
  } else {
    arrange_ncol <- length(pdec)
    arrange_nrow <- 1
  }
  
  h_value <- 3*arrange_nrow
  
  outfile <- paste0(FIGDIR,"/variance_modeling.pdf")
  pdf(file = outfile, paper = "a4r",w=unit(12,'cm'),h=unit(h_value,"cm"),useDingbats = F)
  print(ggpubr::ggarrange(plotlist = pdec, nrow = arrange_nrow,ncol=arrange_ncol))
  dev.off()
}
outfile  <- paste0(RESDIR,"/dec.rds")
message(" -- saving: ",outfile)
saveRDS(dec, file = outfile)

# Select top variable genes ---

# dec_id <- intersect(rownames(dec)[dec$mean>=1], names(which(rowSums(counts(sce)>=25)>=1)))
# Get the top 10% of genes.
# top.hvgs <- getTopHVGs(dec[dec_id,], prop=0.1)
top.hvgs <- getTopHVGs(dec, prop=0.1)

# Get the top 2000 genes.
# top.hvgs2 <- getTopHVGs(dec[dec_id,], n=2000)
top.hvgs2 <-  getTopHVGs(dec, n=2000)

# Get all genes with positive biological components.
# top.hvgs3 <- getTopHVGs(dec[dec_id,], var.threshold=0)
top.hvgs3 <- getTopHVGs(dec, var.threshold=0)

# Get all genes with FDR below 5%.
# top.hvgs4 <- getTopHVGs(dec[dec_id,], fdr.threshold=0.05)
top.hvgs4 <- getTopHVGs(dec, fdr.threshold=0.05)

# vargenes <- top.hvgs2
vargenes <- top.hvgs2
varmodel <- dec[vargenes,]

outfile  <- paste0(RESDIR,"/variance_modeling.xlsx")
varmodel <- cbind("gene_name"=rownames(varmodel),varmodel)
RNAseqRtools::saveXLSresEdgeR2(varmodel, outfile = outfile)

dec_top <- as.data.frame(dec)
dec_top$top <- FALSE
dec_top[vargenes,"top"] <- TRUE
pdec_top <- ggplot(dec_top,aes(x=mean,y=total,col=top)) + 
  ggrastr::rasterise(geom_point(size=0.5),dpi=300) +
  geom_line(aes(x=mean,y=tech),col="red") +
  theme_bw() + my_theme + xlab("Mean of log-expression")+
  ylab("Variance of log-expression") + scale_color_manual(values=c("black","#2171B5"))

outfile <- paste0(FIGDIR,"/variance_modeling_top.pdf")
pdf(file = outfile, paper = "a4",w=unit(3,'cm'),h=unit(3,"cm"),useDingbats = F)
print(pdec_top)
dev.off()

rowData(sce)$vargenes <- FALSE
rowData(sce)$vargenes[match(vargenes,rowData(sce)$gene_name)] <- TRUE

# rowData(sce)$vargenes.seurat <- FALSE
# seurat_obj <- readRDS("../analysis_seurat/results/seurat-standard_human_embryo/seurat_obj.rds")
# rowData(sce)$vargenes.seurat[match(vargenes,Seurat::VariableFeatures(seurat_obj))] <- TRUE

# 5. Batch correction ----
if(!is.null(batch_var)) {
  # Correction strategy for uninteresting effects 
  # For more complex dataset integration, use scrnaseq-correct-batches.r script
  if(batch_corr_method=="limma") {
    if(batch_corr_use_design) {
      assay(sce, "corrected") <- limma::removeBatchEffect(logcounts(sce)
                                                          , design = model.matrix(as.formula(paste0("~",design_var)), data = colData(sce))
                                                          , batch = sce[[batch_var]])
    } else {
      assay(sce, "corrected") <- limma::removeBatchEffect(logcounts(sce)
                                                          , batch = sce[[batch_var]])
    }
    assay.type <- "corrected"
  } else if(batch_corr_method=="regressBatches") {
    if(batch_corr_use_design) {
      sce <- correctExperiments(sce, 
                                PARAM = RegressParam(
                                  design = model.matrix(~sce[[design_var]] + sce[[batch_var]]),
                                  keep = 1:length(unique(sce[[design_var]])))
      )
    }
    else {
      sce <- correctExperiments(sce, batch = sce[[batch_var]],
                                PARAM = RegressParam()
      )
    }
    assay.type <- "corrected"
  } else if(batch_corr_method=="rescaleBatches") {
    sce <- correctExperiments(sce, batch = sce[[batch_var]],
                              PARAM = RescaleParam()
                              )
    assay.type <- "corrected"
  } else if(batch_corr_method=="fastMNN") {
    f.out <- fastMNN(sce, batch=sce[[batch_var]], subset.row=vargenes)
    reducedDim(sce, "corrected") <- reducedDim(f.out, "corrected")
  }
}
# 6. Dimensionality reduction ----
# PCA ---
sce <- scater::runPCA(
  sce,
  ncomponents=100,
  subset_row=rowData(sce)$vargenes,
  exprs_values=assay.type,
  BSPARAM=BiocSingular::ExactParam()
  )

# Number of PCs - Elbow method ---
percent.var <- attr(reducedDim(sce), "percentVar")
chosen.elbow <- PCAtools::findElbowPoint(percent.var)
reducedDim(sce, "PCA.elbow") <- reducedDim(sce)[,1:chosen.elbow]
reducedDimNames(sce)

# Number of PCs - Population structure: clustered PCs ---
pcs <- reducedDim(sce, "PCA")
choices <- getClusteredPCs(pcs)
val <- metadata(choices)$chosen
reducedDim(sce, "PCA.clust") <- pcs[,1:val]

ppca <- list()
ppca_multi <- list()
for(feature in metadata_features) {
  ppca[[feature]] <- plotsceReducedDim(
    sce = sce,
    dimred = "PCA",
    col_by = feature,
    shape_by = shape_by,
    point_size = 0.5,
    pal = palette_features[[feature]]
    )
  
  outfile <- paste0(FIGDIR,"/ppca_",feature,".pdf")
  pdf(file = outfile, paper = "a4r",w=unit(3.5,'cm'),h=unit(3.5,"cm"))
  print(ppca[[feature]])
  dev.off()
  
  ppca_multi[[feature]] <- scater::plotPCA(
    sce,
    ncomponents = 4,
    colour_by = feature
    ) + 
    theme_bw() + 
    my_theme + 
    scale_color_manual(values=palette_features[[feature]])
  
  outfile <- paste0(FIGDIR,"/ppca_multi_",feature,".pdf")
  pdf(file = outfile, paper = "a4r",w=unit(8,'cm'),h=unit(8,"cm"))
  print(ppca_multi[[feature]])
  dev.off()
}
# dimred selection ---
npcs_sel <- list()

for(i in grep("PCA.",reducedDimNames(sce),value = T)) {
  npcs_sel[[i]] <- ncol(reducedDim(sce,i))
}
if(!is.null(batch_var) & batch_corr_method=="fastMNN") {
  dimred_method <- "corrected"
} else if(npc_sel=="min"){
  dimred_method <- names(which.min(unlist(npcs_sel)))
} else if(npc_sel=="max") {
  dimred_method <- names(which.max(unlist(npcs_sel)))
} else if(is.integer(npc_sel)) {
  reducedDim(sce, "PCA.selected") <- pcs[,1:npc_sel]
  dimred_method <- "PCA.selected"
}
# tSNE ---
set.seed(s33d)
sce <- scater::runTSNE(
  sce,
  dimred=dimred_method,
  perplexity=20
  )

# UMAP ---
set.seed(s33d)
sce <- scater::runUMAP(
  sce,
  dimred=dimred_method
  )

# Plots ---
# tSNE ---
ptsne <- list()
for(feature in metadata_features) {
  
  ptsne[[feature]] <- plotsceReducedDim(
    sce = sce,
    dimred = "TSNE",
    col_by = feature,
    shape_by = shape_by,
    point_size = 0.5,
    pal = palette_features[[feature]]
    )
  
  outfile <- paste0(FIGDIR,"/ptsne_",feature,".pdf")
  pdf(file = outfile, paper = "a4r",w=unit(3.5,'cm'),h=unit(3.5,"cm"),useDingbats = F)
  print(ptsne[[feature]])
  dev.off()
}

ptsne_list1 <- patchwork::wrap_plots(
  plotlist = ptsne,
  nrow = 1,
  ncol = NULL
)
outfile <- paste0(FIGDIR,"/ptsne_",paste0(metadata_features,collapse = "_"),".pdf")
pdf(file = outfile, paper = "a4r",w=unit(3.5*length(metadata_features),'cm'),h=unit(6,"cm"))
print(ptsne_list1)
dev.off()

# UMAP ---
pumap <- list()
for(feature in metadata_features) {
  
  pumap[[feature]] <- plotsceReducedDim(
    sce = sce,
    dimred = "UMAP",
    col_by = feature,
    shape_by = shape_by,
    point_size = 0.5,
    pal = palette_features[[feature]]
    )
  
  outfile <- paste0(FIGDIR,"/pumap_",feature,".pdf")
  pdf(file = outfile, paper = "a4r",w=unit(3.5,'cm'),h=unit(3.5,"cm"))
  print(pumap[[feature]])
  dev.off()
}

pumap_list1 <- patchwork::wrap_plots(
  plotlist = pumap,
  nrow = 1,
  ncol = NULL
)
outfile <- paste0(FIGDIR,"/pumap_",paste0(metadata_features,collapse = "_"),".pdf")
pdf(file = outfile, paper = "a4r",w=unit(3.5*length(metadata_features),'cm'),h=unit(6,"cm"))
print(pumap_list1)
dev.off()

# 7. Diagnose technical effects ----
tech_quant_feature <- c("assigned_reads","ngenes","Unmapped","Assigned_rate_total","nr_mitochondrial","nr_rRNA")
for(i in tech_quant_feature) {
  sce <- quantize_tech_feature(sce, tech_quant_feature = i)
}

techn_effects_level1 <- grep("discard_fixed_th_",colnames(colData(sce)),value = T)
ptsne_techn_level1 <- lapply(techn_effects_level1, function(x) 
  {
  plotsceReducedDim(sce = sce, dimred = "TSNE", col_by = x, shape_by = design_var, point_size=0.6) +
    ggtitle(x) + guides(col=guide_legend(title="Discard fixed threshold"))
})
names(ptsne_techn_level1) <- techn_effects_level1

plist_ptsne_techn_level1 <- ggpubr::ggarrange(plotlist = ptsne_techn_level1, ncol = 2,nrow=3, common.legend = T, legend = "right",align = "hv")

outfile <- paste0(FIGDIR,"/technical_effects_ptsne_level1.pdf")
pdf(file = outfile, paper = "a4",w=unit(8,'cm'),h=unit(8,"cm"),useDingbats = F)
print(plist_ptsne_techn_level1)
dev.off()

techn_effects_level2 <- grep("quant_",colnames(colData(sce)),value = T)
ptsne_techn_level2 <- lapply(techn_effects_level2, function(x) 
{
  pal <- RColorBrewer::brewer.pal(9,"Reds")[c(2,4,6,9)]
  plotsceReducedDim(sce = sce, dimred = "TSNE", col_by = x, shape_by = design_var, point_size=0.6,pal = pal) + 
    ggtitle(x) + guides(col=guide_legend(title="Quartile"))
})
names(ptsne_techn_level2) <- techn_effects_level2

plist_ptsne_techn_level2 <- ggpubr::ggarrange(plotlist = ptsne_techn_level2, ncol = 2,nrow=3, common.legend = T, legend = "right",align = "hv")

outfile <- paste0(FIGDIR,"/technical_effects_ptsne_level2.pdf")
pdf(file = outfile, paper = "a4",w=unit(8,'cm'),h=unit(8,"cm"),useDingbats = F)
print(plist_ptsne_techn_level2)
dev.off()

# Variance explained ---
vars <- getVarianceExplained(sce, variables=explained_var_features)

pvars <- plotExplanatoryVariables(vars) + theme_bw() + my_theme 

outfile <- paste0(FIGDIR,"/technical_effects_lmvariance.pdf")
pdf(file = outfile, paper = "a4r",w=unit(5,'cm'),h=unit(5,"cm"),useDingbats = F)
pvars
dev.off()

# Highest expressed features bias ---
phighest <- plotHighestExprs(sce, exprs_values = "counts") + theme_bw() + my_theme + xlab("counts [%]")

outfile <- paste0(FIGDIR,"/technical_effects_highestexprs.pdf")
pdf(file = outfile, paper = "a4r",w=unit(4,'cm'),h=unit(8,"cm"),useDingbats = F)
phighest
dev.off()

# 9. Cell cycle analysis ----
palette_features[["Cell_cycle"]] <- ggsci::pal_jco()(3)
names(palette_features[["Cell_cycle"]]) <- c("G1", "G2M", "S")

if(organism=="mmusculus") {
  cc.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))  
} else if(organism=="hsapiens") {
  cc.pairs <- readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))
}

assignments <- cyclone(sce, cc.pairs, gene.names=gsub("\\..*","",rowData(sce)$gene_id))

sce[["Cell_cycle"]]   <- assignments$phases
ptsne[["Cell_cycle"]] <- plotsceReducedDim(
  sce = sce,
  dimred = "TSNE",
  col_by = "Cell_cycle",
  shape_by = shape_by,
  point_size = 0.5,
  pal = palette_features[["Cell_cycle"]]
  )


outfile <- paste0(FIGDIR,"/ptsne_Cell_cycle.pdf")
pdf(file = outfile, paper = "a4r",w=unit(3.5,'cm'),h=unit(3.5,"cm"))
print(ptsne[["Cell_cycle"]])
dev.off()

pumap[["Cell_cycle"]] <- plotsceReducedDim(
  sce = sce,
  dimred = "UMAP",
  col_by = "Cell_cycle",
  shape_by = shape_by,
  point_size = 0.5,
  pal = palette_features[["Cell_cycle"]]
  )

outfile <- paste0(FIGDIR,"/pumap_Cell_cycle.pdf")
pdf(file = outfile, paper = "a4r",w=unit(3.5,'cm'),h=unit(3.5,"cm"))
print(pumap[["Cell_cycle"]])
dev.off()
# 8. Additional filters ----
# Additional required QC filters

# 9. Mark epiblast lineages ----
epi_markers <- openxlsx::read.xlsx(
  "data/xiang_et_al_human_embryo/41586_2019_1875_MOESM7_ESM.Genes_EPI_Subtypes.xlsx",
  sheet = "Formatted"
)

epi_markers <- lapply(
  split(epi_markers,epi_markers$subtype),
  "[[",
  "gene_name"
)
epi_expr_mat <- assay(
  sce[,sce$cell.type.ch1=="EPI"],
  "logcounts"
  )
epi_rankings <- AUCell::AUCell_buildRankings(
  epi_expr_mat
)
epi_markers_auc <- AUCell::AUCell_calcAUC(
  geneSets = epi_markers,
  rankings = epi_rankings
)

# epi_cells_assignment <- AUCell::AUCell_exploreThresholds(
#   epi_markers_auc,
#   plotHist=TRUE,
#   assign=TRUE
# )

sce$CellType <- sce$cell.type.ch1
sce$CellType[
  sce$Sample %in% names(which(epi_markers_auc@assays@data$AUC["Post-EPIs",]>epi_markers_auc@assays@data$AUC["Pre-EPIs",]))
] <- "Post-EPIs"
sce$CellType[
  sce$Sample %in% names(which(epi_markers_auc@assays@data$AUC["Pre-EPIs",]>epi_markers_auc@assays@data$AUC["Post-EPIs",]))
] <- "Pre-EPIs"

ptsne <- list()
metadata_features <- c("CellType","age.ch1","cell.type.ch1")
for(feature in metadata_features) {
  
  ptsne[[feature]] <- plotsceReducedDim(
    sce = sce,
    dimred = "TSNE",
    col_by = feature,
    shape_by = shape_by,
    point_size = 0.5,
    pal = palette_features[[feature]]
  )
}

outfile <- paste0(FIGDIR,"/ptsne_CellType.pdf")
pdf(file = outfile, paper = "a4r",w=unit(3.5,'cm'),h=unit(3.5,"cm"),useDingbats = F)
print(ptsne[["CellType"]])
dev.off()

ptsne_list1 <- patchwork::wrap_plots(
  plotlist = ptsne,
  nrow = 1,
  ncol = NULL
)
outfile <- paste0(FIGDIR,"/ptsne_",paste0(metadata_features,collapse = "_"),".pdf")
pdf(file = outfile, paper = "a4r",w=unit(3.5*length(metadata_features),'cm'),h=unit(6,"cm"))
print(ptsne_list1)
dev.off()
# 9. Save pre-processed data ----
outfile <- paste0(RESDIR,"/sce.rds")
message(" -- saving: ",outfile)
saveRDS(sce, file = outfile)

if(exists("sce_filtered")) {
  outfile <- paste0(RESDIR,"/sce_filtered.rds")
  message(" -- saving: ",outfile)
  saveRDS(sce_filtered, file = outfile)
}

if(save_additional_files) {
  unfiltered_metadata <- colData(unfiltered)
  outfile <- paste0(RESDIR,"/unfiltered_metadata.txt")
  message(" -- saving: ",outfile)
  write.table(unfiltered_metadata, file = outfile, sep = "\t",row.names = F, quote = F)
  
  sce_metadata_results <- colData(sce)
  outfile <- paste0(RESDIR,"/sce_metadata_results.txt")
  message(" -- saving: ",outfile)
  write.table(sce_metadata_results, file = outfile, sep = "\t",row.names = F, quote = F)
  
  sce_normalized_counts <- as.data.frame(assay(sce,assay.type))
  sce_normalized_counts <- cbind.data.frame("gene_name"=rownames(sce_normalized_counts),sce_normalized_counts)
  outfile <- paste0(RESDIR,"/sce_normalized_counts.",assay.type,".txt")
  message(" -- saving: ",outfile)
  write.table(sce_normalized_counts, file = outfile, sep = "\t",row.names = F, quote = F)
}

