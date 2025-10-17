setwd("/home/epigen/Martello_3D_human_epiblast/Martello_Polo_MIXL1/run240820_run241003/analysis_oscab")
# Analyze scRNAseq - Bioconductor workflow
# Pre-processing
# 0. Resources ----
source("scripts/resources.R")

# paths ---
path_metadata = "data/metadata.txt.gz"
path_counts = "data/raw_counts.csv.gz"
path_gene_info = "data/gene_info.gencode.hsapiens.32.primary_assembly.annotation.txt.gz"

# output dirs ---
path_results  = "results"

if(!dir.exists(path_results)) dir.create(path_results, recursive = T)

# Global params ---
metadata_features <- c("Plate","Time","batch")
organism <- "hsapiens" # hsapiens | mmusculus

assay.type  <- "logcounts"
norm_method <- "pool" # pool | pool_cluster | libsize
norm_filter <- NULL # rRNA_mito | NULL
npc_sel     <- "max" # "min", "max"

batch_var  <- NULL # "batch"
design_var <- "Time"
batch_corr_method <- "fastMNN" # limma | regressBatches | rescaleBatches | fastMNN
batch_corr_use_design <- FALSE
shape_by <- NULL

explained_var_features <- c(
  metadata_features,
  "Assigned_rate_total",
  "assigned_reads",
  "perc_mitochondrial",
  "ngenes"
  )

additional_filters <- FALSE
save_additional_files <- FALSE
s33d <- 1991
# 1. Load and prepare data -----
# Load data ---
counts <- read.delim(
  path_counts,
  stringsAsFactors = F,
  header = T,
  sep = ",",
  row.names = 1
  )
colnames(counts) <- gsub("^X","",colnames(counts))

# metadata
metadata <- read.delim(
  path_metadata,
  stringsAsFactors = F,
  header = T
  )
rownames(metadata) <- metadata[,1]
metadata <- metadata[colnames(counts),]
message("cell metadata names check: ",
        sum(rownames(metadata) == colnames(counts)) == ncol(counts)
)
metadata$Plate <- factor(metadata$Plate)

# gene metadata
gene_info  <- read.delim(
  path_gene_info,
  stringsAsFactors = F,
  header = T
)
rownames(gene_info) <- gene_info[,1]
gene_info <- gene_info[rownames(counts),]
message("gene metadata names check: ",
        sum(rownames(gene_info) == rownames(counts)) == nrow(counts)
)

# Create SCE object ---
sce <- SingleCellExperiment::SingleCellExperiment(
  assays = list(counts = as.matrix(counts)),
  colData = metadata,
  rowData=gene_info
  )

# Filter empty|minibulk|zero-read cells
if(length(which(sce$assigned_reads==0))!=0) sce <- sce[,-which(sce$assigned_reads==0)]
if(length(grep("MBULK|EMPTY|MINIBULK",colnames(sce))!=0)) sce <- sce[,-grep("MBULK|EMPTY|MINIBULK",colnames(sce))]

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
qc_df <- scater::perCellQCMetrics(sce, subsets=list(Mito=is.mito), percent.top=50)
sce <- scater::addPerCellQC(sce, subsets=list(Mito=is.mito), percent.top=50)

qc.assigned <- sce$Assigned_rate_total < 20
qc.lib <- qc_df$sum < 3e4#1e5
qc.nexprs <- qc_df$detected < 5e2#5e3
# qc.spike <- qc_df$altexps_ERCC_percent > 10
qc.mito <- qc_df$subsets_Mito_percent > 75 #10
qc.mito[which(is.na(qc.mito))] <- FALSE
qc.top <- qc_df$percent.top_50 > 90 #10
qc.top[which(is.na(qc.top))] <- FALSE

# default: qc.lib | qc.nexprs | qc.spike | qc.mito | qc.top | qc.assigned
discard_fixed_th <- qc.lib | qc.nexprs | qc.mito | qc.top | qc.assigned
# custom
# discard_fixed_th <- qc.lib | qc.nexprs 

fixed_th_df <- data.frame(
  discard_fixed_th_LibSize=qc.lib,
  discard_fixed_th_NExprs=qc.nexprs,
  discard_fixed_th_Mito_percent=qc.mito,
  discard_fixed_th_percent.top_50=qc.top,
  discard_fixed_th_assigned_rate=qc.assigned,
  discard_fixed_th=discard_fixed_th
  )
rownames(fixed_th_df) <- rownames(qc_df)

discard_fixed_th_df_stats <- data.frame(
  LibSize=sum(qc.lib),
  NExprs=sum(qc.nexprs),
  # SpikeProp=sum(qc.spike),
  MitoProp=sum(qc.mito),
  TopProp=sum(qc.top),
  Total_cells=length(discard_fixed_th),
  Total_discard=sum(discard_fixed_th),
  Total_keep=sum(!discard_fixed_th)
  )
colData(sce) <- cbind.DataFrame(colData(sce),fixed_th_df)

qc_fixed_th_plots <- list()
for(feature in metadata_features) {
  qc_fixed_th_plots[[feature]] <- plot_qc_th(
    sce,
    x_var = feature,
    discard_var = "discard_fixed_th"
    ) 
  outfile <- paste0(path_results,"/QC_plots_fixed_th_",feature,".pdf")
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
# adaptive_th_df <- scater::quickPerCellQC(qc_df, percent_subsets = c("subsets_Mito_percent","percent.top_50"))
adaptive_th_df <- scater::quickPerCellQC(qc_df)
colnames(adaptive_th_df) <- gsub("discard","discard_adaptive_th",colnames(adaptive_th_df))
discard_adaptive_th_df_stats <- colSums(as.matrix(adaptive_th_df))
colData(sce) <- cbind.DataFrame(colData(sce),adaptive_th_df)

qc_adaptive_th_plots_time_treatment <- plot_qc_th(sce, x_var = "Time",discard_var = "discard_adaptive_th")

qc_adaptive_th_plots <- list()
for(feature in metadata_features) {
  qc_adaptive_th_plots[[feature]] <- plot_qc_th(
    sce,
    x_var = feature,
    discard_var = "discard_adaptive_th")  
  outfile <- paste0(path_results,"/QC_plots_MAD_adaptive_th_",feature,".pdf")
  message(" -- saving figure: ", outfile)
  pdf(file = outfile, paper = "a4",w=unit(6,'cm'),h=unit(10,"cm"),useDingbats = F)
  print(qc_adaptive_th_plots[[feature]])
  dev.off()
}

# 2.3 Discard -----
discard_var <- "discard_adaptive_th"
sce$discard <- sce[[discard_var]]
unfiltered  <- sce # save unfiltered copy
table(unfiltered$discard)
sce <- sce[,!sce$discard]

# The hard way
table(!(unfiltered$assigned_reads<1e5 | unfiltered$ngenes < 2000 | unfiltered$perc_mitochondrial >= 25))
# 3. Normalization ----
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

size_factors_df <- data.frame(
  "librarySizeFactors"=librarySizeFactors(sce),
  "DeconvolutionFactors"=sizeFactors(sce)
  )
p_size_factor <- list()
for(feature in metadata_features) {
  size_factors_df[[feature]] <- colData(sce)[,feature]
  p_size_factor[[feature]] <- ggplot(
    size_factors_df,
    aes(
      x = librarySizeFactors,
      y = DeconvolutionFactors,
      col = .data[[feature]]
      )
    ) +
    ggrastr::rasterise(geom_point(size=1)) +
    theme_bw() + my_theme + 
    scale_x_log10() + scale_y_log10() + 
    scale_color_manual(values = palette_features[[feature]])
  
  outfile <- paste0(path_results,"/size_factors_col_feature_",feature,".pdf")
  pdf(file = outfile, paper = "a4",w=unit(3,'cm'),h=unit(3,"cm"),useDingbats = F)
  print(p_size_factor[[feature]])
  dev.off()
}

# 4. Variance Modeling for feature selection ----
if(is.null(batch_var)) {
  if(is.null(norm_filter)) {
    dec <- modelGeneVar(sce)  
  } else {
    dec <- modelGeneVar(sce, subset.row=rowData(sce)$use_for_norm)
  }
  
  pdec <- ggplot(as.data.frame(dec),aes(x=mean,y=total)) + 
    ggrastr::rasterise(geom_point(size=0.5),dpi=300) +
    geom_line(aes(x=mean,y=tech),col="red") +
    theme_bw() + my_theme + xlab("Mean of log-expression")+
    ylab("Variance of log-expression")
  
  outfile <- paste0(path_results,"/variance_modeling.pdf")
  pdf(file = outfile, paper = "a4",w=unit(3,'cm'),h=unit(3,"cm"),useDingbats = F)
  print(pdec)
  dev.off()
  
} else {
  if(is.null(norm_filter)) {
    dec <- modelGeneVar(sce, block=sce[[batch_var]],subset.row=rowData(sce)$use_for_norm) # block for batch  
  } else {
    dec <- modelGeneVar(sce, block=sce[[batch_var]]) # block for batch
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
  
  outfile <- paste0(path_results,"/variance_modeling.pdf")
  pdf(file = outfile, paper = "a4r",w=unit(12,'cm'),h=unit(h_value,"cm"),useDingbats = F)
  print(ggpubr::ggarrange(plotlist = pdec, nrow = arrange_nrow,ncol=arrange_ncol))
  dev.off()
}
outfile  <- paste0(path_results,"/dec.rds")
message(" -- saving: ",outfile)
saveRDS(dec, file = outfile)

dec_id <- intersect(rownames(dec)[dec$mean>=1], names(which(rowSums(counts(sce)>=25)>=1)))
# Get the top 10% of genes.
top.hvgs <- getTopHVGs(dec[dec_id,], prop=0.1)

# Get the top 2500 genes.
top.hvgs2 <- getTopHVGs(dec[dec_id,], n=2500)

# Get all genes with positive biological components.
top.hvgs3 <- getTopHVGs(dec[dec_id,], var.threshold=0)

# Get all genes with FDR below 5%.
top.hvgs4 <- getTopHVGs(dec[dec_id,], fdr.threshold=0.05)

vargenes <- top.hvgs2
varmodel <- dec[vargenes,]

outfile  <- paste0(path_results,"/variance_modeling.xlsx")
varmodel <- cbind("gene_name"=rownames(varmodel),varmodel)
saveXLSresEdgeR2(varmodel, outfile = outfile)

dec_top <- as.data.frame(dec)
dec_top$top <- FALSE
dec_top[vargenes,"top"] <- TRUE
pdec_top <- ggplot(dec_top,aes(x=mean,y=total,col=top)) + 
  ggrastr::rasterise(geom_point(size=0.5),dpi=300) +
  geom_line(aes(x=mean,y=tech),col="red") +
  theme_bw() + my_theme + xlab("Mean of log-expression")+
  ylab("Variance of log-expression") + scale_color_manual(values=c("black","#2171B5"))

outfile <- paste0(path_results,"/variance_modeling_top.pdf")
pdf(file = outfile, paper = "a4",w=unit(3,'cm'),h=unit(3,"cm"),useDingbats = F)
print(pdec_top)
dev.off()

rowData(sce)$vargenes <- FALSE
rowData(sce)$vargenes[match(vargenes,rowData(sce)$gene_name)] <- TRUE
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
  subset_row=vargenes,
  exprs_values=assay.type,
  BSPARAM=BiocSingular::ExactParam()
  )

# Number of PCs - Elbow method ---
percent.var <- attr(SingleCellExperiment::reducedDim(sce), "percentVar")
chosen.elbow <- PCAtools::findElbowPoint(percent.var)
SingleCellExperiment::reducedDim(sce, "PCA.elbow") <- SingleCellExperiment::reducedDim(sce)[,1:chosen.elbow]
reducedDimNames(sce)

# Number of PCs - Population structure: clustered PCs ---
pcs <- SingleCellExperiment::reducedDim(sce, "PCA")
choices <- getClusteredPCs(pcs)
val <- metadata(choices)$chosen
SingleCellExperiment::reducedDim(sce, "PCA.clust") <- pcs[,1:val]

ppca <- list()
ppca_multi <- list()
for(feature in metadata_features) {
  ppca[[feature]] <- plotsceReducedDim(
    sce = sce,
    dimred = "PCA",
    col_by = feature,
    shape_by = shape_by,
    point_size = 0.6,
    pal = palette_features[[feature]]
    )
  
  outfile <- paste0(path_results,"/ppca_",feature,".pdf")
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
  
  outfile <- paste0(path_results,"/ppca_multi_",feature,".pdf")
  pdf(file = outfile, paper = "a4r",w=unit(8,'cm'),h=unit(8,"cm"),useDingbats = F)
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
perplexity <- 25
set.seed(s33d)
sce <- scater::runTSNE(
  sce,
  dimred=dimred_method,
  perplexity=perplexity
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
    point_size = 0.6,
    pal = palette_features[[feature]]
    )
  
  outfile <- paste0(path_results,"/ptsne_",feature,".pdf")
  pdf(file = outfile, paper = "a4r",w=unit(3.5,'cm'),h=unit(3.5,"cm"))
  print(ptsne[[feature]])
  dev.off()
}
ptsne_list1 <- ggpubr::ggarrange(plotlist = ptsne, nrow=1, align = "hv")

outfile <- paste0(path_results,"/ptsne_",paste0(metadata_features,collapse = "_"),".pdf")
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
    point_size = 0.6,
    pal = palette_features[[feature]]
  )
  
  outfile <- paste0(path_results,"/pumap_",feature,".pdf")
  pdf(file = outfile, paper = "a4r",w=unit(3.5,'cm'),h=unit(3.5,"cm"))
  print(pumap[[feature]])
  dev.off()
}
pumap_list1 <- ggpubr::ggarrange(plotlist = pumap, nrow=1, align = "hv")

outfile <- paste0(path_results,"/pumap_",paste0(metadata_features,collapse = "_"),".pdf")
pdf(file = outfile, paper = "a4r",w=unit(3.5*length(metadata_features),'cm'),h=unit(6,"cm"))
print(pumap_list1)
dev.off()

# 7. Diagnose technical effects ----
tech_quant_feature <- c("assigned_reads","ngenes","Unmapped","Assigned_rate_total","nr_mitochondrial")
for(i in tech_quant_feature) {
  print(i)
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

outfile <- paste0(path_results,"/technical_effects_ptsne_level1.pdf")
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

outfile <- paste0(path_results,"/technical_effects_ptsne_level2.pdf")
pdf(file = outfile, paper = "a4",w=unit(8,'cm'),h=unit(8,"cm"),useDingbats = F)
print(plist_ptsne_techn_level2)
dev.off()

# Variance explained ---
vars <- getVarianceExplained(sce, variables=explained_var_features)

pvars <- plotExplanatoryVariables(vars) + theme_bw() + my_theme 

outfile <- paste0(path_results,"/technical_effects_lmvariance.pdf")
pdf(file = outfile, paper = "a4r",w=unit(5,'cm'),h=unit(5,"cm"),useDingbats = F)
pvars
dev.off()

# Highest expressed features bias ---
phighest <- plotHighestExprs(sce, exprs_values = "counts") + theme_bw() + my_theme + xlab("counts [%]")

outfile <- paste0(path_results,"/technical_effects_highestexprs.pdf")
pdf(file = outfile, paper = "a4r",w=unit(4,'cm'),h=unit(8,"cm"),useDingbats = F)
phighest
dev.off()

# 9. Cell cycle analysis ----
palette_features[["Cell_cycle"]] <- ggsci::pal_jco()(3)
names(palette_features[["Cell_cycle"]]) <- c("G1", "G2M", "S")

if(organism=="mmusculus") {
  cc.pairs <- readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))  
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
  point_size = 0.6,
  pal = palette_features[["Cell_cycle"]]
  )


outfile <- paste0(path_results,"/ptsne_Cell_cycle.pdf")
pdf(file = outfile, paper = "a4r",w=unit(5,'cm'),h=unit(6,"cm"),useDingbats = F)
print(ptsne[["Cell_cycle"]])
dev.off()

pumap[["Cell_cycle"]] <- plotsceReducedDim(
  sce = sce,
  dimred = "UMAP",
  col_by = "Cell_cycle",
  shape_by = shape_by,
  point_size = 0.6,
  pal = palette_features[["Cell_cycle"]]
  )


outfile <- paste0(path_results,"/pumap_Cell_cycle.pdf")
pdf(file = outfile, paper = "a4r",w=unit(5,'cm'),h=unit(6,"cm"),useDingbats = F)
print(pumap[["Cell_cycle"]])
dev.off()
# 8. Additional filters ----
# Additional required QC filters

# 9. Save pre-processed data ----
outfile <- paste0(path_results,"/sce.rds")
message(" -- saving: ",outfile)
saveRDS(sce, file = outfile)

if(exists("sce_filtered")) {
  outfile <- paste0(path_results,"/sce_filtered.rds")
  message(" -- saving: ",outfile)
  saveRDS(sce_filtered, file = outfile)
}

if(save_additional_files) {
  unfiltered_metadata <- colData(unfiltered)
  outfile <- paste0(path_results,"/unfiltered_metadata.txt")
  message(" -- saving: ",outfile)
  write.table(unfiltered_metadata, file = outfile, sep = "\t",row.names = F, quote = F)
  
  sce_metadata_results <- colData(sce)
  outfile <- paste0(path_results,"/sce_metadata_results.txt")
  message(" -- saving: ",outfile)
  write.table(sce_metadata_results, file = outfile, sep = "\t",row.names = F, quote = F)
  
  sce_normalized_counts <- as.data.frame(assay(sce,assay.type))
  sce_normalized_counts <- cbind.data.frame("gene_name"=rownames(sce_normalized_counts),sce_normalized_counts)
  outfile <- paste0(path_results,"/sce_normalized_counts.",assay.type,".txt")
  message(" -- saving: ",outfile)
  write.table(sce_normalized_counts, file = outfile, sep = "\t",row.names = F, quote = F)
}
