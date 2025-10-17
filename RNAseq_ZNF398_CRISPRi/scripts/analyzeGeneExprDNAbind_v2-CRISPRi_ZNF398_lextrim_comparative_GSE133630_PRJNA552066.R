# setup ------
setwd("/home/epigen/Martello_3D_human_epiblast/RNAseq/dataset/v2-CRISPRi_ZNF398_lextrim")

library(RNAseqRtools)

path_binding_info <- "/home/epigen/Martello_3D_human_epiblast/ChIPseq/bit_encode-chip-seq-pipeline/dataset/comparative_GSE133630_PRJNA552066/peak/ChIP-seq_of_ZNF398_in_BG01V_human_embryonic_stem_cells-no_control.ChIP-seq_of_ZNF398_in_H9_human_embryonic_stem_cells-no_control.idr.bchrfilt.binding_info.gz"
path_dge_toptable <- "/home/epigen/Martello_3D_human_epiblast/RNAseq/dataset/v2-CRISPRi_ZNF398_lextrim/DGE_sep_group/edger.toptable_clean.ALL_contrast.mark_seqc.header_added.gz"
path_normcounts <- "results/normcounts_tmm.logcpm.txt.gz"
path_metadata <- "metadata.txt"

path_results <- "results"
# read data ------
binding_info <- read.delim(path_binding_info)
dge_toptable <- read.delim(path_dge_toptable)
normcounts <- read.delim(path_normcounts)
metadata <- read.delim(path_metadata)

# split by comparison ----
binding_info <- binding_info[-which(is.na(binding_info$gene)),]
binding_and_dge <- lapply(
  split(dge_toptable, dge_toptable$contrast),
  function(x) {
    merge(x, binding_info, by.x="GeneID", by.y="gene", all.x=T)
  }
)

# heatmaps ----
pal_analysis <- c(
  "P_NI_2D_Day7" = "#FFCCFF",
  "ZNF398_NI_2D_Day7" = "#C6DBEF",
  "ZNF398_I_2D_Day7" = "#D94801" ,
  "ZNF398_NI_3D_Day4" = "#2171B5",
  "ZNF398_I_3D_Day4" = "#7F2704" 
)
hm_bd_list <- list()
for(i in names(binding_and_dge)) {
  bd_i <- binding_and_dge[[i]]
  
  bd_i <- bd_i[bd_i$significance!=0,]
  expr_mat <- normcounts[normcounts$gene_name %in% bd_i$GeneID,]
  rownames(expr_mat) <- expr_mat$gene_name
  expr_mat <- expr_mat[,-c(1:3)]
  expr_mat <- as.matrix(expr_mat)
  groups <- c(gsub("(.*)_vs_(.*)","\\1",unique(bd_i$contrast)),
              gsub("(.*)_vs_(.*)","\\2",unique(bd_i$contrast)))
  samples <- metadata$sample[metadata$condition_sepgr %in% groups]
  samples <- gsub("\\-","\\.",samples)
  expr_mat <- expr_mat[,samples]
  
  expr_hm <- get_heatmap4(
    m = expr_mat,
    show_row_dend = F,
    show_row_names = F,
    bm = F,
    myLegend = "logCPM",
    myPalette = c("#54278F","white","#D94801"),
    myZscale = c(-1.5,0,1.5),
    retHm = T,
    annotDF = data.frame("group"=metadata$condition_sepgr[match(colnames(expr_mat),gsub("\\-","\\.",metadata$sample))]),
    annotCol = list("group"=pal_analysis),
    annotation_width_heigth = .7,
    show_annotation_name = F,
    split = 2,
    row_title = c("down-regulated", "up-regulated"),
    row_title_gp = gpar(fontsize = 8),
    cluster_columns = F,
    column_title = "ZNF398 CRISPRi RNAseq",
  )
  
  binding_mat <- matrix(
    0,
    nrow = nrow(expr_mat),
    ncol = 2,
    dimnames = list(rownames(expr_mat), c("PLS","ELS"))
  )
  bd_map <- unique(bd_i[,c("GeneID","ccRE")])
  bd_map <- ddply(bd_map, .(GeneID), summarize, ccRE=paste0(ccRE, collapse=","))
  bd_map <- bd_map[match(rownames(expr_mat), bd_map$GeneID),]
  print(nrow(bd_map) == nrow(binding_mat))
  for(j in colnames(binding_mat)) {
    binding_mat[,j] <- grepl(j, bd_map$ccRE)
  }
  binding_hm <- Heatmap(
    binding_mat,
    col=c('white','darkgreen'),
    show_row_dend = F,
    show_row_names = F,
    cluster_columns = F,
    border = T,
    column_title = "ZNF398 DNA binding",
    column_names_gp = gpar(fontsize=8),
    column_title_gp = gpar(fontsize=8),
    show_heatmap_legend = F,
  )
  hm_bd <- expr_hm + binding_hm
  
  outfile  <- paste0(path_results, "/hm_bd.",i,".pdf")
  message('-- saving to: ',outfile)
  pdf(file = outfile, paper = "a4", w = unit(4.2,'cm'), h = unit(4.8,'cm'))
  draw(hm_bd, heatmap_legend_side="bottom")
  dev.off()
  
  hm_bd_list[[i]] <- draw(hm_bd, heatmap_legend_side="bottom")
}


# Markers expression ----
y <- readRDS("results/edger.processed.counts.rds")
y$samples$condition_sepgr <- factor(y$samples$condition_sepgr, levels=names(pal_analysis))
y$samples$sample <- gsub("\\-","\\.",y$samples$sample)

plt.expr <- list()
markers.list <- list(
  "Pluripotency" = c(
    "NANOG","POU5F1","SOX2","ZNF398","PRDM14"
  ),
  "Primitive streak" = c(
    "MIXL1","EOMES","GATA6","KDR","SOX17","TBXT"
  ),
  "TGF-beta" = c(
    "LEFTY1","LEFTY2","SKIL","SMAD7","SMAD3"
  ),
  "Epithelial/Polarity" = c(
    "CDH1","ESRP1","EPCAM","CLDN7","PODXL"
  ),
  "Mesenchymal" = c(
    "CDH2","SNAI1","VIM","FN1","SNAI2"
  )
)

for(i in names(markers.list)) {
  genes <- intersect(markers.list[[i]],rownames(y))
  expr.mat <- y$CPM
  expr.mat <- expr.mat[genes,]
  de.2d.idx <- genes %in% dge_toptable[dge_toptable$contrast=="ZNF398_I_2D_Day7_vs_ZNF398_NI_2D_Day7" & dge_toptable$significance!=0,'GeneID']
  de.3d.idx <- genes %in% dge_toptable[dge_toptable$contrast=="ZNF398_I_3D_Day4_vs_ZNF398_NI_3D_Day4" & dge_toptable$significance!=0,'GeneID']
  bound.idx <- genes %in% binding_info$gene
  genes.lab <- genes
  names(genes.lab) <- genes
  genes.lab[de.2d.idx] <- paste0(genes[de.2d.idx]," *")
  genes.lab[de.3d.idx] <- paste0(genes.lab[de.3d.idx]," °")
  genes.lab[bound.idx] <- paste0(genes.lab[bound.idx]," ^")
  rownames(expr.mat) <- genes.lab[rownames(expr.mat)]
  
  plt.expr[[i]] <- plotExpression(
    expr.mat,
    experimental_info = y$samples,
    gene = rownames(expr.mat),
    group.by = "condition_sepgr",
    pal = pal_analysis,
    expression.unit = "CPM",
    show.names = T,
    names.rot = 45,
    facet_ncol = length(markers.list[[i]]),
  ) +
    ggtitle(i) +
    theme(
      plot.title = element_text(size=8, face="plain"),
      strip.text = element_text(size=8, face="plain")
    )
}

plt.expr.arranged <- patchwork::wrap_plots(
  plt.expr,
  ncol = 1
) + patchwork::plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

figout  <- paste0(path_results, "/plt.expr.marker_genes.bound.pdf")
pdf(file = figout, paper = "a4", w = unit(8.5,'cm'), h = unit(11,'cm'))
print(plt.expr.arranged)
dev.off()

# Markers expression 2 ----
y <- readRDS("results/edger.processed.counts.rds")
y$samples$condition_sepgr <- factor(y$samples$condition_sepgr, levels=names(pal_analysis))
y$samples$sample <- gsub("\\-","\\.",y$samples$sample)

plt.expr <- list()

# Create a data frame for epithelial marker genes
epithelial_markers <- data.frame(
  Category = c(
    # Cytoskeletal and Structural Markers
    rep("Cytoskeletal and Structural Markers", 6),
    # Adhesion Molecules
    rep("Adhesion Molecules", 6),
    # Polarity Regulators
    rep("Polarity Regulators", 4),
    # Transcription Factors
    rep("Transcription Factors", 4),
    # Basement Membrane and ECM Markers
    rep("Basement Membrane and ECM Markers", 4),
    # Functional and Transport Markers
    rep("Functional and Transport Markers", 3)
  ),
  Gene = c(
    # Cytoskeletal and Structural Markers
    "KRT8", "KRT18", "KRT5", "KRT14", "KRT19", "EPCAM",
    # Adhesion Molecules
    "CDH1", "CLDN1", "CLDN3", "CLDN4", "OCLN", "DSG1",
    # Polarity Regulators
    "CRB3", "PARD3", "LGL1", "CRB3",
    # Transcription Factors
    "GRHL2", "FOXA1", "GATA3", "OVOL1",
    # Basement Membrane and ECM Markers
    "LAMA1", "LAMB1", "COL4A1", "COL4A2",
    # Functional and Transport Markers
    "MUC1", "SLC2A1", "ABCG2"
  )
)
markers.list <- lapply(split(epithelial_markers, epithelial_markers$Category),"[[",2)


for(i in names(markers.list)) {
  genes <- intersect(markers.list[[i]],rownames(y))
  expr.mat <- y$CPM
  expr.mat <- expr.mat[genes,]
  de.2d.idx <- genes %in% dge_toptable[dge_toptable$contrast=="ZNF398_I_2D_Day7_vs_ZNF398_NI_2D_Day7" & dge_toptable$significance!=0,'GeneID']
  de.3d.idx <- genes %in% dge_toptable[dge_toptable$contrast=="ZNF398_I_3D_Day4_vs_ZNF398_NI_3D_Day4" & dge_toptable$significance!=0,'GeneID']
  bound.idx <- genes %in% binding_info$gene
  genes.lab <- genes
  names(genes.lab) <- genes
  genes.lab[de.2d.idx] <- paste0(genes[de.2d.idx]," *")
  genes.lab[de.3d.idx] <- paste0(genes.lab[de.3d.idx]," °")
  genes.lab[bound.idx] <- paste0(genes.lab[bound.idx]," ^")
  rownames(expr.mat) <- genes.lab[rownames(expr.mat)]
  
  plt.expr[[i]] <- plotExpression(
    expr.mat,
    experimental_info = y$samples,
    gene = rownames(expr.mat),
    group.by = "condition_sepgr",
    pal = pal_analysis,
    expression.unit = "CPM",
    show.names = T,
    names.rot = 45,
    facet_ncol = length(markers.list[[i]]),
  ) +
    ggtitle(i) +
    theme(
      plot.title = element_text(size=8, face="plain"),
      strip.text = element_text(size=8, face="plain")
    )
}

plt.expr.arranged <- patchwork::wrap_plots(
  plt.expr,
  ncol = 1
) + patchwork::plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

figout  <- paste0(path_results, "/plt.expr.marker_genes.bound.2.pdf")
pdf(file = figout, paper = "a4", w = unit(8.5,'cm'), h = unit(11,'cm'))
print(plt.expr.arranged)
dev.off()

# Selected for figure----
source("/home/epigen/Martello_3D_human_epiblast/public/figures/scripts/resources.R")
sample.id <- y$samples$sample[y$samples$condition_sepgr!="P_NI_2D_Day7"]
y$samples$condition_sepgr_fix <- y$samples$condition_sepgr
levels(y$samples$condition_sepgr_fix) <- gsub("ZNF398_NI_2D_Day7","Ctrl (2D)",levels(y$samples$condition_sepgr_fix))
levels(y$samples$condition_sepgr_fix) <- gsub("ZNF398_I_2D_Day7","KD (2D)",levels(y$samples$condition_sepgr_fix))
levels(y$samples$condition_sepgr_fix) <- gsub("ZNF398_NI_3D_Day4","Ctrl (3D)",levels(y$samples$condition_sepgr_fix))
levels(y$samples$condition_sepgr_fix) <- gsub("ZNF398_I_3D_Day4","KD (3D)",levels(y$samples$condition_sepgr_fix))
names(pal_analysis) <- levels(y$samples$condition_sepgr_fix)
plt.expr.sel <- plotExpression(
  y$CPM[,sample.id],
  experimental_info = y$samples[sample.id,],
  gene = c("ESRP1","PODXL"),
  group.by = "condition_sepgr_fix",
  pal = pal_analysis,
  expression.unit = "CPM",
  show.names = T,
  names.rot = 45,
  facet_ncol = length(markers.list[[i]]),
) +
  theme_classic() +
  my_theme_barplt +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

figout  <- paste0(path_results, "/plt.expr.sel.pdf")
pdf(file = figout, paper = "a4", w = unit(3.65,'cm'), h = unit(1.84,'cm'))
print(plt.expr.sel)
dev.off()

# GSEA bound genes ----
library(fgsea)
unique(dge_toptable$contrast)
i = "ZNF398_I_3D_Day4_vs_ZNF398_NI_3D_Day4"
de_i <- dge_toptable[dge_toptable$contrast==i,]
de_i <- de_i[,-1]
colnames(de_i) <- gsub("Pvalue","PValue", colnames(de_i))
de_rank_i <- getRankingMetric(de_i)
bound_genes <- list("dna_bound"=unique(binding_info$gene))
tmp <- fgsea(bound_genes,de_rank_i)

# Binding by cell state ----
dge_summary <- read.delim("results/dge.summary.gz")
dge_summary$bound <- dge_summary$gene_name %in% binding_info$gene

dge_summary_bound_count <- sapply(
  colnames(dge_summary)[-c(1,8)],
  function(i) {
    sum(rowSums(dge_summary[,c(i, "bound")])==2)
  }
)
dge_summary_bound_frac <- dge_summary_bound_count/colSums(dge_summary[,-c(1,8)])

