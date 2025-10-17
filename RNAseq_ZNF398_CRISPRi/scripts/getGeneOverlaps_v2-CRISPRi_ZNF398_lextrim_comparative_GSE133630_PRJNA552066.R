# setup ------
setwd("/home/epigen/Martello_3D_human_epiblast/RNAseq/dataset/v2-CRISPRi_ZNF398_lextrim")

library(RNAseqRtools)
library(GeneOverlap)
library(enrichR)

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

# generate gene lists ----
binding_info <- binding_info[-which(is.na(binding_info$gene)),]
# binding_info <- binding_info[which(abs(binding_info$distTSS)<=1e3),]
gene_lists <- list(
  "dna_bound"=unique(binding_info$gene),
  "up_reg_2d"=dge_toptable[dge_toptable$contrast=="ZNF398_I_2D_Day7_vs_ZNF398_NI_2D_Day7" & dge_toptable$significance==1,'GeneID'],
  "down_reg_2d"=dge_toptable[dge_toptable$contrast=="ZNF398_I_2D_Day7_vs_ZNF398_NI_2D_Day7" & dge_toptable$significance==-1,'GeneID'],
  "up_reg_3d"=dge_toptable[dge_toptable$contrast=="ZNF398_I_3D_Day4_vs_ZNF398_NI_3D_Day4" & dge_toptable$significance==1,'GeneID'],
  "down_reg_3d"=dge_toptable[dge_toptable$contrast=="ZNF398_I_3D_Day4_vs_ZNF398_NI_3D_Day4" & dge_toptable$significance==-1,'GeneID']
  )

gene_lists_len <- unlist(lapply(gene_lists, length))
# gene overlap ----
gene_overlap_mat <- GeneOverlap::newGOM(gene_lists, spec = "hg19.gene")

# Visualize gene overlap: heatmap DNA bound ----
overlap.pvalue.mat <- getMatrix(gene_overlap_mat,'pval')[1, drop=F,]
overlap.intersect.mat <- getMatrix(gene_overlap_mat,'intersection')[1, drop=F,]
overlap.int.frac.mat <- overlap.intersect.mat/gene_lists_len[colnames(overlap.intersect.mat)]

column_ha = HeatmapAnnotation(
  frac = anno_barplot(overlap.int.frac.mat[1,],
                      bar_width = .8,
                      gp = gpar(fill='skyblue'),
                      height = unit(1, "cm"),
                      axis_param = list(gp=gpar(fontsize=6))),
  annotation_name_gp = gpar(fontsize = 7),
  annotation_name_side = 'left',
  annotation_name_rot = 90,
  annotation_label = 'bound\ngenes\n(%)'
  )

overlap.pvalue.mat <- signif(overlap.pvalue.mat,digits = 3)
overlap.pvalue.mat[overlap.pvalue.mat<2.2e-16] <- 2.2e-16
overlap.pvalue.mat.draw <- overlap.intersect.mat
# overlap.pvalue.mat.draw[overlap.pvalue.mat.draw==2.2e-16] <- "< 2.2e-16"

mat.to.plot <- -log10(overlap.pvalue.mat[1,,drop=F])
colnames(mat.to.plot) <- gsub("_reg_"," in ",colnames(mat.to.plot))
rownames(mat.to.plot) <- "DNA bound"
overlap.pvalue.hm <- ComplexHeatmap::Heatmap(
  mat.to.plot,
  col = RColorBrewer::brewer.pal(8,"Oranges"),
  show_row_dend = T,
  show_row_names = T,
  cluster_columns = T,
  row_names_side = "left",
  row_names_gp = gpar(fontsize=8),
  row_title_gp = gpar(fontsize=8),
  column_names_gp = gpar(fontsize=8),
  column_title_gp = gpar(fontsize=8),
  column_dend_gp = gpar(lwd=0.5),
  row_dend_gp = gpar(lwd=0.5),
  column_names_rot = 45,
  border = T,
  show_column_dend = F,
  cell_fun = function(j, i, x, y, w, h, col) {
    grid.text(
      overlap.pvalue.mat.draw[i, j],
      x, y,
      gp = gpar(col = "black", fontsize = 7)
    )
  },
  heatmap_legend_param = list(
    title = "-log10(P-value)",
    legend_direction = "horizontal",
    title_gp = gpar(fontsize = 6),
    labels_gp = gpar(fontsize = 6)
  ),
  use_raster = F,
  column_title = "Binding & regulation (ZNF398)",
  top_annotation = column_ha
)

outfile  <- paste0(path_results, "/overlap.dna.bound.pvalue.hm.pdf")
pdf(file = outfile, paper = "a4", w = unit(2.4,'cm'), h = unit(2,'cm'))
draw(overlap.pvalue.hm, heatmap_legend_side="bottom")
dev.off()

# Intersection: bound & regulated ----
bound_and_reg_list <- list(
  "bound_down_2d" = gene_overlap_mat@go.nested.list$down_reg_2d$dna_bound@intersection,
  "bound_down_3d" = gene_overlap_mat@go.nested.list$down_reg_3d$dna_bound@intersection,
  "bound_up_2d" = gene_overlap_mat@go.nested.list$up_reg_2d$dna_bound@intersection,
  "bound_up_3d" = gene_overlap_mat@go.nested.list$up_reg_3d$dna_bound@intersection
)
enrichr_bound_and_reg <- list()
for(i in names(bound_and_reg_list)) {
  enrichr_bound_and_reg[[i]] <- enrichR::enrichr(
    bound_and_reg_list[[i]],
    databases = c(
      "GO_Molecular_Function_2023",
      "GO_Cellular_Component_2023",
      "GO_Biological_Process_2023",
      "MSigDB_Hallmark_2020"
    )
  )
  outfile <- paste0(path_results,"/enrichr.",i,".xlsx")
  saveXLSresEdgeR2(enrichr_bound_and_reg[[i]], outfile = outfile)
}