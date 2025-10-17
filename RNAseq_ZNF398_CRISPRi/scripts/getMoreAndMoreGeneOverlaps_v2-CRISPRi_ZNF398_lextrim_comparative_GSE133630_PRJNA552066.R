# setup ------
setwd("/home/epigen/Martello_3D_human_epiblast/RNAseq/dataset/v2-CRISPRi_ZNF398_lextrim")

library(RNAseqRtools)
library(GeneOverlap)
library(enrichR)
library(ggpubr)

path_binding_info <- "/home/epigen/Martello_3D_human_epiblast/ChIPseq/bit_encode-chip-seq-pipeline/dataset/comparative_GSE133630_PRJNA552066/peak/ChIP-seq_of_ZNF398_in_BG01V_human_embryonic_stem_cells-no_control.ChIP-seq_of_ZNF398_in_H9_human_embryonic_stem_cells-no_control.idr.bchrfilt.binding_info.gz"
path_dge_toptable <- "/home/epigen/Martello_3D_human_epiblast/RNAseq/dataset/v2-CRISPRi_ZNF398_lextrim/DGE_sep_group/edger.toptable_clean.ALL_contrast.mark_seqc.header_added.gz"
path_dge_oe <- "/home/epigen/Martello_3D_human_epiblast/public/figures/data/RNA-seq_bulk POLO/DE analysis POLO ZNF398.xlsx"
path_dge_toptable_oe <- "/home/epigen/Martello_3D_human_epiblast/public/figures/data/RNA-seq_bulk POLO/results/edger.de_contrasts.table.xlsx"

path_normcounts <- "results/normcounts_tmm.logcpm.txt.gz"
path_metadata <- "metadata.txt"

path_results <- "results/more_and_more_gene_overlaps"
if(!dir.exists(path_results)) dir.create(path_results)

# read data ------
binding_info <- read.delim(path_binding_info)
dge_toptable <- read.delim(path_dge_toptable)
normcounts <- read.delim(path_normcounts)
metadata <- read.delim(path_metadata)

# DESeq results
deseq_results <- openxlsx::read.xlsx(
  path_dge_oe,
  sheet = 1
)
contrast_list <- colnames(deseq_results)[c(3,5,7,9,11)]
deseq_results_l <- lapply(
  contrast_list,
  function(i) {
    pos_i <- which(colnames(deseq_results) == i)
    cidx <- c(pos_i,pos_i+1)
    res <- deseq_results[-1,c(1,2,cidx)]
    colnames(res) <- deseq_results[1,c(1,2,cidx)]
    res$log2FoldChange <- as.numeric(res$log2FoldChange)
    res$padj <- as.numeric(res$padj)
    res <- res[order(res$padj),]
    return(res)
  }
)
names(deseq_results_l) <- gsub("\\.","\\_",contrast_list)

dge_oe <- deseq_results_l$ZNF398_3D_SB_VS_EMPTY_3D_SB
dge_toptable_oe <- openxlsx::read.xlsx(path_dge_toptable_oe, sheet = '3D_ZNF398_SB43_vs_3D_Empty_SB43')

dge_toptable <- dge_toptable[dge_toptable$contrast!='ZNF398_NI_2D_Day7_vs_P_NI_2D_Day7',]

# generate gene lists ----
binding_info <- binding_info[-which(is.na(binding_info$gene)),]
binding_info <- binding_info[binding_info$gene!="ZNF398",]
# binding_info <- binding_info[which(abs(binding_info$distTSS)<=1e3),]
gene_lists <- list(
  "dna_bound"=unique(binding_info$gene),
  "promoter_bound"=unique(binding_info$gene[grepl("PLS",binding_info$ccRE)]),
  "enhancer_bound"=unique(binding_info$gene[grepl("ELS",binding_info$ccRE)]),
  "up_reg_2d"=dge_toptable[dge_toptable$contrast=="ZNF398_I_2D_Day7_vs_ZNF398_NI_2D_Day7" & dge_toptable$significance==1,'GeneID'],
  "down_reg_2d"=dge_toptable[dge_toptable$contrast=="ZNF398_I_2D_Day7_vs_ZNF398_NI_2D_Day7" & dge_toptable$significance==-1,'GeneID'],
  "up_reg_3d"=dge_toptable[dge_toptable$contrast=="ZNF398_I_3D_Day4_vs_ZNF398_NI_3D_Day4" & dge_toptable$significance==1,'GeneID'],
  "down_reg_3d"=dge_toptable[dge_toptable$contrast=="ZNF398_I_3D_Day4_vs_ZNF398_NI_3D_Day4" & dge_toptable$significance==-1,'GeneID'],
  # "up_reg_3d_oe"=dge_oe[which(dge_oe$log2FoldChange>=.5 & dge_oe$padj<0.05), "Symbols"],
  # "down_reg_3d_oe"=dge_oe[which(dge_oe$log2FoldChange<=(-.5) & dge_oe$padj<0.05), "Symbols"]
  "up_reg_3d_oe"=dge_toptable_oe[which(dge_toptable_oe$logFC>=.5 & dge_toptable_oe$FDR<0.05), "gene_name"],
  "down_reg_3d_oe"=dge_toptable_oe[which(dge_toptable_oe$logFC<=(-.5) & dge_toptable_oe$FDR<0.05), "gene_name"]
  )

gene_lists[["up_reg_both"]] <- intersect(gene_lists[["up_reg_2d"]], gene_lists[["up_reg_3d"]])
gene_lists[["down_reg_both"]] <- intersect(gene_lists[["down_reg_2d"]], gene_lists[["down_reg_3d"]])
gene_lists_len <- unlist(lapply(gene_lists, length))

# gene overlap ----
gene_overlap_mat <- GeneOverlap::newGOM(gsetA = gene_lists[1:3],gene_lists[-c(1:3)], spec = "hg19.gene")

# scatter plot ----
dge_merged <- merge(
  dge_toptable,
  # dge_oe,
  dge_toptable_oe,
  by.x = "GeneID",
  by.y = "gene_name",
  all = F
)
# dge_merged$log2FoldChange[which(is.na(dge_merged$log2FoldChange))] <- 0
# dge_merged$padj[which(is.na(dge_merged$padj))] <- 1
# dge_merged$logFC[which(is.na(dge_merged$logFC))] <- 0
# dge_merged$Pvalue_adj[which(is.na(dge_merged$Pvalue_adj))] <- 1

# dge_merged$score_ko <- -log10(dge_merged$Pvalue_adj) * dge_merged$logFC
# dge_merged$score_oe <- -log10(dge_merged$padj) * dge_merged$log2FoldChange

dge_merged$score_ko <- -log10(dge_merged$Pvalue_adj) * dge_merged$logFC.x
dge_merged$score_oe <- -log10(dge_merged$FDR) * dge_merged$logFC.y

for(i in names(gene_lists)) {
  dge_merged[[i]] <- dge_merged$GeneID %in% gene_lists[[i]]  
}
dge_merged$status <- 'none'
dge_merged$status[which(dge_merged$up_reg_both & dge_merged$down_reg_3d_oe & dge_merged$dna_bound)] <- 'bound & repressed'
dge_merged$status[which(dge_merged$down_reg_both & dge_merged$up_reg_3d_oe & dge_merged$dna_bound)] <- 'bound & activated'

dge_merged$status_expr <- 'none'
dge_merged$status_expr[which(dge_merged$up_reg_both & dge_merged$down_reg_3d_oe)] <- 'repressed'
dge_merged$status_expr[which(dge_merged$down_reg_both & dge_merged$up_reg_3d_oe)] <- 'activated'

status_lab_nm <- table(dge_merged[dge_merged$contrast==dge_merged$contrast[1],'status'])
status_lab <- paste0(names(status_lab_nm), " [n=",status_lab_nm,"]")
status_lab[3] <- "none"
names(status_lab) <- names(status_lab_nm)
dge_merged$status <- status_lab[dge_merged$status]
pal_status <- c("#54278F", "lightgrey", "#D94801")
names(pal_status) <- unique(dge_merged$status)[c(3,1,2)]
dge_merged$contrast_lab <- dge_merged$contrast
dge_merged$contrast_lab <- gsub(".*(2D|3D).*","\\1",dge_merged$contrast_lab)

dge_scatter_plt <- ggplot(
  # mapping = aes(x=log2FoldChange, y=logFC, col=status)
  mapping = aes(x=logFC.y, y=logFC.x, col=status)
) + geom_point(data=dge_merged[which(dge_merged$status=='none'),], size=0.4) +
  # geom_point(data=dge_merged[which(dge_merged$status_expr=='repressed'),], size=0.4, col='#CCCCFF', alpha=0.6) +
  # geom_point(data=dge_merged[which(dge_merged$status_expr=='activated'),], size=0.4, col='#FFCCCC', alpha=0.6) +
  # ggpubr::stat_cor(data=dge_merged, aes(x=log2FoldChange, y=logFC), inherit.aes = F) +
  geom_point(data=dge_merged[which(dge_merged$status!='none'),], size=0.4) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  facet_grid(~contrast_lab) +
  scale_color_manual(values = pal_status) +
  theme_bw() + my_theme +
  guides(col=guide_legend(override.aes = list(size=3))) +
  theme(legend.position = 'right', aspect.ratio = 1) +
  ylab("logFC KO") + xlab("logFC OE")
dge_scatter_plt <- ggrastr::rasterise(dge_scatter_plt, dpi=600)

outfile <- paste0(path_results,"/dge_scatter_plt.pdf")
message(" -- saving to: ", outfile)
pdf(file = outfile, paper = "a4r", width = unit(5.5,'cm'), height = unit(2.2,'cm'))
print(dge_scatter_plt)
dev.off()

# dge_merged_sel <- unique(dge_merged[,c(1,12:23)])
dge_merged_sel <- unique(dge_merged[,c(1,15:26)])
dge_merged_sel <- dge_merged_sel[dge_merged_sel$status!="none",]
dge_merged_sel <- dge_merged_sel[order(dge_merged_sel$status),]

outfile <- paste0(path_results, "/dge_merged_sel.xlsx")
RNAseqRtools::saveXLSresEdgeR2(dge_merged_sel, outfile, name="Bound and regulated")

# file for motifs ----
binding_info_motifs <- unique(binding_info[,c("idx","gene")])
binding_info_motifs <- binding_info_motifs[binding_info_motifs$gene%in%c(gene_lists$up_reg_3d, gene_lists$down_reg_3d),]
binding_info_motifs$cluster <- NA
binding_info_motifs$cluster[binding_info_motifs$gene %in% gene_lists$up_reg_3d] <- 'repressed'
binding_info_motifs$cluster[binding_info_motifs$gene %in% gene_lists$down_reg_3d] <- 'activated'
binding_info_motifs <- binding_info_motifs[,-2]
colnames(binding_info_motifs)[1] <- 'loc'
outfile <- paste0(path_results, "/binding_info_motifs.txt")

write.table(
  binding_info_motifs,
  file = outfile,
  row.names = F,
  col.names = T,
  sep = "\t",
  quote = F
)

# heatmap fc ----
gene.idx <- dge_merged_sel$GeneID
# exclude if no binding in prom/enhancer
gene.idx <- gene.idx[gene.idx %in% unique(c(gene_lists$promoter_bound,gene_lists$enhancer_bound))]

fc.mat <- lapply(
  split(dge_toptable, dge_toptable$contrast),
  function(x) {
    x <- x[,c("GeneID", "logFC", "Pvalue_adj")]
    colnames(x)[3] <- "FDR"
    x[x$GeneID %in% gene.idx, ]
  }
)
names(fc.mat) <- gsub("ZNF398_I_2D_Day7_vs_ZNF398_NI_2D_Day7", "KO_2D", names(fc.mat))
names(fc.mat) <- gsub("ZNF398_I_3D_Day4_vs_ZNF398_NI_3D_Day4", "KO_3D", names(fc.mat))
fc.mat[["OE_3D"]] <- dge_toptable_oe[,c("gene_name", "logFC", "FDR")]
colnames(fc.mat[["OE_3D"]])[1] <- "GeneID"
fc.mat[["OE_3D"]] <- fc.mat[["OE_3D"]][fc.mat[["OE_3D"]]$GeneID %in% gene.idx, ]

for(i in names(fc.mat)) {
  colnames(fc.mat[[i]])[-1] <- paste0(colnames(fc.mat[[i]])[-1], "_", i)
}
fc.mat <- Reduce(function(a, b) merge(a, b, all=T, by="GeneID"), fc.mat)
rownames(fc.mat) <- fc.mat[,1]
fc.mat <- fc.mat[,-1]
fc.mat <- fc.mat[,c(1,3,5)]
colnames(fc.mat) <- gsub("logFC_", "", colnames(fc.mat))
fc.mat <- as.matrix(fc.mat)

color_palette <- circlize::colorRamp2(
  c(-1.5,0,1.5),
  colors = c("#54278F","white","#D94801")
)

expr_fc_hm <- Heatmap(
  fc.mat,
  col=color_palette,
  show_row_dend = F,
  show_row_names = T,
  cluster_columns = F,
  border = T,
  split = 2,
  row_title = c("Activated", "Repressed"),
  row_title_gp = gpar(fontsize=8),
  column_split = c("KO", "KO", "OE"),
  # column_title = "DNA binding",
  column_names_gp = gpar(fontsize=8),
  column_title_gp = gpar(fontsize=8),
  row_names_gp = gpar(fontsize=8),
  show_heatmap_legend = T,
  width = 4,
  heatmap_legend_param = list(
    title = "logFC",
    title_gp = gpar(fontsize=8, fontface="plain"),
    labels_gp = gpar(fontsize=8, fontface="plain"),
    values_gp = gpar(fontsize=8, fontface="plain")
  )
)

# binding ---
binding_mat <- matrix(
  0,
  nrow = nrow(fc.mat),
  ncol = 3,
  dimnames = list(rownames(fc.mat), c("PLS","pELS","dELS"))
)
bd_map <- unique(binding_info[,c("gene","ccRE")])
bd_map <- ddply(bd_map, .(gene), summarize, ccRE=paste0(ccRE, collapse=","))
bd_map <- bd_map[match(rownames(fc.mat), bd_map$gene),]
print(nrow(bd_map) == nrow(binding_mat))
for(j in colnames(binding_mat)) {
  binding_mat[,j] <- grepl(j, bd_map$ccRE)
}
binding_hm <- Heatmap(
  binding_mat,
  col=c('white','lightblue'),
  show_row_dend = F,
  show_row_names = T,
  cluster_columns = F,
  border = T,
  column_title = "DNA binding",
  column_names_gp = gpar(fontsize=8),
  column_title_gp = gpar(fontsize=8),
  row_names_gp = gpar(fontsize=8),
  show_heatmap_legend = F,
  width = 3,
)
hm_bd <- binding_hm + expr_fc_hm
draw(hm_bd, main_heatmap=2)

outfile <- paste0(path_results,"/hm_bd.fc.pdf")
message('-- saving to: ',outfile)
pdf(file = outfile, paper = "a4", w = unit(2.8,'cm'), h = unit(4.8,'cm'))
draw(hm_bd, main_heatmap=2, heatmap_legend_side="bottom")
dev.off()
