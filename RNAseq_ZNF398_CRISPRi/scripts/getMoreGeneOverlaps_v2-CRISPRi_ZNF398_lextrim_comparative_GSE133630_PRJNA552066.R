# setup ------
setwd("/home/epigen/Martello_3D_human_epiblast/RNAseq/dataset/v2-CRISPRi_ZNF398_lextrim")

library(RNAseqRtools)
library(GeneOverlap)
library(enrichR)
library(ggpubr)

path_binding_info <- "/home/epigen/Martello_3D_human_epiblast/ChIPseq/bit_encode-chip-seq-pipeline/dataset/comparative_GSE133630_PRJNA552066/peak/ChIP-seq_of_ZNF398_in_BG01V_human_embryonic_stem_cells-no_control.ChIP-seq_of_ZNF398_in_H9_human_embryonic_stem_cells-no_control.idr.bchrfilt.binding_info.gz"
path_dge_toptable <- "/home/epigen/Martello_3D_human_epiblast/RNAseq/dataset/v2-CRISPRi_ZNF398_lextrim/DGE_sep_group/edger.toptable_clean.ALL_contrast.mark_seqc.header_added.gz"
path_normcounts <- "results/normcounts_tmm.logcpm.txt.gz"
path_metadata <- "metadata.txt"

path_results <- "results/more_gene_overlaps"
if(!dir.exists(path_results)) dir.create(path_results)

# read data ------
binding_info <- read.delim(path_binding_info)
dge_toptable <- read.delim(path_dge_toptable)
normcounts <- read.delim(path_normcounts)
metadata <- read.delim(path_metadata)

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
  "down_reg_3d"=dge_toptable[dge_toptable$contrast=="ZNF398_I_3D_Day4_vs_ZNF398_NI_3D_Day4" & dge_toptable$significance==-1,'GeneID']
  )

gene_lists[["up_reg_both"]] <- intersect(gene_lists[["up_reg_2d"]], gene_lists[["up_reg_3d"]])
gene_lists[["down_reg_both"]] <- intersect(gene_lists[["down_reg_2d"]], gene_lists[["down_reg_3d"]])
gene_lists_len <- unlist(lapply(gene_lists, length))
# gene overlap ----
gene_overlap_mat <- GeneOverlap::newGOM(gsetA = gene_lists[1:3],gene_lists[-c(1:3)], spec = "hg19.gene")

# Visualize gene overlap: bar/dot plot ----
overlap.pvalue.df <- as.data.frame(getMatrix(gene_overlap_mat,'pval'))
overlap.pvalue.df$bound_type <- rownames(overlap.pvalue.df)
overlap.pvalue.df$var_type <- 'pval'

overlap.intersect.df <- as.data.frame(getMatrix(gene_overlap_mat,'intersection'))
overlap.int.frac.df <- as.data.frame(t(apply(overlap.intersect.df,1,function(x) x/gene_lists_len[colnames(overlap.intersect.df)])))
overlap.or.df <- as.data.frame(getMatrix(gene_overlap_mat,'odds.ratio'))

overlap.intersect.df$bound_type <- rownames(overlap.intersect.df)
overlap.intersect.df$var_type <- 'intersect'

overlap.int.frac.df$bound_type <- rownames(overlap.int.frac.df)
overlap.int.frac.df$var_type <- 'int_frac'

overlap.or.df$bound_type <- rownames(overlap.or.df)
overlap.or.df$var_type <- 'odds_ratio'

overlap.df <- rbind(
  overlap.pvalue.df,
  overlap.intersect.df,
  overlap.int.frac.df,
  overlap.or.df
)
rownames(overlap.df) <- NULL

outfile  <- paste0(path_results, "/overlap.df.xlsx")
RNAseqRtools::saveXLSresEdgeR2(
  overlap.df,
  outfile = outfile,
  name = "Bound and regulated"
  )

# plt1: barplot frac bound dge ---
overlap.df <- reshape2::melt(overlap.df, id.vars=c("bound_type","var_type"))
overlap.df$direction <- gsub("^(.*)_reg_(.*)$","\\1",overlap.df$variable)
overlap.df$cond <- gsub("^(.*)_reg_(.*)$","\\2",overlap.df$variable)
overlap.df$cond <- factor(overlap.df$cond, levels=c("both","2d","3d"))

pal.dge.summary <- c("#D94801","#D94801","#D94801",
                     "#54278F", "#54278F", "#54278F")
names(pal.dge.summary) <- c("up_reg_both","up_reg_2d","up_reg_3d",
                            "down_reg_both","down_reg_2d","down_reg_3d")

overlap.int.up.plt <- ggplot(
  overlap.df[overlap.df$var_type=="int_frac" & overlap.df$bound_type=="dna_bound" & overlap.df$direction=="up",],
  aes(x=value, y=cond, fill=variable)
) +
  # facet_wrap(~direction, ncol=1) +
  geom_col(aes(col=cond), fill='lightblue',show.legend = F) +
  ylab("Regulated in") + xlab("% of bound genes") +
  scale_fill_manual(values = pal.dge.summary) +
  scale_color_manual(values = c("2d" = 'white', "3d" = 'white', "both" = 'black')) +
  scale_alpha_manual(values = c("2d" = 0.6, "3d" = 1, "both" = 0.3)) +
  theme_classic() + my_theme +
  theme(panel.grid.major = element_blank())

overlap.int.down.plt <- ggplot(
  overlap.df[overlap.df$var_type=="int_frac" & overlap.df$bound_type=="dna_bound" & overlap.df$direction=="down",],
  aes(x=value, y=cond, fill=variable)
) +
  # facet_wrap(~direction, ncol=1) +
  geom_col(aes(col=cond), fill='lightblue',show.legend = F) +
  ylab("Regulated in") + xlab("% of bound genes") +
  scale_fill_manual(values = pal.dge.summary) +
  scale_color_manual(values = c("2d" = 'white', "3d" = 'white', "both" = 'black')) +
  scale_alpha_manual(values = c("2d" = 0.6, "3d" = 1, "both" = 0.3)) +
  theme_classic() + my_theme +
  theme(panel.grid.major = element_blank())

# plt2: dot plot bound dge overlaps ---
to.unmelt <- overlap.df[overlap.df$direction=="up",]
to.unmelt$bound_type <- gsub("dna_bound","all",to.unmelt$bound_type)
to.unmelt$bound_type <- gsub("_bound","",to.unmelt$bound_type)
to.unmelt$bound_type <- gsub("(\\w)(\\w*)", "\\U\\1\\L\\2", to.unmelt$bound_type, perl=TRUE)
to.unmelt$bound_type <- factor(to.unmelt$bound_type, levels=c("All","Promoter","Enhancer"))
overlap.unm.up <- dcast(to.unmelt, variable+direction+cond+bound_type~var_type, value.var = 'value', fun.aggregate = mean)
overlap.unm.up$pval[-log10(overlap.unm.up$pval) > 15] <- 1e-15
overlap.unm.up.plt <- ggplot(
  overlap.unm.up,
  aes(x=bound_type, y=cond, fill=-log10(pval), size=intersect)
) +
  # facet_wrap(~direction, ncol=1) +
  scale_size(range = c(2, 7)) +
  geom_point(shape=21) +
  xlab("Binding region") + ylab(NULL) +
  scale_fill_gradient(
    low = RColorBrewer::brewer.pal(9,"Reds")[2],
    high = "#D94801",
    limits = c(1, 15),
    na.value = 'white'
  ) +
  guides(
    fill = guide_colorbar(barwidth = 4.5, barheight = .6, order = 2),
    size = guide_legend(order = 1)
  ) +
  theme_bw() + my_theme

to.unmelt <- overlap.df[overlap.df$direction=="down",]
to.unmelt$bound_type <- gsub("dna_bound","all",to.unmelt$bound_type)
to.unmelt$bound_type <- gsub("_bound","",to.unmelt$bound_type)
to.unmelt$bound_type <- gsub("(\\w)(\\w*)", "\\U\\1\\L\\2", to.unmelt$bound_type, perl=TRUE)
to.unmelt$bound_type <- factor(to.unmelt$bound_type, levels=c("All","Promoter","Enhancer"))
overlap.unm.down <- dcast(to.unmelt, variable+direction+cond+bound_type~var_type, value.var = 'value', fun.aggregate = mean)

overlap.unm.down.plt <- ggplot(
  overlap.unm.down,
  aes(x=bound_type, y=cond, fill=-log10(pval), size=intersect)
) +
  # facet_wrap(~direction, ncol=1) +
  scale_size(range = c(2, 7)) +
  geom_point(shape=21) +
  xlab("Binding region") + ylab(NULL) +
  scale_fill_gradient(
    low = RColorBrewer::brewer.pal(9,"Purples")[2],
    high = "#54278F",
    limits = c(1,15),
    na.value = 'white'
  ) +
  guides(
    fill = guide_colorbar(barwidth = 4.5, barheight = .6, order = 2),
    size = guide_legend(order = 1)
  ) +
  theme_bw() + my_theme

# arranged plt ---
overlap.plt.up <- ggpubr::ggarrange(overlap.int.up.plt,overlap.unm.up.plt,
                  align="hv", common.legend = T, widths = c(1,1.4))
overlap.plt.up <- ggpubr::annotate_figure(overlap.plt.up,
                                          top = text_grob('ZNF398 DNA binding & up-regulated genes (overlap)', size = 8))
outfile  <- paste0(path_results, "/overlap.plt.up.pdf")
pdf(file = outfile, paper = "a4", w = unit(4,'cm'), h = unit(1.8,'cm'))
print(overlap.plt.up)
dev.off()

overlap.plt.down <- ggpubr::ggarrange(overlap.int.down.plt,overlap.unm.down.plt,
                            align="hv", common.legend = T, widths = c(1,1.4))
overlap.plt.down <- ggpubr::annotate_figure(overlap.plt.down,
                                            top = text_grob('ZNF398 DNA binding & down-regulated genes (overlap)', size = 8))

outfile  <- paste0(path_results, "/overlap.plt.down.pdf")
pdf(file = outfile, paper = "a4", w = unit(4,'cm'), h = unit(1.8,'cm'))
print(overlap.plt.down)
dev.off()

# DGE counts ----
dge.down <- melt(gene_lists_len[grep('down',names(gene_lists_len))])
dge.down$cond_full <- rownames(dge.down)
dge.down$cond <- gsub(".*_","",rownames(dge.down))
dge.down$cond <- factor(dge.down$cond, levels=c("both","2d","3d"))

dge.down.plt <- ggplot(dge.down, aes(y=cond, x=value, fill=cond_full)) +
  geom_col(aes(col=cond), show.legend = F,fill='#BCBDDC') +
  ylab("Regulated in") + xlab("N. regulated genes") +
  scale_color_manual(values = c("2d" = 'white', "3d" = 'white', "both" = 'black')) +
  # scale_alpha_manual(values = c("2d" = 0.6, "3d" = 1, "both" = 0.3)) +
  scale_fill_manual(values = pal.dge.summary) +
  theme_classic() + my_theme +
  theme(panel.grid.major = element_blank())

overlap.plt.down2 <- ggpubr::ggarrange(dge.down.plt, overlap.int.down.plt+ylab(NULL),overlap.unm.down.plt,
                                      align="hv", common.legend = T, widths = c(1,1,1.5),nrow = 1)
overlap.plt.down2 <- ggpubr::annotate_figure(overlap.plt.down2,
                                            top = text_grob('ZNF398 DNA binding & down-regulated genes (overlap)', size = 8))

outfile  <- paste0(path_results, "/overlap.plt.down.2.pdf")
pdf(file = outfile, paper = "a4", w = unit(5.2,'cm'), h = unit(1.8,'cm'))
print(overlap.plt.down2)
dev.off()

# DGE and bound counts ----
dge.bound.down.plt <- ggplot(overlap.unm.down[overlap.unm.down$bound_type=="All",], aes(y=cond, x=intersect, fill=variable)) +
  geom_col(aes(col=cond), show.legend = F,fill='lightblue') +
  ylab("Regulated & bound in") + xlab("N. of genes") +
  scale_color_manual(values = c("2d" = 'white', "3d" = 'white', "both" = 'black')) +
  # scale_alpha_manual(values = c("2d" = 0.6, "3d" = 1, "both" = 0.3)) +
  scale_fill_manual(values = pal.dge.summary) +
  theme_classic() + my_theme +
  theme(panel.grid.major = element_blank())

overlap.unm.or.down.plt <- ggplot(
  overlap.unm.down,
  aes(x=bound_type, y=cond, fill=-log10(pval), size=odds_ratio)
) +
  # facet_wrap(~direction, ncol=1) +
  scale_size(range = c(2, 7)) +
  geom_point(shape=21) +
  xlab("Binding region") + ylab(NULL) +
  scale_fill_gradient(
    low = RColorBrewer::brewer.pal(9,"Purples")[2],
    high = "#54278F",
    limits = c(1,15),
    na.value = 'white'
  ) +
  guides(
    fill = guide_colorbar(barwidth = 4.5, barheight = .6, order = 2),
    size = guide_legend(order = 1)
  ) +
  theme_bw() + my_theme +
  theme(axis.text.y = element_blank())

overlap.plt.down3 <- ggpubr::ggarrange(dge.bound.down.plt, overlap.unm.or.down.plt+ylab(NULL),
                                       align="hv", common.legend = T, widths = c(1,1.4),nrow = 1)
overlap.plt.down3 <- ggpubr::annotate_figure(overlap.plt.down3,
                                             top = text_grob('ZNF398 DNA binding & down-regulated genes (overlap)', size = 8))

outfile  <- paste0(path_results, "/overlap.plt.down.3.pdf")
pdf(file = outfile, paper = "a4", w = unit(4.1,'cm'), h = unit(1.9,'cm'))
print(overlap.plt.down3)
dev.off()

dge.bound.up.plt <- ggplot(overlap.unm.up[overlap.unm.up$bound_type=="All",], aes(y=cond, x=intersect, fill=variable)) +
  geom_col(aes(col=cond), show.legend = F,fill='lightblue') +
  ylab("Regulated & bound in") + xlab("N. of genes") +
  scale_color_manual(values = c("2d" = 'white', "3d" = 'white', "both" = 'black')) +
  # scale_alpha_manual(values = c("2d" = 0.6, "3d" = 1, "both" = 0.3)) +
  scale_fill_manual(values = pal.dge.summary) +
  theme_classic() + my_theme +
  theme(panel.grid.major = element_blank())

overlap.unm.or.up.plt <- ggplot(
  overlap.unm.up,
  aes(x=bound_type, y=cond, fill=-log10(pval), size=odds_ratio)
) +
  # facet_wrap(~direction, ncol=1) +
  scale_size(range = c(2, 7)) +
  geom_point(shape=21) +
  xlab("Binding region") + ylab(NULL) +
  scale_fill_gradient(
    low = RColorBrewer::brewer.pal(9,"Reds")[2],
    high = "#D94801",
    limits = c(1,15),
    na.value = 'white'
  ) +
  guides(
    fill = guide_colorbar(barwidth = 4.5, barheight = .6, order = 2),
    size = guide_legend(order = 1)
  ) +
  theme_bw() + my_theme +
  theme(axis.text.y = element_blank())

overlap.plt.up3 <- ggpubr::ggarrange(dge.bound.up.plt, overlap.unm.or.up.plt+ylab(NULL),
                                       align="hv", common.legend = T, widths = c(1,1.4),nrow = 1)
overlap.plt.up3 <- ggpubr::annotate_figure(overlap.plt.up3,
                                             top = text_grob('ZNF398 DNA binding & up-regulated genes (overlap)', size = 8))

outfile  <- paste0(path_results, "/overlap.plt.up.3.pdf")
pdf(file = outfile, paper = "a4", w = unit(4.1,'cm'), h = unit(1.9,'cm'))
print(overlap.plt.up3)
dev.off()

# up and down
to.unmelt <- overlap.df[overlap.df$variable %in% c("up_reg_both", "down_reg_both"), ]
to.unmelt$bound_type <- gsub("dna_bound","all",to.unmelt$bound_type)
to.unmelt$bound_type <- gsub("_bound","",to.unmelt$bound_type)
to.unmelt$bound_type <- gsub("(\\w)(\\w*)", "\\U\\1\\L\\2", to.unmelt$bound_type, perl=TRUE)
to.unmelt$bound_type <- factor(to.unmelt$bound_type, levels=c("All","Promoter","Enhancer"))
overlap.unm <- dcast(to.unmelt, variable+direction+cond+bound_type~var_type, value.var = 'value', fun.aggregate = mean)
overlap.unm$sign <- -1
overlap.unm$sign[which(overlap.unm$direction=="up")] <- 1

dge.bound.plt <- ggplot(overlap.unm[overlap.unm$bound_type=="All",], aes(y=variable, x=intersect, fill=variable)) +
  geom_col(aes(col=cond), show.legend = F,fill='lightblue', width = 0.8) +
  ylab("Regulated & bound in") + xlab("N. of genes") +
  scale_color_manual(values = c("2d" = 'white', "3d" = 'white', "both" = 'black')) +
  # scale_alpha_manual(values = c("2d" = 0.6, "3d" = 1, "both" = 0.3)) +
  scale_fill_manual(values = pal.dge.summary) +
  theme_classic() + my_theme +
  theme(panel.grid.major = element_blank())

overlap.unm.or.plt <- ggplot(
  overlap.unm,
  aes(x=bound_type, y=variable, fill=-log10(pval) * sign, size=odds_ratio)
) +
  # facet_wrap(~direction, ncol=1) +
  scale_size(range = c(2, 7)) +
  geom_point(shape=21) +
  xlab("Binding region") + ylab(NULL) +
  scale_fill_gradient2(
    low = "#54278F",
    mid = "white",
    high = "#D94801",
    # limits = c(1,15),
    na.value = 'white'
  ) +
  guides(
    fill = guide_colorbar(barwidth = 4.5, barheight = .6, order = 2),
    size = guide_legend(order = 1)
  ) +
  theme_bw() + my_theme +
  theme(axis.text.y = element_blank())


overlap.plt.both <- ggpubr::ggarrange(dge.bound.plt, overlap.unm.or.plt+ylab(NULL),
                                     align="hv", common.legend = T, widths = c(1, 1.4), nrow = 1)
overlap.plt.both <- ggpubr::annotate_figure(overlap.plt.both,
                                           top = text_grob('ZNF398 DNA binding & regulated genes (overlap)', size = 8))

outfile  <- paste0(path_results, "/overlap.plt.both.pdf")
pdf(file = outfile, paper = "a4r", w = unit(5,'cm'), h = unit(1.8,'cm'))
print(overlap.plt.both)
dev.off()
# EnrichR intersection: bound & regulated ----
bound_and_reg_list <- list(
  "bound_down_2d" = gene_overlap_mat@go.nested.list$down_reg_2d$dna_bound@intersection,
  "bound_down_3d" = gene_overlap_mat@go.nested.list$down_reg_3d$dna_bound@intersection,
  "bound_down_both" = gene_overlap_mat@go.nested.list$down_reg_both$dna_bound@intersection,
  "bound_up_2d" = gene_overlap_mat@go.nested.list$up_reg_2d$dna_bound@intersection,
  "bound_up_3d" = gene_overlap_mat@go.nested.list$up_reg_3d$dna_bound@intersection,
  "bound_up_both" = gene_overlap_mat@go.nested.list$up_reg_both$dna_bound@intersection
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

# heatmap down ----
pal_analysis <- c(
  "P_NI_2D_Day7" = "#FFCCFF",
  "ZNF398_NI_2D_Day7" = "#C6DBEF",
  "ZNF398_I_2D_Day7" = "#D94801" ,
  "ZNF398_NI_3D_Day4" = "#2171B5",
  "ZNF398_I_3D_Day4" = "#7F2704" 
)

# both ---
sample.idx <- gsub("\\-","\\.",metadata$sample[metadata$condition_sepgr!="P_NI_2D_Day7"])
gene.idx <- bound_and_reg_list$bound_down_both
# exclude if no binding in prom/enhancer
gene.idx <- gene.idx[gene.idx %in% unique(c(gene_lists$promoter_bound,gene_lists$enhancer_bound))]

expr_mat <- normcounts[normcounts$gene_name %in% gene.idx,]
rownames(expr_mat) <- expr_mat$gene_name
expr_mat <- expr_mat[,-c(1:3)]
expr_mat <- as.matrix(expr_mat)
expr_mat <- expr_mat[,sample.idx]

expr_hm <- get_heatmap4(
  m = expr_mat,
  show_row_dend = F,
  show_row_names = T,
  bm = F,
  myLegend = "logCPM",
  myPalette = c("#54278F","white","#D94801"),
  myZscale = c(-1.5,0,1.5),
  retHm = T,
  annotDF = data.frame("group"=metadata$condition_sepgr[match(colnames(expr_mat),gsub("\\-","\\.",metadata$sample))]),
  annotCol = list("group"=pal_analysis),
  annotation_width_heigth = .7,
  show_annotation_name = F,
  row_title_gp = gpar(fontsize = 8),
  cluster_columns = F,
  column_title = "ZNF398 CRISPRi RNAseq",
)

# 2d ---
sample.idx <- gsub("\\-","\\.",metadata$sample[metadata$Culture=="2D" & metadata$condition_sepgr!="P_NI_2D_Day7"])
# gene.idx <- bound_and_reg_list$bound_down_both

expr_mat <- normcounts[normcounts$gene_name %in% gene.idx,]
rownames(expr_mat) <- expr_mat$gene_name
expr_mat <- expr_mat[,-c(1:3)]
expr_mat <- as.matrix(expr_mat)
expr_mat <- expr_mat[,sample.idx]

expr_hm_2d <- get_heatmap4(
  m = expr_mat,
  show_row_dend = F,
  show_row_names = T,
  bm = F,
  myLegend = "logCPM",
  myPalette = c("#54278F","white","#D94801"),
  myZscale = c(-1.5,0,1.5),
  retHm = T,
  annotDF = data.frame("group"=metadata$condition_sepgr[match(colnames(expr_mat),gsub("\\-","\\.",metadata$sample))]),
  annotCol = list("group"=pal_analysis),
  annotation_width_heigth = .7,
  show_annotation_name = F,
  row_title_gp = gpar(fontsize = 8),
  cluster_columns = F,
  column_title = "RNAseq (2D)",
)

# 3d ---
sample.idx <- gsub("\\-","\\.",metadata$sample[metadata$Culture=="3D" & metadata$condition_sepgr!="P_NI_2D_Day7"])
# gene.idx <- bound_and_reg_list$bound_down_both

expr_mat <- normcounts[normcounts$gene_name %in% gene.idx,]
rownames(expr_mat) <- expr_mat$gene_name
expr_mat <- expr_mat[,-c(1:3)]
expr_mat <- as.matrix(expr_mat)
expr_mat <- expr_mat[,sample.idx]

expr_hm_3d <- get_heatmap4(
  m = expr_mat,
  show_row_dend = F,
  show_row_names = T,
  bm = F,
  myLegend = "logCPM",
  myPalette = c("#54278F","white","#D94801"),
  myZscale = c(-1.5,0,1.5),
  retHm = T,
  annotDF = data.frame("group"=metadata$condition_sepgr[match(colnames(expr_mat),gsub("\\-","\\.",metadata$sample))]),
  annotCol = list("group"=pal_analysis),
  annotation_width_heigth = .7,
  show_annotation_name = F,
  row_title_gp = gpar(fontsize = 8),
  cluster_columns = F,
  column_title = "RNAseq (3D)",
)

# binding ---
binding_mat <- matrix(
  0,
  nrow = nrow(expr_mat),
  ncol = 3,
  dimnames = list(rownames(expr_mat), c("PLS","pELS","dELS"))
)
bd_map <- unique(binding_info[,c("gene","ccRE")])
bd_map <- ddply(bd_map, .(gene), summarize, ccRE=paste0(ccRE, collapse=","))
bd_map <- bd_map[match(rownames(expr_mat), bd_map$gene),]
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
)
hm_bd <- binding_hm + expr_hm

outfile  <- paste0(path_results, "/hm_bd.bound_down_both.pdf")
message('-- saving to: ',outfile)
pdf(file = outfile, paper = "a4", w = unit(5.5,'cm'), h = unit(5.5,'cm'))
draw(hm_bd, heatmap_legend_side="bottom")
dev.off()

hm_bd_sep <- binding_hm + expr_hm_2d + expr_hm_3d

outfile  <- paste0(path_results, "/hm_bd.bound_down_both.sep.pdf")
message('-- saving to: ',outfile)
pdf(file = outfile, paper = "a4", w = unit(5.5,'cm'), h = unit(5.5,'cm'))
draw(hm_bd_sep, heatmap_legend_side="bottom")
dev.off()
