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

path_results <- "results/cellmarker2"
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

# read cell markers ----
cell_marker <- openxlsx::read.xlsx("Cell_marker_Human.xlsx")
cell_marker <- cell_marker[cell_marker$cell_type=="Normal cell",]
cell_marker <- cell_marker[-which(is.na(cell_marker$Symbol)),]
cell_name_set <- lapply(split(cell_marker,cell_marker$cell_name),"[[","Symbol")

# gene overlap ----
cell_name_set_gom <- GeneOverlap::newGOM(gene_lists[4:7], cell_name_set, spec = "hg19.gene")
pval_mat <- GeneOverlap::getMatrix(cell_name_set_gom, name = 'pval')
pval_mat <- melt(pval_mat)
pval_mat$fdr <- p.adjust(pval_mat$value,method='fdr')
colnames(pval_mat)[3] <- 'pval'

or_mat <- GeneOverlap::getMatrix(cell_name_set_gom, name = 'odds.ratio')
or_mat <- melt(or_mat)
colnames(or_mat)[3] <- 'or'
emat <- merge(pval_mat, or_mat, by=c("Var1","Var2"))


pval_mat[grep("Neural|Neuron",pval_mat$Var2, ignore.case = T),]
pval_mat[grep("Epithelial",pval_mat$Var2, ignore.case = T),]
pval_mat[grep("Mesenchymal",pval_mat$Var2, ignore.case = T),]
pval_mat[grep("Stem cell",pval_mat$Var2, ignore.case = T),]

sel_name <- c(
  "Neural stem cell",
  "Neural progenitor cell",
  "Epithelial cell",
  "Epithelial progenitor cell",
  "Mesenchymal cell",
  "Pluripotent stem cell"
)
pval_mat_sel <- pval_mat[pval_mat$Var2%in%sel_name,]
emat_sel <- emat[emat$Var2%in%sel_name,]

int_nl <- GeneOverlap::getNestedList(cell_name_set_gom, name = 'intersection')
int_nl <- melt(int_nl[sel_name])

emat_sel$direction <- gsub("^(.*)_reg_(.*)$","\\1",emat_sel$Var1)
emat_sel$cond <- gsub("^(.*)_reg_(.*)$","\\2",emat_sel$Var1)
emat_sel <- emat_sel[emat_sel$or!=0,]

emat_sel_up <- emat_sel[emat_sel$direction=='up',]
emat_sel_up$fdr[emat_sel_up$fdr>0.05] <- 1
emat_up_plt <- ggplot(emat_sel_up, aes(x=cond,y=Var2, size=or, fill=-log10(fdr))) +
  geom_point(shape=21) +
  scale_size(range = c(2, 6), breaks = c(5, 15, 25), limits = c(0,53)) +
  xlab(NULL) + ylab(NULL) +
  scale_fill_gradient(
    low = RColorBrewer::brewer.pal(9,"Reds")[2],
    high = "#D94801",
    limits = c(1, 5),
    na.value = 'white'
  ) +
  guides(
    fill = guide_colorbar(barwidth = 3, barheight = .6, order = 2),
    size = guide_legend(order = 1, ncol = 1)
  ) +
  theme_bw() + my_theme + theme(legend.box = "vertical",
                                legend.box.just = 'left',
                                legend.justification = c(0, 0))

emat_sel_down <- emat_sel[emat_sel$direction=='down',]
emat_sel_down$fdr[emat_sel_down$fdr>0.05] <- 1
emat_down_plt <- ggplot(emat_sel_down, aes(x=cond,y=Var2, size=or, fill=-log10(fdr))) +
  geom_point(shape=21) +
  scale_size(range = c(2, 6), breaks = c(5, 15, 25), limits = c(0,53)) +
  xlab(NULL) + ylab(NULL) +
  scale_fill_gradient(
    low = RColorBrewer::brewer.pal(9,"Purples")[2],
    high = "#54278F",
    limits = c(1, 5),
    na.value = 'white'
  ) +
  guides(
    fill = guide_colorbar(barwidth = 3, barheight = .6, order = 2),
    size = guide_legend(order = 1, ncol = 1)
  ) +
  theme_bw() + my_theme + theme(legend.box = "vertical",
                                legend.box.just = 'left',
                                legend.justification = c(0, 0))

emat_plt <- ggpubr::ggarrange(emat_up_plt, emat_down_plt, nrow=1, align = 'hv')
emat_plt <- annotate_figure(emat_plt, top = text_grob("Cell markers enrichment (CellMarker 2.0)", size = 8))

outfile  <- paste0(path_results, "/emat_plt.pdf")
pdf(file = outfile, paper = "a4", w = unit(4,'cm'), h = unit(3.14,'cm'))
print(emat_plt)
dev.off()

outfile  <- paste0(path_results, "/emat_sel.gz")
write.table(emat_sel, file = gzfile(outfile), row.names = F,
            col.names = T, sep = "\t", quote = F)

outfile  <- paste0(path_results, "/emat.xlsx")
saveXLSresEdgeR2(emat, name = "CellMarker2", outfile = outfile)

# new gene overlap plt ----
emat_sel$sign <- gsub("down","-1",emat_sel$direction)
emat_sel$sign <- gsub("up","1",emat_sel$sign)
emat_sel$sign <- as.integer(emat_sel$sign)
# emat_sel$fdr[emat_sel$fdr>0.05] <- 1
emat_sel$score <- with(emat_sel, -log10(fdr)*sign)
emat_sel$Var2 <- factor(emat_sel$Var2, levels=sel_name[c(1,2,5,3,4,6)])
emat_sel$direction <- gsub("up","Up-regulated",emat_sel$direction)
emat_sel$direction <- gsub("down","Down-regulated",emat_sel$direction)
emat_sel$direction <- factor(emat_sel$direction, levels = c("Down-regulated","Up-regulated"))
emat_sel$cond <- gsub("d","D",emat_sel$cond)
emat_sel$sig <- emat_sel$fdr<0.05
emat_plt2 <- ggplot(emat_sel, aes(x=cond,y=Var2, size=or, fill=score, col=sig)) +
  geom_point(shape=21) + facet_wrap(~direction) +
  scale_size(range = c(2, 6), breaks = c(5, 15, 25), limits = c(0,53)) +
  xlab(NULL) + ylab(NULL) +
  scale_fill_gradient2(
    low = "#54278F",
    mid = 'white',
    high = "#D94801",
    # limits = c(-5, 5),
    # na.value = 'white'
  ) +
  scale_color_manual(values = c('white','black')) +
  guides(
    fill = guide_colorbar(barwidth = 3, barheight = .6, order = 2),
    size = guide_legend(order = 1, ncol = 1)
  ) +
  theme_bw() + my_theme + theme(strip.background = element_blank())

outfile  <- paste0(path_results, "/emat_plt2.pdf")
pdf(file = outfile, paper = "a4", w = unit(2.8,'cm'), h = unit(2.58,'cm'))
print(emat_plt2)
dev.off()
# plot dge expression ----
markers_sel <- list(
  "Pluripotent stem cell" = c("SOX2","POU5F1","NANOG","LIN28A"),
  "Epithelial cell" = c("CDH1","CDH3","EPCAM","CLDN6","CLDN3","ESRP1","ESRP2","TSPAN7","ITGA6"),
  "Epithelial progenitor cell" = c("PODXL","CD24"),
  "Mesenchymal cell" = c("CDH2","TWIST1","CDH11","VIM","FN1","ZEB1", "ZEB2"),
  "Neural progenitor cell" = c("SOX9","NES","RBFOX2","PAX3")
)
int_nl[int_nl$L1 == "Epithelial cell",]

markers_sel <- unlist(markers_sel)
names(markers_sel) <- gsub("\\d+$","",names(markers_sel))
fc_mat <- dge_toptable[dge_toptable$GeneID%in%markers_sel,]
fc_mat <- reshape2::dcast(fc_mat, contrast~GeneID, value.var = 'logFC', fun.aggregate = mean)
fc_mat$contrast <- gsub("ZNF398_I_2D_Day7_vs_ZNF398_NI_2D_Day7","2D",fc_mat$contrast)
fc_mat$contrast <- gsub("ZNF398_I_3D_Day4_vs_ZNF398_NI_3D_Day4","3D",fc_mat$contrast)
fc_mat <- fc_mat[1:2,]
rownames(fc_mat) <- fc_mat[,1]
fc_mat <- fc_mat[,-1]
fc_mat <- as.matrix(fc_mat)

fc_hm <- ComplexHeatmap::Heatmap(
  fc_mat[,markers_sel],
  col = c("#54278F","white","#D94801"),
  show_row_dend = F,
  show_column_dend = F,
  border = T,
  column_split = factor(names(markers_sel), levels=unique(names(markers_sel))[c(5,4,3,2,1)]),
  cluster_column_slices = F,
  row_names_gp = gpar(fontsize=8),
  column_names_rot = 45,
  column_names_gp = gpar(fontsize=8),
  row_title_gp = gpar(fontsize=8),
  column_title_gp = gpar(fontsize=8),
  heatmap_legend_param = list(title='log2FC[KDvsCtrl]',
                              title_gp=gpar(fontsize=8,face='plain'),
                              labels_gp=gpar(fontsize=8,face='plain'),
                              legend_direction='horizontal')
)

outfile  <- paste0(path_results, "/fc_hm.pdf")
message('-- saving to: ',outfile)
pdf(file = outfile, paper = "a4r", w = unit(6,'cm'), h = unit(1.65,'cm'))
draw(fc_hm, heatmap_legend_side='bottom')
dev.off()

# summary table significant ----
int_nl_all <- GeneOverlap::getNestedList(cell_name_set_gom, name = 'intersection')
int_nl_all <- reshape2::melt(int_nl_all)
colnames(int_nl_all) <- c("gene_name","condition","gene_set")
int_nl_all <- ddply(
  int_nl_all,
  .(condition, gene_set),
  summarize,
  genes = paste0(gene_name, collapse = ";")
)
colnames(emat)[1:2] <- c("condition","gene_set")
emat_sig <- emat[emat$fdr<0.05,]
emat_sig <- merge(emat_sig, int_nl_all, by = c("condition","gene_set"), all.x=T)

outfile  <- paste0(path_results, "/emat.fdr_0.05.xlsx")
saveXLSresEdgeR2(emat_sig, name = "CellMarker2", outfile = outfile)
