# Analyze RNAseq QuantSeq
setwd("/home/epigen/Martello_3D_human_epiblast/RNAseq/dataset/v2-CRISPRi_ZNF398_lextrim")
# 0. Resources ----
# Tools -- 
suppressPackageStartupMessages(library(RNAseqRtools))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(enrichR))
suppressPackageStartupMessages(library(ggrastr))

# params --
expression.unit = "cpm"
group = "condition_sepgr"
group_ref = "P_NI_2D_Day7"
shape_by = "Culture"

nvar = 2000

# paths --
path.counts = "GEP.count.gz"
path.metadata = "metadata.txt"
path.gene.info = "/home/reference_data/bioinfotree/task/gencode/dataset/hsapiens/32/primary_assembly.annotation.ensg2gene_symbol2biotype.map.header_added"

path.results = "results"

if(!dir.exists(path.results)) dir.create(path.results)

s33d = 1991
# 1. Load data ----
# counts
counts <- read.delim(
  path.counts,
  stringsAsFactors = F,
  header = T,
  sep = "\t",
  row.names = 1
)

# metadata
metadata <- read.delim(
  path.metadata,
  stringsAsFactors = F,
  header = T,
  sep = "\t",
)
rownames(metadata) <- gsub("\\-","\\.",metadata[,1])
metadata <- metadata[colnames(counts),]
message("sample metadata names check: ",
        sum(rownames(metadata) == colnames(counts)) == ncol(counts)
)

# gene metadata
gene_info  <- read.delim(
  path.gene.info,
  stringsAsFactors = F,
  header = T,
  sep = "\t",
)
gene_info <- plyr::ddply(
  gene_info,
  .(gene_name),
  summarize,
  gene_id=paste0(gene_id, collapse = ";"),
  gene_type=paste0(gene_type, collapse = ";")
)
rownames(gene_info) <- gene_info$gene_name
gene_info <- gene_info[rownames(counts),]
message("gene metadata names check: ",
        sum(rownames(gene_info) == rownames(counts)) == nrow(counts)
)

# create edger object
y <- RNAseqRtools::processRNAseqEdgeR(
  m = counts,
  experimental_info = metadata,
  gene_info = gene_info,
  group = group,
  reference = group_ref,
  filter.expr.th = 1,
  filter.sample.th = min(table(metadata$Group))
  )

message(" -- writing to: ",paste0(path.results,"/normcounts_tmm.logcpm.txt.gz"))
write.table(
  cbind.data.frame(y$genes,y$logCPM),
  file = gzfile(paste0(path.results,"/normcounts_tmm.logcpm.txt.gz")),
  row.names = F,
  col.names = T,
  sep = "\t",
  quote = F
  )

message(" -- writing to: ",paste0(path.results,"/edger.processed.counts.rds"))
saveRDS(y,paste0(path.results,"/edger.processed.counts.rds"))

pal_analysis <- c(
    "P_NI_2D_Day7" = "#FFCCFF",
    "ZNF398_NI_2D_Day7" = "#C6DBEF",
    "ZNF398_I_2D_Day7" = "#D94801" ,
    "ZNF398_NI_3D_Day4" = "#2171B5",
    "ZNF398_I_3D_Day4" = "#7F2704" 
    )

# 2. Sample clustering ----
expr.m <- y$logCPM
colnames(expr.m) <- paste0(y$samples$GroupRed," #",c(1,2))

sample.cor.cluster.hm <- RNAseqRtools::plotCorrelation(
  expr.m,
  reorder_cormat_by_cor_dist = T,
  draw_cor_values = T,
  cor_values_size = 4.5
)

figout  <- paste0(path.results, "/sample.cor.cluster.hm.pdf")
pdf(file = figout, paper = "a4", w = unit(4,'cm'), h = unit(4.5,'cm'))
draw(sample.cor.cluster.hm)
dev.off()

sample.cor.cluster.hm2 <- RNAseqRtools::plotCorrelation(
  expr.m,
  reorder_cormat_by_cor_dist = T,
  draw_cor_values = F,
  column_split = 4,
  split = 4
)

figout  <- paste0(path.results, "/sample.cor.cluster.hm2.pdf")
pdf(file = figout, paper = "a4", w = unit(4,'cm'), h = unit(4.5,'cm'))
draw(sample.cor.cluster.hm2)
dev.off()

# 3. PCA ----
y$samples$condition_sepgr <- factor(y$samples$condition_sepgr, levels=names(pal_analysis))
pca <- RNAseqRtools::plotPCA(
  y$logCPM,
  experimental_info = y$samples,
  col_by = group,
  shape_by = shape_by,
  labels = F,
  scree_plot = F,
  scree_plot_type = "standard_deviation",
  point_size = 2,
  legend_position = "right",
  default.title = "PCA: RNAseq [CRISPRi ZNF398]",
  pal = pal_analysis,
  col_legend_nrow = 5,
  shape_legend_nrow = 2
  ) + theme(plot.title = element_text(face="plain",size=8))

figout  <- paste0(path.results, "/pca.pdf")
pdf(file = figout, paper = "a4", w = unit(4.35,'cm'), h = unit(4,'cm'))
print(pca)
dev.off()

vargenes <- RNAseqRtools::find.variable.genes(y$CPM, n = nvar)
pca.vargenes <- RNAseqRtools::plotPCA(
  y$logCPM[vargenes$genes,],
  experimental_info = y$samples,
  col_by = group,
  shape_by = shape_by,
  labels = F,
  scree_plot = F,
  scree_plot_type = "standard_deviation",
  point_size = 2,
  legend_position = "right",
  pal = pal_analysis,
  col_legend_nrow = 5,
  shape_legend_nrow = 2
) +
  ggtitle(
    label = "PCA: RNAseq [CRISPRi ZNF398]",
    subtitle = paste0("n. variable genes = ",nvar)) +
  theme(
    plot.title = element_text(face="plain",size=8),
    plot.subtitle = element_text(face="plain",size=7, hjust = 1),
    panel.grid.major = element_blank(),
    ) +
  guides(col=guide_legend(title = "Condition"))

figout  <- paste0(path.results, "/pca.vargenes.n_",nvar,".pdf")
pdf(file = figout, paper = "a4", w = unit(4.,'cm'), h = unit(3.6,'cm'))
print(pca.vargenes)
dev.off()

vargenes.sel <- RNAseqRtools::find.variable.genes(y$CPM[,y$samples$condition_sepgr!="P_NI_2D_Day7"], n = nvar)
pca.vargenes.sel <- RNAseqRtools::plotPCA(
  y$logCPM[vargenes$genes,y$samples$condition_sepgr!="P_NI_2D_Day7"],
  experimental_info = y$samples,
  col_by = group,
  shape_by = shape_by,
  labels = F,
  scree_plot = F,
  scree_plot_type = "standard_deviation",
  point_size = 2,
  legend_position = "right",
  pal = pal_analysis,
  col_legend_nrow = 5,
  shape_legend_nrow = 2
) +
  ggtitle(
    label = "PCA: RNAseq [CRISPRi ZNF398]",
    subtitle = paste0("n. variable genes = ",nvar)) +
  theme(
    plot.title = element_text(face="plain",size=8),
    plot.subtitle = element_text(face="plain",size=7, hjust = 1),
    panel.grid.major = element_blank(),
  ) +
  guides(col=guide_legend(title = "Condition"))

figout  <- paste0(path.results, "/pca.vargenes.n_",nvar,".sel.pdf")
pdf(file = figout, paper = "a4", w = unit(4.,'cm'), h = unit(3.6,'cm'))
print(pca.vargenes.sel)
dev.off()

# 4. DGE ----
path.dge.toptable <- "DGE_sep_group/edger.toptable_clean.ALL_contrast.mark_seqc.header_added.gz"
dge.toptable <- read.delim(path.dge.toptable)

volcano.plt <- list()
for(i in unique(dge.toptable$contrast)) {
  de.table <- dge.toptable[dge.toptable$contrast==i,]
  colnames(de.table)[5] <- "FDR"
  de.table$logCPM <- 1
  rownames(de.table) <- de.table$GeneID
  volcano.plt[[i]] <- RNAseqRtools::plotDiffExprRes(
    de.table,
    type = "volcano",
    pal=c("#54278F","lightgrey","#D94801"),
    gtitle = unique(de.table$contrast),
    lfcTh = 0.5,
    point.size = .2,
  ) +
    theme(plot.title = element_text(face="plain"),
          panel.grid.major = element_blank()) +
    guides(col=guide_legend(override.aes = list(size=3)))
  volcano.plt[[i]] <- ggrastr::rasterise(volcano.plt[[i]], dpi=400)
}

volcano.plt.arr <- ggpubr::ggarrange(plotlist = volcano.plt,
                                     nrow=1, align = 'hv', common.legend = T)
figout  <- paste0(path.results, "/volcano.plt.arr.pdf")
pdf(file = figout, paper = "a4r", w = unit(7.6,'cm'), h = unit(2.8,'cm'))
print(volcano.plt.arr)
dev.off()

# 5. DGE overlap ----
dge.lists <- list(
  "up_reg_2d"=dge.toptable[dge.toptable$contrast=="ZNF398_I_2D_Day7_vs_ZNF398_NI_2D_Day7" & dge.toptable$significance==1,'GeneID'],
  "down_reg_2d"=dge.toptable[dge.toptable$contrast=="ZNF398_I_2D_Day7_vs_ZNF398_NI_2D_Day7" & dge.toptable$significance==-1,'GeneID'],
  "up_reg_3d"=dge.toptable[dge.toptable$contrast=="ZNF398_I_3D_Day4_vs_ZNF398_NI_3D_Day4" & dge.toptable$significance==1,'GeneID'],
  "down_reg_3d"=dge.toptable[dge.toptable$contrast=="ZNF398_I_3D_Day4_vs_ZNF398_NI_3D_Day4" & dge.toptable$significance==-1,'GeneID']
)
# plot_venn_diagram(dge.lists[c(1,3)])
# plot_venn_diagram(dge.lists[c(2,4)])

dge.summary <- data.frame(
  "gene_name" = sort(unique(unlist(dge.lists)))
)
dge.summary$up_in_both <- dge.summary$gene_name %in% intersect(dge.lists$up_reg_2d,dge.lists$up_reg_3d)
dge.summary$up_in_2d <- dge.summary$gene_name %in% setdiff(dge.lists$up_reg_2d,dge.lists$up_reg_3d)
dge.summary$up_in_3d <- dge.summary$gene_name %in% setdiff(dge.lists$up_reg_3d,dge.lists$up_reg_2d)

dge.summary$down_in_both <- dge.summary$gene_name %in% intersect(dge.lists$down_reg_2d,dge.lists$down_reg_3d)
dge.summary$down_in_2d <- dge.summary$gene_name %in% setdiff(dge.lists$down_reg_2d,dge.lists$down_reg_3d)
dge.summary$down_in_3d <- dge.summary$gene_name %in% setdiff(dge.lists$down_reg_3d,dge.lists$down_reg_2d)

dge.summary.count <- reshape2::melt(colSums(dge.summary[,-1]))
dge.summary.count$status <- rownames(dge.summary.count)
dge.summary.count$direction <- paste0(gsub("_.*","",rownames(dge.summary.count)),"-regulated")

dge.summary.count$status <- factor(dge.summary.count$status, levels = c(
  "up_in_both","up_in_2d","up_in_3d",
  "down_in_both", "down_in_2d","down_in_3d"
))
dge.summary.count$status_alpha <- gsub(".*_in_","",dge.summary.count$status)
dge.summary.count$status_alpha <- gsub("2d","2d-only",dge.summary.count$status_alpha)
dge.summary.count$status_alpha <- gsub("3d","3d-only",dge.summary.count$status_alpha)
dge.summary.count$status_alpha <- factor(dge.summary.count$status_alpha, levels = c("both","2d-only","3d-only"))

pal.dge.summary <- c("#D94801","#D94801","#D94801",
                     "#54278F", "#54278F", "#54278F")
names(pal.dge.summary) <- levels(dge.summary.count$status)
dge.summary.plt <- ggplot(dge.summary.count, aes(x=status_alpha, y=value,fill=status, label=value)) +
  geom_col(aes(alpha=status_alpha, col=status_alpha), show.legend = F) + coord_flip() +
  facet_wrap(~direction, ncol=2) +
  geom_text(size=2.5, nudge_y = 85, show.legend = F) +
  xlab(NULL) + ylab("N. of genes (x 100)") +
  scale_fill_manual(values = pal.dge.summary) +
  scale_color_manual(values = c("2d-only" = 'white', "3d-only" = 'white', "both" = 'black')) +
  scale_alpha_manual(values = c("2d-only" = 0.6, "3d-only" = 1, "both" = 0.3)) +
  theme_classic() + my_theme +
  scale_y_continuous(labels = function(x) x/100) +
  theme(panel.grid.major = element_blank())
  
figout  <- paste0(path.results, "/dge.summary.plt.pdf")
pdf(file = figout, paper = "a4", w = unit(3.6,'cm'), h = unit(1.5,'cm'))
print(dge.summary.plt)
dev.off()

outfile  <- paste0(path.results, "/dge.summary.gz")
write.table(dge.summary, gzfile(outfile), row.names = F,
            col.names = T, sep = "\t", quote = F)
# 6. Markers expression ----
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
  de.2d.idx <- genes %in% dge.toptable[dge.toptable$contrast=="ZNF398_I_2D_Day7_vs_ZNF398_NI_2D_Day7" & dge.toptable$significance!=0,'GeneID']
  de.3d.idx <- genes %in% dge.toptable[dge.toptable$contrast=="ZNF398_I_3D_Day4_vs_ZNF398_NI_3D_Day4" & dge.toptable$significance!=0,'GeneID']
  genes.lab <- genes
  names(genes.lab) <- genes
  genes.lab[de.2d.idx] <- paste0(genes[de.2d.idx]," *")
  genes.lab[de.3d.idx] <- paste0(genes.lab[de.3d.idx]," °")
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

figout  <- paste0(path.results, "/plt.expr.marker_genes.pdf")
pdf(file = figout, paper = "a4", w = unit(8.5,'cm'), h = unit(11,'cm'))
print(plt.expr.arranged)
dev.off()

# 7. Common volcano ----
dge_lab_sel <- c(
  head(intersect(dge.lists$up_reg_2d,dge.lists$up_reg_3d),n=4),
  head(intersect(dge.lists$down_reg_2d,dge.lists$down_reg_3d),n=4)
)
volcano.comm.plt <- list()
for(i in unique(dge.toptable$contrast)[-1]) {
  de.table <- dge.toptable[dge.toptable$contrast==i,]
  colnames(de.table)[5] <- "FDR"
  de.table$logCPM <- 1
  rownames(de.table) <- de.table$GeneID
  volcano.comm.plt[[i]] <- RNAseqRtools::plotDiffExprRes(
    de.table,
    type = "volcano",
    pal=c("#54278F","lightgrey","#D94801"),
    gtitle = unique(de.table$contrast),
    label_selection = dge_lab_sel,
    lfcTh = 0.5,
    point.size = .2,
  ) +
    theme(plot.title = element_text(face="plain"),
          panel.grid.major = element_blank()) +
    guides(col=guide_legend(override.aes = list(size=3)))
  volcano.comm.plt[[i]]$layers[[5]] <- NULL
  volcano.comm.plt[[i]] <- volcano.comm.plt[[i]] +
    ggrepel::geom_text_repel(show.legend = F, size = 2,
                             segment.size = 0.1, max.overlaps = 1000)
  volcano.comm.plt[[i]] <- ggrastr::rasterise(volcano.comm.plt[[i]], dpi=600)
}

volcano.comm.plt.arr <- ggpubr::ggarrange(plotlist = volcano.comm.plt,
                                     nrow=1, align = 'hv', common.legend = T)
figout  <- paste0(path.results, "/volcano.comm.plt.arr.pdf")
pdf(file = figout, paper = "a4r", w = unit(5.2,'cm'), h = unit(2.7,'cm'))
print(volcano.comm.plt.arr)
dev.off()
