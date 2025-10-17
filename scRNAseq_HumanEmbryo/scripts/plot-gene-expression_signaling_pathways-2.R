setwd("/home/epigen/Martello_3D_human_epiblast/Martello_Polo_HumanEmbryo/analysis_oscab")
# 0. Resources ----
source("scripts/resources.R")

# paths
path_results = "results/xiang_et_al_human_embryo/signaling_pathways"
path_gene_list = "data/xiang_et_al_human_embryo/Ligands_receptors_signalling_pathways_human.xlsx"
path_gene_list_sel = "data/xiang_et_al_human_embryo/Ligands_receptors_signalling_pathways_expr_selection.txt"

if(!dir.exists(path_results)) dir.create(path_results)

# 1. Read data ----
sce <- readRDS("results/xiang_et_al_human_embryo/preprocess/sce.rds")

gene_list <- openxlsx::read.xlsx(path_gene_list, sheet = 1)
colnames(gene_list)[1] <- "type"
gene_list_l <- list()
for(i in colnames(gene_list)[-1]) {
  gene_list_l[[i]] <- gene_list[,c("type",i)]
  gene_list_l[[i]] <- gene_list_l[[i]][!is.na(gene_list_l[[i]][,2]),]
  dup_idx <- which(duplicated(gene_list_l[[i]][,2]))
  if(length(dup_idx)!=0) gene_list_l[[i]] <- gene_list_l[[i]][-dup_idx,]
  colnames(gene_list_l[[i]])[2] <- "gene_name"
  gene_list_l[[i]][["pathway"]] <- i
}
gene_list_l <- do.call(rbind,gene_list_l)
gene_list_l <- gene_list_l[gene_list_l$pathway %in% c("FGF","INSULIN","TGF-beta"),]
rownames(gene_list_l) <- NULL
gene_list_l <- split(
  gene_list_l,
  gene_list_l$type
)

to_fill <- data.frame(
  "type"="Receptors",
  "gene_name"="TDGF1",
  "pathway"="TGF-beta"
)
gene_list_l$Receptors <- rbind(gene_list_l$Receptors, to_fill)

# 2. Retrieve gene expression ----
cells_for_receptors <- sce$Sample[sce$CellType=="Post-EPIs"]
cells_for_ligands <- sce$Sample[sce$CellType%in%c("Post-EPIs","Hypoblast","CTBs")]
m <- assay(sce,"logcounts")

process_mat <- function(m, genes, cells) 
{
  mat <- as.data.frame(m[rownames(m) %in% genes$gene_name,colnames(m) %in% cells])
  mat <- cbind("gene_name"=rownames(mat),mat)
  mat <- reshape2::melt(
    mat,
    id.var="gene_name",
    variable.name = "Sample",
    value.name = "expression")
  mat <- as.data.frame(merge(mat,colData(sce),by="Sample"))
  mat <- merge(mat, genes, by="gene_name")
  mat$gene_name <- factor(
    mat$gene_name,
    levels = unique(
      stringr::str_sort(mat$gene_name, numeric = T)
    )
  )
  return(mat)
}

# receptors
mat_receptors <- process_mat(
  m = m,
  genes = gene_list_l$Receptors,
  cells = cells_for_receptors
)
mat_receptors$gene_name <- factor(
  mat_receptors$gene_name,
  levels = unique(mat_receptors$gene_name)[c(1:12,14:16,13)]
  )

# ligands
mat_ligands <- process_mat(
  m = m,
  genes = gene_list_l$Ligands,
  cells = cells_for_ligands
)


# 3. Plot expression by grid ----
plot_by_grid <- function(mat, title,pal=NULL) 
{
  expr_plt <- ggplot(
    mat,
    aes(x=gene_name,y=expression,col=CellType)
    ) + 
    geom_boxplot(outlier.shape = NA, fill="white", show.legend = F) +
    ggrastr::rasterise(geom_jitter(width = .2,height = 0,size=0.3),dpi=300) +
    facet_grid(CellType~pathway,scales = "free",space = "free_x" ) +
    theme_bw() +
    my_theme + 
    xlab(NULL) + ylab("Expression") +
    theme(
      axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
      aspect.ratio = NULL
      ) +
    ggtitle(title) +
    scale_y_continuous(
      labels = scales::number_format(accuracy = 0.1)
    ) +
    guides(col=guide_legend(override.aes=list(size = 3))) +
    scale_color_manual(values = pal)
  
  return(expr_plt)
}

pal_embryo <- c(RColorBrewer::brewer.pal(9,"Purples")[c(4,6,9)]
                ,RColorBrewer::brewer.pal(8,"Dark2")[c(1,4)]
                ,RColorBrewer::brewer.pal(9,"Oranges")[c(3,6,8)])
names(pal_embryo) <- c("Pre-EPIs", "Post-EPIs","PSA-EPI"
                       ,"ICM","Hypoblast"
                       ,"CTBs","STBs","EVTs")
# receptors
expr_receptors_grid_plt <- plot_by_grid(
  mat = mat_receptors,
  title = "Signaling : Receptors",
  pal = pal_embryo
)
outfile <- paste0(path_results,"/expr_receptors_grid_plt_2.pdf")
message(" -- saving to: ", outfile)
pdf(file = outfile, paper = "a4r", width = unit(5.3,'cm'), height = unit(1.8,'cm'))
print(expr_receptors_grid_plt)
dev.off()
