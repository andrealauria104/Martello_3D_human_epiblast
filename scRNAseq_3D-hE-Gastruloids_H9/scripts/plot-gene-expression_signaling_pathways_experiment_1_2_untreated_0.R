setwd("/home/epigen/Martello_3D_human_epiblast/Martello_Polo_TGFb/analysis_oscab")
# 0. Resources ----
source("scripts/resources.R")

# paths
path_results = "results/correct-batches_experiment_1_2_untreated_0/signaling_pathways"
path_gene_list = "data/xiang_et_al_human_embryo/Ligands_receptors_signalling_pathways_human.xlsx"
path_gene_list_sel = "data/xiang_et_al_human_embryo/Ligands_receptors_signalling_pathways_expr_selection.txt"
path_extended_data = "results/correct-batches_experiment_1_2_untreated_0/cluster-cells/extended_data.rds"

if(!dir.exists(path_results)) dir.create(path_results)

# 1. Read data ----
sce <- readRDS("results/correct-batches_experiment_1_2_untreated_0/cluster-cells/sce.rds")
extended_data <- readRDS(path_extended_data)
palette_features <- extended_data$palette_features

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

rownames(gene_list_l) <- NULL
gene_list_l <- split(
  gene_list_l,
  gene_list_l$type
)
gene_list_l$Ligands <- rbind(gene_list_l$Ligands,
      data.frame("type"="Ligands",
                 "gene_name"= "TDGF1",
                 "pathway" = "TGF-beta")
      )
# 2. Retrieve gene expression ----
# cells_for_receptors <- sce$Sample[sce$CellType=="Post-EPIs"]
# cells_for_ligands <- sce$Sample[sce$CellType%in%c("Post-EPIs","Hypoblast","CTBs")]
m <- assay(sce,"logcounts")

process_mat <- function(m, genes, cells=NULL) 
{
  mat <- as.data.frame(m[rownames(m) %in% genes$gene_name,])
  if(!is.null(cells)) mat <- m[,colnames(m) %in% cells]
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
  # cells = cells_for_receptors
)

# ligands
mat_ligands <- process_mat(
  m = m,
  genes = gene_list_l$Ligands,
  # cells = cells_for_ligands
)

# 3. Plot expression by grid ----
plot_by_grid <- function(mat, title,pal=NULL,group.by="CellType") 
{
  expr_plt <- ggplot(
    mat,
    aes(x=gene_name,y=expression,col=.data[[group.by]])
    ) + 
    geom_boxplot(outlier.shape = NA, fill="white", show.legend = F) +
    ggrastr::rasterise(geom_jitter(width = .2,height = 0,size=0.3),dpi=300) +
    facet_grid(as.formula(paste0(group.by,"~pathway")),scales = "free",space = "free_x" ) +
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

# receptors
expr_receptors_grid_plt <- plot_by_grid(
  mat = mat_receptors,
  title = "Signaling : Receptors",
  pal = palette_features$ClusterNum,
  group.by = "ClusterNum"
)
outfile <- paste0(path_results,"/expr_receptors_grid_plt.pdf")
message(" -- saving to: ", outfile)
pdf(file = outfile, paper = "a4r", width = unit(15,'cm'), height = unit(15,'cm'))
print(expr_receptors_grid_plt)
dev.off()

# ligands
expr_ligands_grid_plt <- plot_by_grid(
  mat = mat_ligands,
  title = "Signaling : Ligands",
  pal = palette_features$ClusterNum,
  group.by = "ClusterNum"
)
outfile <- paste0(path_results,"/expr_ligands_grid_plt.pdf")
message(" -- saving to: ", outfile)
pdf(file = outfile, paper = "a4r", width = unit(15,'cm'), height = unit(15,'cm'))
print(expr_ligands_grid_plt)
dev.off()

# 4. Plot expression by wrap ----
plot_by_wrap <- function(mat, title,pal=NULL) 
{
  expr_plt <- ggplot(
    mat,
    aes(x=ClusterNum,y=expression,col=ClusterNum)
  ) + 
    geom_boxplot(outlier.shape = NA, fill="white", show.legend = F) +
    ggrastr::rasterise(geom_jitter(width = .2,height = 0,size=0.3),dpi=300) +
    facet_wrap(~gene_name,ncol=5) +
    theme_bw() +
    my_theme + 
    xlab(NULL) + ylab("Expression") +
    theme(
      aspect.ratio = 1
    ) +
    ggtitle(title) +
    scale_y_continuous(
      labels = scales::number_format(accuracy = 0.1)
    ) +
    guides(col=guide_legend(override.aes=list(size = 3))) +
    scale_color_manual(values = pal)
  
  return(expr_plt)
}


for(i in unique(mat_ligands$pathway)) {
  mat.to.plot <- mat_ligands[mat_ligands$pathway==i,]
  outfile <- paste0(path_results,"/expr_ligands_wrap.",i,".pdf")
  message(" -- saving to: ", outfile)
  pdf(file = outfile, paper = "a4r", width = unit(8,'cm'), height = unit(8,'cm'))
  print(
    plot_by_wrap(
      mat.to.plot,
      title = paste0(i," : Ligands"),
      pal = palette_features$ClusterNum
    )
  )
  dev.off()
}

for(i in unique(mat_receptors$pathway)) {
  mat.to.plot <- mat_receptors[mat_receptors$pathway==i,]
  outfile <- paste0(path_results,"/expr_receptors_wrap.",i,".pdf")
  message(" -- saving to: ", outfile)
  pdf(file = outfile, paper = "a4r", width = unit(8,'cm'), height = unit(8,'cm'))
  print(
    plot_by_wrap(
      mat.to.plot,
      title = paste0(i," : Receptors"),
      pal = palette_features$ClusterNum
    )
  )
  dev.off()
}


