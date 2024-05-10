#
# Find markers by comparing cells in one group vs all the remaining cells
#
# 0. Resources ----
source("scripts/resources.R")
suppressPackageStartupMessages(library(matrixTests))

# paths
path_sce = "results/correct-batches_experiment_1_2_untreated_human_embryo/cluster-cells/sce.rds"
path_results = "results/correct-batches_experiment_1_2_untreated_human_embryo/cluster-cells/markers-one-vs-rest"

if(!dir.exists(path_results)) dir.create(path_results, recursive = T)

# 0.1 Global params ----
assay.type = "logcounts"
group.by = "Cluster"

padj.th = 0.05
pvalue.th = 0.01
min.pct = 0.1
# 1. Read and prepare data ----
sce <- readRDS(path_sce)

groups <- stringr::str_sort(
  as.character(
    unique(sce[[group.by]])
  ),
  numeric = T
)
cells.by.group <- sce[[group.by]]
names(cells.by.group) <- colnames(sce)

mat <- assay(sce, assay.type)

# 2. Run test ----
de.genes.by.group <- list()
for(i in 1:length(groups)) {
  
  # Filter lowly detected genes
  x = mat[,cells.by.group==groups[i]]
  # x.which.detected <- rowSums(x>0)>=ceiling(min.pct*ncol(x))
  perc.x <- round(rowSums(x>0)/ncol(x),3)
  
  y = mat[,cells.by.group!=groups[i]]
  # y.which.detected <- rowSums(y>0)>=ceiling(min.pct*ncol(y))
  perc.y <- round(rowSums(y>0)/ncol(y),3)
  
  # min.pct in at least 1 group
  # detected.genes.by.group <- which(x.which.detected | y.which.detected)
  detected.genes.by.group <- which(pmax(perc.x,perc.y) >= min.pct)
  
  # Run wilcoxon
  de.genes.by.group[[i]] <- matrixTests::row_wilcoxon_twosample(
    x = x[detected.genes.by.group,],
    y = y[detected.genes.by.group,],
    alternative = "g"
  )
  de.genes.by.group[[i]] <- de.genes.by.group[[i]][
    order(de.genes.by.group[[i]]$pvalue, decreasing = F),
  ]
  de.genes.by.group[[i]]$fdr <- p.adjust(
    de.genes.by.group[[i]]$pvalue,
    method = "fdr",
    n = nrow(mat)
    )
  de.genes.by.group[[i]]$pvalue_adj <- p.adjust(
    de.genes.by.group[[i]]$pvalue,
    method = "bonferroni",
    n = nrow(mat)
  )
}
names(de.genes.by.group) <- groups

# 3. Save results ----
outfile <- paste0(path_results,"/de.genes.by.group.rds")
message(" -- saving: ",outfile)
saveRDS(de.genes.by.group, file = outfile)

# Filter by pvalue
de.genes.by.group.sig <- lapply(
  de.genes.by.group,
  function(x) {
    x <- x[x$pvalue<pvalue.th,] 
    x <- cbind.data.frame("gene_name"=rownames(x),x)
    return(x)
    }
  )
outfile <- paste0(path_results,"/de.genes.by.group.pvalue_",pvalue.th,".xlsx")
saveXLSresEdgeR2(res = de.genes.by.group.sig, outfile = outfile)

# Filter by pvalue adjusted
de.genes.by.group.sig <- lapply(
  de.genes.by.group,
  function(x) {
    x <- x[x$pvalue_adj<padj.th,]
    x <- cbind.data.frame("gene_name"=rownames(x),x)
    return(x)
  }
)
outfile <- paste0(path_results,"/de.genes.by.group.pvalue_adj_",padj.th,".xlsx")
saveXLSresEdgeR2(res = de.genes.by.group.sig, outfile = outfile)
