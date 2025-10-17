# setup ----
setwd("/home/epigen/Martello_3D_human_epiblast/public_datasets/Argelaguet_et_al_Nature_2019")
library(Seurat)

# read data ----
counts <- read.delim("GSE133725_embryo_cell_counts_matrix_030719.tsv.gz", check.names = FALSE)
metadata <- read.delim("sample_metadata.txt.gz")
rownames(metadata) <- metadata$sample
all(colnames(counts)[-1] %in% metadata$sample)

gene_info <- read.delim("/home/reference_data/bioinfotree/task/gencode/dataset/mmusculus/M23/primary_assembly.annotation.ensg2gene_symbol2biotype.map.header_added")
gene_info$ens_id <- gsub("\\.\\d+$","", gene_info$gene_id)
table(counts$ens_id %in% gene_info$ens_id)
counts <- merge(gene_info, counts, by = "ens_id", all.y = T)
counts[1:10,1:10]

# create object ----
counts_mat <- as.matrix(counts[, 5:ncol(counts)])
rownames(counts_mat) <- counts$ens_id

all(colnames(counts_mat) %in% metadata$sample)
metadata <- metadata[which(metadata$sample %in% colnames(counts_mat)),]
metadata <- metadata[metadata$pass_rnaQC, ]
counts_mat <- counts_mat[, metadata$sample]

identical(
  colnames(counts_mat),
  metadata$sample
)

identical(
  colnames(counts_mat),
  rownames(metadata)
)

seurat_obj <- Seurat::CreateSeuratObject(
  counts = counts_mat,
  project = "GSE133725",
  meta.data = metadata,
  min.cells = 3,
  min.features = 200
)

# run seurat pipeline ----
seurat_obj <- Seurat::NormalizeData(seurat_obj)
seurat_obj <- Seurat::FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
seurat_obj <- Seurat::ScaleData(seurat_obj, features = rownames(seurat_obj))
seurat_obj <- Seurat::RunPCA(seurat_obj)
seurat_obj <- Seurat::RunUMAP(seurat_obj, dims = 1:50)

# save ----
saveRDS(seurat_obj, "GSE133725.seurat_obj.rds")
