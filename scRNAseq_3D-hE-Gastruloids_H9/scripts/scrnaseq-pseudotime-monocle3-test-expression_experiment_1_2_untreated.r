#
# Monocle3 Graph test
#
# 0. Resources ----
source("scripts/resources.R")
suppressPackageStartupMessages(library(monocle3))

# paths
path_cds     = "results/correct-batches_experiment_1_2_untreated/cluster-cells/pseudotime-monocle3/cds.rds"

# parameters
q_value_cutoff  = 0.01
branch_1_vector = c(1,4,3)
branch_1_label  = "PSC_to_PS-like"
branch_2_vector = c(1,4,5,6)
branch_2_label  = "PSC_to_PostEpi-like"

# 1. Read data ----
cds <- readRDS(path_cds)

# 2. Perform test ----
colData(cds)$Branch_1 <- "other"
colData(cds)$Branch_1[colData(cds)$Cluster%in%branch_1_vector] <- branch_1_label
colData(cds)$Branch_2 <- "other"
colData(cds)$Branch_2[colData(cds)$Cluster%in%branch_2_vector] <- branch_2_label

pr_graph_test_res_1 <- monocle3::graph_test(cds[,colData(cds)$Branch_1==branch_1_label], neighbor_graph="principal_graph", cores=8)
rownames(pr_graph_test_res_1) <- pr_graph_test_res_1$gene_short_name
pr_graph_test_res_1_sig <- subset(pr_graph_test_res_1, q_value < q_value_cutoff)
pr_graph_test_res_2 <- monocle3::graph_test(cds[,colData(cds)$Branch_2==branch_2_label], neighbor_graph="principal_graph", cores=8)
rownames(pr_graph_test_res_2) <- pr_graph_test_res_2$gene_short_name
pr_graph_test_res_2_sig <- subset(pr_graph_test_res_2, q_value < q_value_cutoff)

# save results
pr_graph_test_res <- list(pr_graph_test_res_1,pr_graph_test_res_2)
names(pr_graph_test_res) <- c(branch_1_label,branch_2_label)

outfile <- "seurat/results_pseudotime-monocle3/pr_graph_test_res.rds"
message(" -- saving to: ", outfile)
saveRDS(pr_graph_test_res, file = outfile)