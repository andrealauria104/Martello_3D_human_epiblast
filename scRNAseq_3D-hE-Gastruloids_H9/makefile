# GNU make setup
SHELL:=/bin/bash

.DELETE_ON_ERROR:

# ------------- #
# Configuration #
# ------------- #
# conda
CONDA_ROOT?=/opt/conda
CONDA_VERSION?=miniconda3
CONDA_ACTIVATE=source $(CONDA_ROOT)/$(CONDA_VERSION)/etc/profile.d/conda.sh; conda activate
ENV_PREFIX?=$(CONDA_ROOT)/$(CONDA_VERSION)

define RUN_RSCRIPT
	$(CONDA_ACTIVATE) rstudio_Rv4.0.3;\
	Rscript --vanilla $<
	touch $@
endef

# -------------- #
# OSCAB analysis #
# -------------- #

# Preprocessing
results/experiment_1/sce.rds: scripts/preprocess-scrnaseq-experiment_1.R \
			      data/experiment_1/raw_counts.csv.gz data/experiment_1/metadata_full.txt \
			      data/gene_info.gencode.hsapiens.32.primary_assembly.annotation.txt.gz
	$(RUN_RSCRIPT)
results/experiment_2_untreated/sce.rds: scripts/preprocess-scrnaseq-experiment_2_untreated.R \
			                data/experiment_2_untreated/raw_counts.csv.gz data/experiment_2_untreated/metadata_full.txt \
					data/gene_info.gencode.hsapiens.32.primary_assembly.annotation.txt.gz
	$(RUN_RSCRIPT)

# Batch correction
results/correct-batches_experiment_1_2_untreated_0/sce.rds: scripts/scrnaseq-correct-batches_experiment_1_2_untreated_0.R \
							    results/experiment_1/sce.rds results/experiment_2_untreated/sce.rds
	$(RUN_RSCRIPT)

results/correct-batches_experiment_1_2_untreated/sce.rds: scripts/scrnaseq-correct-batches_experiment_1_2_untreated_0.R \
							  results/experiment_1/sce.rds results/experiment_2_untreated/sce.rds \
							  results/correct-batches_experiment_1_2_untreated_0/cluster-cells/sce.rds
	$(RUN_RSCRIPT)

# Clustering 
results/correct-batches_experiment_1_2_untreated_0/cluster-cells/sce.rds: scripts/scrnaseq-cluster-cells_experiment_1_2_untreated_0.R \
									  results/correct-batches_experiment_1_2_untreated_0/sce.rds
	$(RUN_RSCRIPT)

results/correct-batches_experiment_1_2_untreated/cluster-cells/sce.rds: scripts/scrnaseq-cluster-cells_experiment_1_2_untreated.R \
									  results/correct-batches_experiment_1_2_untreated/sce.rds \
								
	$(RUN_RSCRIPT)

# Pseudotime
results/correct-batches_experiment_1_2_untreated_0/cluster-cells/pseudotime-monocle3/cds.rds: scripts/scrnaseq-pseudotime-monocle3_experiment_1_2_untreated_0.R \
											      results/correct-batches_experiment_1_2_untreated_0/cluster-cells/sce.rds
	$(RUN_RSCRIPT)
results/correct-batches_experiment_1_2_untreated/cluster-cells/pseudotime-monocle3/cds.rds: scripts/scrnaseq-pseudotime-monocle3_experiment_1_2_untreated.R \
											      results/correct-batches_experiment_1_2_untreated/cluster-cells/sce.rds
	$(RUN_RSCRIPT)

# ------------------- #
# Integrated analysis #
# ------------------- #

