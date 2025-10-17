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

# -------------- #
# OSCAB analysis #
# -------------- #
results/xiang_et_al_human_embryo/preprocess/sce.rds: scripts/preprocess-scrnaseq_xiang_et_al_human_embryo.R \
                                           	     data/xiang_et_al_human_embryo/raw_counts.csv.gz data/xiang_et_al_human_embryo/metadata.txt.gz 
	$(CONDA_ACTIVATE) rstudio_Rv4.0.3;\
	Rscript --vanilla $<
	touch $@

# ------------------- #
# Integrated analysis #
# ------------------- #

