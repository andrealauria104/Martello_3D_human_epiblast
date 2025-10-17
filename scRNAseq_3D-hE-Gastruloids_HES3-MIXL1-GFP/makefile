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

# ------------------------------------------------------------------------------
# Collect data
# ------------------------------------------------------------------------------
# Raw counts
data/raw_counts.csv.gz: ../data/raw_counts.csv.gz
	mkdir -p `dirname $@`
	cd `dirname $@`;\
	ln -s  ../$< .

# Metadata
data/metadata.txt.gz: ../data/metadata.txt.gz
	mkdir -p `dirname $@`
	cd `dirname $@`;\
	ln -s  ../$< .

data/gene_info.gencode.hsapiens.32.primary_assembly.annotation.txt.gz:
	mkdir -p `dirname $@`
	cd `dirname $@`;\
	ln -s ../../../../Martello_Polo_TGFb/data/gene_info.gencode.hsapiens.32.primary_assembly.annotation.txt.gz .
