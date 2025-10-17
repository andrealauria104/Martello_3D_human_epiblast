## 3D Human Epiblast Model

Code repository for the paper: "A human epiblast model reveals dynamic TGF-beta mediated control of epithelial identity during mammalian epiblast development".

* Publication: 
* BioRxiv: TGF-beta dynamically controls epithelial identity in a 3D model of human epiblast (https://doi.org/10.1101/2023.12.07.570575)
* GEO Series: [GSE248567](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE248567)

### Overview

This repository contains analyses corresponding to different sequencing datasets and figures used in the study. Each subfolder includes description of the pre-proccessing steps (alignment, quantification, etc.), code, processed data, and results specific to a figure or dataset.

#### Folder Structure

* `public_datasets/`
Single-cell RNA-seq analysis of publicly available human and mouse embryo datasets.

* `RNAseq_ZNF398_CRISPRi/`
Analysis of RNA-seq data from ZNF398 CRISPRi experiments.

* `scRNAseq_3D-hE-Gastruloids_H9/`
Single-cell RNA-seq analysis of 3D-hE-Gastruloids derived from H9 cell line.

* `scRNAseq_3D-hE-Gastruloids_HES3-MIXL1-GFP/`
Single-cell RNA-seq analysis of 3D-hE-Gastruloids derived from HES3-MIXL1-GFP line.

* `scRNAseq_HumanEmbryo/`
Single-cell RNA-seq analysis of human embryo datasets for comparison and integration with 3D-hE-Gastruloids data.


#### Data availability

The processed human embryo reference scRNA-seq dataset was retrieved from: https://petropoulos-lanner-labs.clintec.ki.se/dataset.download.html. The mouse embryo scRNA-seq dataset was retrieved from the Gene Expression Omnibus (GEO) database, with accession code: GSE133725. The human embryo scRNA-seq dataset for comparative analyses was retrieved from the GEO database, with accession code: GSE136447. The ZNF398 ChIP-seq dataset was retrieved from the GEO database with accession code: GSE133630. The datasets generated in this study are available as raw and processed data in the GEO database with accession code GSE248567.
