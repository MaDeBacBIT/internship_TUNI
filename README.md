# Internship Tampere University
## Integrating and clustering scRNAseq data to reveal cell types in primary prostate cancer tumors
Internship project at Tampere University. The scripts are used for integrating and clustering scRNAseq data to reveal cell types in primary prostate cancer tumors.
The dataset used is a 10x RNA and ATAC single cell multiomics dataset from primary prostate cancer and metastasis samples, published by Keshavarzian et al. in Nature communications. DOI: [https://doi.org/10.1038/s41467-025-67856-5](https://doi.org/10.1038/s41467-025-67856-5)

After proving that the batch effect is visible, multiple integration methods (Harmony and scVI) were tested.

The data, figures and scripts are all stored on the TCSC server and afterwards the scripts are synced to GitHub. The data is stored as a Seurat object and saved to an .Rdata file format. Figures are stored in .pdf and scripts in .R format.

The structure of the repository is as following:
```
├── README.md
└── scripts
    ├── analysis_base_RNA_MDB.R
    ├── analysis_cluster_RNA_MDB.R
    ├── analysis_integration_RNA_MDB.R
    ├── analysis_RNA_stacked_barplot_cluster.R
    ├── qc_RNA_MDB.R
    └── run_integration.sh
```
## Scripts
### Preprocessing - script not included

In advance, partial preprocessing of the dataset was already done. The data was loaded in, filtered for RNA and minimal amount of reads. The different Seurat objects of the samples were merged together.

### Preprocessing and PCA script `/scripts/analysis_base_RNA_MDB.R`
**Features:** Preprocessing, PCA

First, the required packages are installed and loaded. The Seurat object is further preprocessed. It is normalized and the variable features per sample are searched, layers joined and the object are scaled. The data is also subset to exclude cells with high mitochondrial read percentages. Next, the data is explored and dimensions are reduced with PCA.
Requires filtered and merged data.

### Integration `/scripts/analysis_integration_RNA_MDB.R`
**Features:** Harmony integration / scVI integration

After PCA is ran, the integration can be done. Using IntegrateLayers(), the selected method Harmony or scVI integration is used. For running the scVIIntegration, the `run_integration.sh` script should be used and a conda environment with scvi-tools is required. 
Requires preprocessing and PCA to be done in advance.

### Integration SLURM script `/scripts/run_integration.sh`
**Features:** create job / activate environment / library overrides

When started, a sbatch job is created. The terminal output is redirected to a log file. R and a conda environment with scvi-tools installed are loaded. For running the integration Rscript, library overrides are needed. 
Requires R and conda with scvi-tools to be installed and integration Rscript. 

### Clustering and visualization `/scripts/analysis_cluster_RNA_MDB.R`
**Features:** Join layers / KNN / cluster / UMAP / CellCycleScores / Marker expression

The different layers are joined. Neighbors and clusters are computed. The Seurat clusters are visualized using UMAP and colored/grouped based on the different metadata columns. CellCycleScores are calculated and the different phases are visualized. Using celltype markers, marker expression is shown on the FeaturePlot and a DotPlot. Clustree can be used to determine the resolution.

Can use the data with and without integration.

### Stacked barplot `/scripts/analysis_RNA_stacked_barplot_cluster.R`
**Features:** stacked barplot

Visualizes cluster composition across samples/patients using a stacked barplot.
Requires clustering to be done in advance.

### Calculating QC metrics `/scripts/qc_RNA_MDB.R`
**Features:** QC metrics / export to .xlsx workbook

Using the metadata, a new row is created for each patient/sample with columns containing the amount of cells and the average / median / minimum / maximum of reads/cell, features/cell and percent mitochondrial features. The data is exported to an .xlsx workbook, with a separate worksheet for patients and samples.

## Dependencies
R packages: [`Seurat`](https://satijalab.org/seurat/), [`ggplot2`](https://ggplot2.tidyverse.org/), [`Harmony`](https://github.com/immunogenomics/harmony), [`reticulate`](https://rstudio.github.io/reticulate/), [`SeuratWrappers`](https://satijalab-seurat-wrappers.mintlify.app/), [`viridis`](https://cran.r-project.org/web/packages/viridis/vignettes/intro-to-viridis.html), [`clustree`](https://cran.r-project.org/web/packages/clustree/refman/clustree.html), [`openxlsx`](https://cran.r-project.org/web/packages/openxlsx/refman/openxlsx.html)

Conda: [`scvi-tools`](https://docs.scvi-tools.org/en/stable/installation.html)

## Author
Scripts created by Maarten De Backere for Tampere University (2026)
