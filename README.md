# Internship Tampere University
## Integrating and clustering scRNAseq data to reveal cell types in primary prostate cancer tumors
Internship project at Tampere University. The scripts are used for integrating and clustering scRNAseq data to reveal cell types in primary prostate cancer tumors.
The dataset used is a 10x RNA and ATAC single cell multiomics dataset from primary prostate cancer and metastasis samples, published by Keshavarzian et al. in Nature communications. DOI: [https://doi.org/10.1038/s41467-025-67856-5](https://doi.org/10.1038/s41467-025-67856-5)

After proving that the batch effect is visible, one or more integration methods will be tested.

The data, figures and scripts are all stored on the TCSC server and afterwards synced to GitHub. The data is stored in a Seurat object and saved to an .Rdata file format. Figures are stored in .pdf and scripts in .R.

The structure of the repository is as following:
```
├── README.md
└── scripts
    ├── analysis_base_RNA_MDB.R
    ├── analysis_cluster_RNA_MDB.R
    ├── analysis_integration_RNA_MDB.R
    ├── analysis_RNA_stacked_barplot_cluster.R
    └── qc_RNA_MDB.R
```
## Scripts
### Preprocessing - script not included
In advance, preprocessing of the dataset was already done. The data was loaded in, filtered for RNA and minimal amount of reads. The different Seurat objects of the samples were merged together and normalized. Next, the variable features per sample were searched, layers joined and the object was scaled. The data was also subset to exclude cells with high mitochondrial read percentages. 

### PCA script `/scripts/analysis_base_RNA_MDB.R`
**Features:** PCA

First, the required packages are installed and loaded. Next, the data is explored and dimensions are reduced with PCA. Requires preprocessed data.

### Integration `/scripts/analysis_integration_RNA_MDB.R`
**Features:** Harmony integration / join layers

After PCA is ran, the integration can be done. Using IntegrateLayers(), Harmony integration is used. Temporarily this script also includes the preprocessing and PCA, but this will be moved in the future. Will require PCA and preprocessing to be done in advance.

### Clustering and visualization `/scripts/analysis_cluster_RNA_MDB.R`
**Features:** KNN / cluster / UMAP / CellCycleScores / Marker expression

Neighbors and clusters are computed. The Seurat clusters are visualized using UMAP and colored/grouped based on the different metadata columns. CellCycleScores are calculated and the different phases are visualized. Using celltype markers, marker expression was shown on the FeaturePlot and by using a DotPlot. Can use the data with and without integration.

### Stacked barplot `/scripts/analysis_RNA_stacked_barplot_cluster.R`
**Features:** stacked barplot

Visualizes cluster composition across samples/patients using a stacked barplot. Requires clustering to be done in advance.

### Calculating QC metrics `/scripts/qc_RNA_MDB.R`
**Features:** QC metrics / export to .xlsx workbook

Using the metadata, a new row is created for each patient/sample with columns containing the amount of cells and the average / median / minimum / maximum of reads/cell, features/cell and percent mitochondrial features. The data is exported to an .xlsx workbook, with a separate worksheet for patients and samples.

## Dependencies
R packages: `Seurat`, `ggplot2`, `harmony`, `readxl`, `openxlsx`

## Author
Scripts created by Maarten De Backere for Tampere University (2026)
