# scRNA-seq
# written by Maarten De Backere

#######################################
# Packages
#######################################
# Check if the required packages are installed - if not, install them

req_packages <- c("Seurat","ggplot2","harmony")

my_packages <- installed.packages()[, "Package"] #list installed packages
for (pkg in req_packages) {
  if(!(pkg %in% my_packages)) {
    install.packages(pkg)
  } 
}
# print("All packages installed.")

library(Seurat)
library(ggplot2)
library(harmony)
set.seed(1234)

#######################################
# load in the preprocessed data
#######################################
load("/scratch/svc_td_compbio/users/MaDeBa/data/seurat_objects_list_RNAonly_filtered_merged_normalized_HVG_joined.RData")

seurat_obj_merged_joined
# table(seurat_obj_merged_joined$orig.ident) #cell counts?
# head(seurat_obj_merged_joined[[]], 2)
# seurat_obj_merged_joined@
# seurat_obj_merged_joined@meta.data
# head(colnames(seurat_obj_merged_joined@meta.data))

# seurat_obj_merged_joined@assays
# Assays(seurat_obj_merged_joined)
VariableFeatures(seurat_obj_merged_joined)[1:20]

#######################################
# PCA
#######################################

# ?RunPCA
# Check size of seurat object
object.size(seurat_obj_merged_joined)

# Run PCA (linear dimensionality reduction)
seurat_obj_merged_joined <- RunPCA(seurat_obj_merged_joined)

# Plots stored in pdf format in /figures/explore_seurat/...

# Create scatter plot of PC1 and PC2 
# pdf(file = "/scratch/svc_td_compbio/users/MaDeBa/figures/explore_seurat/rna_PCA_PC1-2_plot.pdf")
# DimPlot(seurat_obj_merged_joined, reduction = "pca") + 
#   NoLegend() + 
#   labs(title = "PCA plot") +
#   theme(plot.title = element_text(hjust=0.5))
# dev.off()

# Create elbowplot to see what PC capture most data variance
pdf(file = "/scratch/svc_td_compbio/users/MaDeBa/figures/explore_seurat/rna_PCA_elbowplot.pdf")
ElbowPlot(seurat_obj_merged_joined)+
  geom_line(aes(x=dims,y=y_data)) +
  labs(title = "Elbow plot of PCA variance") +
  theme(plot.title = element_text(hjust=0.5))
dev.off()

# Visualize top genes associated with reduction components
# pdf(file = "/scratch/svc_td_compbio/users/MaDeBa/figures/explore_seurat/rna_PCA_vizdim1-2.pdf")
# VizDimLoadings(seurat_obj_merged_joined, dims = 1:2, reduction = "pca")
# dev.off()

# Draws a heatmap focusing on a principal component (displays top 15 genes with highest and lowest pc scores )
pdf(file = "/scratch/svc_td_compbio/users/MaDeBa/figures/explore_seurat/rna_PCA_heatmap1-15.pdf", width=10, height=10)
DimHeatmap(seurat_obj_merged_joined, dims = 1:9, cells = 5000, balanced = TRUE, reduction = "pca")
dev.off()

# -------------------------- build in save method here -----------------

#######################################
# Integration
#######################################

# Integration will be located here in the analysis

# -------------------------- build in save method here -----------------

#######################################
# Neighbors/clustering
#######################################

# Find nearest neighbors, within pca dims 1:30 
seurat_obj_merged_joined <- FindNeighbors(seurat_obj_merged_joined, dims = 1:30)

# Find clusters
seurat_obj_merged_joined <- FindClusters(seurat_obj_merged_joined, resolution = 0.5)

# Check what cluster the first 5 cells belong to
# head(Idents(seurat_obj_merged_joined), 5)

#######################################
# UMAP Visualization
#######################################

# Run umap, with PCA dimensions 1:30
seurat_obj_merged_joined <- RunUMAP(seurat_obj_merged_joined, dims = 1:30)

# Visualize clusters (to pdf file)
pdf(file = "/scratch/svc_td_compbio/users/MaDeBa/figures/explore_seurat/rna_umap_v1.pdf", width=12, height=9)
DimPlot(seurat_obj_merged_joined, reduction = "umap")
dev.off()

# check possible groupings on according to metadata
# colnames(seurat_obj_merged_joined[[]])

# group by sample name
# DimPlot(seurat_obj_merged_joined, reduction = "umap", group.by = "orig.ident") + 
#   NoLegend()

# Visualize clusters (to pdf file) grouped by cluster, with label and higher resolution
pdf(file = "/scratch/svc_td_compbio/users/MaDeBa/figures/explore_seurat/rna_umap_v2.pdf", width=12, height=9)
DimPlot(seurat_obj_merged_joined, reduction = "umap", group.by = "seurat_clusters", label = TRUE, raster.dpi = c(900,900)) + 
  NoLegend()
dev.off()


#Cluster composition per sample
table(seurat_obj_merged_joined$seurat_clusters, seurat_obj_merged_joined$orig.ident)

# -------------------------- build in save method here -----------------

#######################################
# Markers
#######################################

# load in markers (provided)
markers <- list(
  Tumor = c("AMACR", "NKX3-1", "FOLH1", "ERG", "SCHLAP1", "SPINK1",
            "PCA3"),
  T_cells = c("CD3D", "CD3E", "CD8A", "CD4", "FOXP3", "TRAC"),
  Luminal = c("KRT8", "KRT18", "KLK3", "AR", "ACPP", "HPN", "MSMB"),
  Stroma = c("DCN", "LUM", "COL1A1", "COL3A1", "VIM", "PDGFRA", "FAP",
             "ACTA2"),
  Basal = c("KRT5", "KRT14", "TP63", "CD44", "KRT15"),
  Myeloid = c("CD68", "LYZ", "AIF1", "CD14", "CSF1R", "MRC1"),
  B_cells = c("CD79A", "CD79B", "MS4A1", "CD19", "BANK1", "PAX5"),
  Club = c("SCGB1A1", "SCGB3A1", "KRT7", "PIGR"),
  Mast = c("CPA3", "TPSAB1", "TPSB2", "KIT", "MS4A2"))

# visualize marker expression on umap
pdf(file = "/scratch/svc_td_compbio/users/MaDeBa/figures/explore_seurat/rna_umap_luminal.pdf", width=12, height=9)
FeaturePlot(
  seurat_obj_merged_joined,
  features = markers$Luminal
)
dev.off()

#find markers
markers2 <- FindAllMarkers(seurat_obj_merged_joined)

