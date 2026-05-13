# scRNA-seq
# written by Maarten De Backere

#######################################
# Packages
#######################################
# Check if the required packages are installed - if not, install them

req_packages <- c("Seurat","ggplot2","harmony","openxlsx","viridis")

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
library(readxl)
library(openxlsx)
library(viridis)
set.seed(1234)

#######################################
# load in the preprocessed data
#######################################

# store location where the data, figures and scripts are located, without / on the end
base_dir <- "/scratch/svc_td_compbio/users/MaDeBa"


# load("/scratch/svc_td_compbio/users/MaDeBa/data/seurat_objects_list_RNAonly_filtered_merged_normalized_HVG_joined.RData")

# have a look at the seurat object
seurat_obj_merged_joined

# ?save()

# save the seurat object to a different file and then continue to use the new file
# save(seurat_obj_merged_joined, file ="/scratch/svc_td_compbio/users/MaDeBa/data/seurat_objects_list_RNAonly_filtered_merged_normalized_HVG_joined_copy_analysis.RData")

# load in the copy of preprocessed data, this will be overwritten further in the analysis
load("/scratch/svc_td_compbio/users/MaDeBa/data/seurat_objects_list_RNAonly_filtered_merged_normalized_HVG_joined_copy_analysis.RData")

# table(seurat_obj_merged_joined$orig.ident) #cell counts?
# head(seurat_obj_merged_joined[[]], 2)
# seurat_obj_merged_joined@
# seurat_obj_merged_joined@meta.data
# head(colnames(seurat_obj_merged_joined@meta.data))

# seurat_obj_merged_joined@assays
# Assays(seurat_obj_merged_joined)
# VariableFeatures(seurat_obj_merged_joined)[1:20]
# DefaultAssay(seurat_obj_merged_joined)

# create patient column with patient id, extracted from orig.ident
# seurat_obj_merged_joined@meta.data$patient <- sub(".*-(Patient[0-9]+).*", "\\1", seurat_obj_merged_joined@meta.data$orig.ident)

#######################################
# QC => scripts/qc_RNA_MDB.R
#######################################

#######################################
# PCA
#######################################

# ?RunPCA
# Check size of seurat object
# object.size(seurat_obj_merged_joined)

# Run PCA (linear dimensionality reduction)
seurat_obj_merged_joined <- RunPCA(seurat_obj_merged_joined)

# save data
# save(seurat_obj_merged_joined, file = "/scratch/svc_td_compbio/users/MaDeBa/data/seurat_objects_list_RNAonly_filtered_merged_normalized_HVG_joined_pca.RData")

# load in saved data
load("/scratch/svc_td_compbio/users/MaDeBa/data/seurat_objects_list_RNAonly_filtered_merged_normalized_HVG_joined_pca.RData")


# Plots stored in pdf format in /figures/explore_seurat/...

# Create scatter plot of PC1 and PC2 
# pdf(file = "/scratch/svc_td_compbio/users/MaDeBa/figures/explore_seurat/rna_PCA_PC1-2_plot.pdf")
# DimPlot(seurat_obj_merged_joined, reduction = "pca") + 
#   NoLegend() + 
#   labs(title = "PCA plot") +
#   theme(plot.title = element_text(hjust=0.5))
# dev.off()
# option to color bypatient: group by

png(filename = "/scratch/svc_td_compbio/users/MaDeBa/figures/explore_seurat/rna_PCA_PC1-2_plot_patient_legend.png",
    width = 2000,
    height = 2000,
    res = 300,
    type = "cairo")
DimPlot(seurat_obj_merged_joined,
        dims = c(1,2),
        reduction = "pca",
        group.by = "patient") +
  # NoLegend() +
  labs(title = "PCA plot 1-2") +
  theme(plot.title = element_text(hjust=0.5),
        legend.position = "bottom")
dev.off()

png(filename = "/scratch/svc_td_compbio/users/MaDeBa/figures/explore_seurat/rna_PCA_PC1-27_plot_patient.png",
    width = 2000,
    height = 2000,
    res = 300,
    type = "cairo")
DimPlot(seurat_obj_merged_joined,
        dims = c(1,27),
        reduction = "pca",
        group.by = "patient") +
  # NoLegend() +
  labs(title = "PCA plot 1-27") +
  theme(plot.title = element_text(hjust=0.5),
        legend.position = "bottom")
dev.off()


# Create elbowplot to see what PC capture most data variance
# ?ElbowPlot
pdf(file = "/scratch/svc_td_compbio/users/MaDeBa/figures/explore_seurat/rna_PCA_elbowplot_40dim.pdf", width = 10, height = 6)
ElbowPlot(seurat_obj_merged_joined, ndims = 40)+
  geom_line(aes(x=dims,y=y_data)) +
  labs(title = "Elbow plot of PCA variance") +
  theme(plot.title = element_text(hjust=0.5))
dev.off()

# Visualize top genes associated with reduction components
# pdf(file = "/scratch/svc_td_compbio/users/MaDeBa/figures/explore_seurat/rna_PCA_vizdim1-2.pdf")
# VizDimLoadings(seurat_obj_merged_joined, dims = 1:2, reduction = "pca")
# dev.off()

# Draws a heatmap focusing on a principal component (displays top 15 genes with highest and lowest pc scores )
pdf(file = "/scratch/svc_td_compbio/users/MaDeBa/figures/explore_seurat/rna_PCA_heatmap1-15.pdf", width=30, height=30)
DimHeatmap(seurat_obj_merged_joined,
           dims = 1:15,
           cells = 5000,
           balanced = TRUE,
           reduction = "pca")
dev.off()

#######################################
# Integration
#######################################

# Integration will be located here in the analysis


# -------------------------- build in save method here -----------------

#######################################
# Neighbors/clustering
#######################################

# Find nearest neighbors, within pca dims 1:30 
# seurat_obj_merged_joined <- FindNeighbors(seurat_obj_merged_joined, dims = 1:30)
seurat_obj_merged_joined <- FindNeighbors(seurat_obj_merged_joined, dims = 1:27)

# Find clusters
# seurat_obj_merged_joined <- FindClusters(seurat_obj_merged_joined, resolution = 0.5)
seurat_obj_merged_joined <- FindClusters(seurat_obj_merged_joined, resolution = 0.6)


# Check what cluster the first 5 cells belong to
# head(Idents(seurat_obj_merged_joined), 5)

#######################################
# UMAP Visualization
#######################################

# Run umap reduction method, with PCA dimensions 1:30
# seurat_obj_merged_joined <- RunUMAP(seurat_obj_merged_joined, dims = 1:30)
seurat_obj_merged_joined <- RunUMAP(seurat_obj_merged_joined, dims = 1:27)

# Visualize clusters (to pdf file)
# pdf(file = "/scratch/svc_td_compbio/users/MaDeBa/figures/explore_seurat/rna_umap_v1.pdf", width=12, height=9)
# DimPlot(seurat_obj_merged_joined, reduction = "umap")
# dev.off()

# check possible groupings on according to metadata
# colnames(seurat_obj_merged_joined[[]])

# group by sample name
# DimPlot(seurat_obj_merged_joined, reduction = "umap", group.by = "orig.ident") + 
#   NoLegend()

# # Visualize clusters (to pdf file) grouped by cluster, with label and higher resolution
# pdf(file = "/scratch/svc_td_compbio/users/MaDeBa/figures/explore_seurat/rna_umap_1-27dim_06res_cluster.pdf",
#     width=12,
#     height=9)
# DimPlot(seurat_obj_merged_joined,
#         reduction = "umap",
#         group.by = "seurat_clusters",
#         label = TRUE,
#         raster.dpi = c(600,600)) + 
#   NoLegend() +
#   labs(title = "RNA - cluster (PCA dims = 1:27, res = 0.6)") +
#   theme(
#     plot.title = element_text(hjust=0.5)
#   )
# dev.off()
# 
# ## BY PATIENT
# 
# pdf(file = "/scratch/svc_td_compbio/users/MaDeBa/figures/explore_seurat/rna_umap_1-27dim_06res_patient.pdf",
#     width=12,
#     height=9)
# DimPlot(seurat_obj_merged_joined,
#         reduction = "umap",
#         group.by = "patient",
#         label = FALSE,
#         raster.dpi = c(600,600)) +
#   labs(title = "RNA - patient (PCA dims = 1:27, res = 0.6)") +
#   theme(
#     plot.title = element_text(hjust=0.5),
#     legend.position = "bottom"
#         )
# dev.off()
# 
# ## BY SAMPLE
# pdf(file = "/scratch/svc_td_compbio/users/MaDeBa/figures/explore_seurat/rna_umap_1-27dim_06res_sample.pdf",
#     width=18,
#     height=9)
# DimPlot(seurat_obj_merged_joined,
#         reduction = "umap",
#         group.by = "orig.ident",
#         label = FALSE,
#         raster.dpi = c(900,900)) +
#   labs(title = "RNA - sample (PCA dims = 1:27, res = 0.6)") +
#   theme(
#     plot.title = element_text(hjust=0.5),
#     legend.position = "right"
#   )
# dev.off()
# 
# ## BY ncount
# pdf(file = "/scratch/svc_td_compbio/users/MaDeBa/figures/explore_seurat/rna_umap_1-27dim_06res_ncount.pdf",
#     width=12,
#     height=9)
# FeaturePlot(seurat_obj_merged_joined,
#         reduction = "umap",
#         features = "nCount_RNA") +
#   # NoLegend() +
#   labs(title = "RNA - nCount (PCA dims = 1:27, res = 0.6)") +
#   theme(
#     plot.title = element_text(hjust=0.5)
#   )
# dev.off()
# 
# ## BY nfeature
# pdf(file = "/scratch/svc_td_compbio/users/MaDeBa/figures/explore_seurat/rna_umap_1-27dim_06res_nfeature.pdf",
#     width=12,
#     height=9)
# FeaturePlot(seurat_obj_merged_joined,
#             reduction = "umap",
#             features = "nFeature_RNA") +
#   # NoLegend() +
#   labs(title = "RNA - nFeature (PCA dims = 1:27, res = 0.6)") +
#   theme(
#     plot.title = element_text(hjust=0.5)
#   )
# dev.off()
# 
# ## BY percentmt
# pdf(file = "/scratch/svc_td_compbio/users/MaDeBa/figures/explore_seurat/rna_umap_1-27dim_06res_percentmt.pdf",
#     width=12,
#     height=9)
# FeaturePlot(seurat_obj_merged_joined,
#             reduction = "umap",
#             features = "percent.mt") +
#   # NoLegend() +
#   labs(title = "RNA - percent.mt (PCA dims = 1:27, res = 0.6)") +
#   theme(
#     plot.title = element_text(hjust=0.5)
#   )
# dev.off()


################################## CREATE ALL UMAP plots with 1 loop

# # create subdirectory (not ending on /), chooses what to group by (also renamed) , visualize marker expression on umap, use color viridis and improve for all celltypes
subdir_umap <- "/figures/rna_umap"
umap_groups <- list(
  "cluster" = "seurat_clusters",
  "patient" = "patient",
  "sample"= "orig.ident",
  "nCount" = "nCount_RNA",
  "nFeature" = "nFeature_RNA",
  "percent_mt" = "percent.mt"
  )

# go over list of group names
for (grp in names(umap_groups)) {
  # select groups
  orig_grp <- umap_groups[[grp]]

  # create pdf with renamed group in name, change width if sample (because of high amount of samples)
  if (grp == "sample") {
    pdf_umap_width <- 18
  } else {pdf_umap_width <- 12}
  pdf(file = paste0(base_dir,subdir_umap,"/rna_umap_",grp,".pdf"),
      width=pdf_umap_width,
      height=9)
  
  # create DimPlot for categorical variables and featureplot for numerical
  if (grp %in% c("cluster","patient","sample")) {
    umap_p <- DimPlot(seurat_obj_merged_joined,
                        reduction = "umap",
                        group.by = orig_grp,
                        label = FALSE,
                        raster.dpi = c(600,600))
  } else {
      umap_p <- FeaturePlot(seurat_obj_merged_joined,
                            reduction = "umap",
                            features = orig_grp,
                            cols = viridis(100))
  }
  
  # apply styling and add title
  umap_p <- umap_p +
    labs(title = paste0("RNA - ",grp," (PCA dims = 1:27, res = 0.6)")) +
    theme(
      plot.title = element_text(hjust=0.5)
      # ,legend.position = "bottom"
        )
  
  # print plot and close pdf
  print(umap_p)
  dev.off()
  # to check what has already been done: (prints in console)
  print(grp)
    }

#Cluster composition per sample
table(seurat_obj_merged_joined$seurat_clusters, seurat_obj_merged_joined$orig.ident)

# save data
# save(seurat_obj_merged_joined, file = "/scratch/svc_td_compbio/users/MaDeBa/data/seurat_objects_list_RNA_umap.RData")

# load in saved data
load("/scratch/svc_td_compbio/users/MaDeBa/data/seurat_objects_list_RNA_umap.RData")


#######################################
# CellCycleScoring
#######################################

seurat_obj_merged_joined <- CellCycleScoring(seurat_obj_merged_joined, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes)

## BY PHASE

pdf(file = "/scratch/svc_td_compbio/users/MaDeBa/figures/explore_seurat/rna_CCS_umap_phases.pdf",
    width = 12, height = 9)
DimPlot(seurat_obj_merged_joined,
        reduction = "umap",
        group.by = "Phase") +
  labs(title = "RNA - phases (PCA dims = 1:27, res = 0.6)") +
  theme(
    plot.title = element_text(hjust=0.5),
    legend.position = "bottom"
  )
dev.off()

pdf(file = "/scratch/svc_td_compbio/users/MaDeBa/figures/explore_seurat/rna_CCS_umap_split_phases.pdf",
    width = 15, height = 9)
DimPlot(seurat_obj_merged_joined,
        reduction = "umap",
        split.by = "Phase",
        group.by = "Phase") +
  labs(title = "RNA - phases (PCA dims = 1:27, res = 0.6)") +
  theme(
    plot.title = element_text(hjust=0.5),
    legend.position = "bottom"
  )
dev.off()

pdf(file = "/scratch/svc_td_compbio/users/MaDeBa/figures/explore_seurat/rna_CCS_umap_phase-cluster.pdf",
    width = 15, height = 9)
DimPlot(seurat_obj_merged_joined,
        reduction = "umap",
        split.by = "Phase",
        group.by = "seurat_clusters",
        label = TRUE) +
  labs(title = "RNA - phases + clusters (PCA dims = 1:27, res = 0.6)") +
  theme(
    plot.title = element_text(hjust=0.5),
    legend.position = "none"
  )
dev.off()

# pdf(file = "/scratch/svc_td_compbio/users/MaDeBa/figures/explore_seurat/rna_CCS_vln_scores.pdf",
#     width = 15, height = 9)
# VlnPlot(seurat_obj_merged_joined,
#         features = c("S.Score", "G2M.Score"),
#         group.by = "orig.ident")
# dev.off()

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

# visualize marker expression on umap, use color viridis and improve for all celltypes
subdir_ftplt <- "/figures/rna_marker_FeatPlot"
for (ct in names(markers)) {
  # select celltype markers for each celltype
  ct_markers <- markers[[ct]]
  
  # create pdf with celltype in name
  pdf(file = paste0(base_dir,subdir_ftplt,"/rna_umap_ftplt_",ct,".pdf"),
      width=12,
      height=9)
  # create feature plot with colors of viridis, for each celltype
  ft_plt <- FeaturePlot(   
    seurat_obj_merged_joined,   
    features = ct_markers, 
    cols = viridis(100) 
    )
  
  # add plot to pdf and close pdf
  print(ft_plt)
  dev.off()
} 


# use dotplot to show expression, change colors and y-axis
pdf(file = "/scratch/svc_td_compbio/users/MaDeBa/figures/explore_seurat/rna_markers_dotplot_cluster.pdf", width = 22, height = 10)
DotPlot(seurat_obj_merged_joined,
        features = markers,
        group.by = "seurat_clusters") + 
  scale_color_viridis_b() +
  labs(title = "Marker expression - clusters",
       x = "Cell type - gene markers",
       y = "Clusters") +
  theme(plot.title = element_text(hjust=0.5))+
  RotatedAxis()
dev.off()

