# scRNA-seq - PCA
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

# save seurat object after pca
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
