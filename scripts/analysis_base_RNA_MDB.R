# scRNA-seq - PCA
# written by Maarten De Backere

#######################################
# Packages
#######################################
# Check if the required packages are installed - if not, install them

req_packages <- c("Seurat","ggplot2")

my_packages <- installed.packages()[, "Package"] #list installed packages
for (pkg in req_packages) {
  if(!(pkg %in% my_packages)) {
    install.packages(pkg)
  } 
}
# load the different packages
library(Seurat)
library(ggplot2)

# set a seed
set.seed(1234)

print("All packages loaded and installed.")

#######################################
# Load in / create preprocessed data
#######################################

# store location where the data, figures and scripts are located, without / on the end
base_dir <- "/scratch/svc_td_compbio/users/MaDeBa"

# preprocessing as layers cannot be joined in advance
# load in non pca data
load(paste0(base_dir, "/data/seurat_objects_list_RNAonly_filtered_merged.RData"))

# Normalize per-sample (operates on each counts.SAMPLE layer separately)
seurat_obj_merged <- NormalizeData(seurat_obj_merged)

# Find variable features per-sample
seurat_obj_merged <- FindVariableFeatures(seurat_obj_merged)

# Scale and continue
seurat_obj_merged <- ScaleData(seurat_obj_merged)

# calculate percent mt and subset to have it lower than 10
seurat_obj_merged[["percent.mt"]] <- PercentageFeatureSet(seurat_obj_merged, pattern = "^MT-")
seurat_obj_merged <- subset(seurat_obj_merged, subset = percent.mt < 10)

# create patient column with patient id, extracted from orig.ident
seurat_obj_merged@meta.data$patient <- sub(".*-(Patient[0-9]+).*", "\\1", seurat_obj_merged@meta.data$orig.ident)

#######################################
# PCA
#######################################

# ?RunPCA
# Check size of seurat object
# object.size(seurat_obj_merged)

# Run PCA (linear dimensionality reduction)
seurat_obj_merged <- RunPCA(seurat_obj_merged)

# Save and join 
pca_file <- paste0(base_dir,"/data/seurat_objects_list_RNA_pca_nojoin.RData")
save(seurat_obj_merged, file = pca_file)
load(pca_file)

# Plots stored in pdf format in /figures/explore_seurat/...

# Create scatter plot of PC1 and PC2 (colored by patient)
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

# Create scatter plot of PC1 and PC27 (colored by patient)
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

print("Plots created and base script finished.")
