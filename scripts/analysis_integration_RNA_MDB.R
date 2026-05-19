# scRNA-seq - Integration
# written by Maarten De Backere

#######################################
# Packages
#######################################
# Check if the required packages are installed - if not, install them

req_packages <- c("Seurat","ggplot2","harmony","viridis")

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
library(viridis)
set.seed(1234)

#######################################
# Load in the data
#######################################

# store location where the data, figures and scripts are located, without / on the end
base_dir <- "/scratch/svc_td_compbio/users/MaDeBa"

# load in PCA data
# PCA
# load(paste0(base_dir,"/data/seurat_objects_list_RNAonly_filtered_merged_normalized_HVG_joined_pca.RData"))

# load in non pca data
load(paste0(base_dir, "/data/seurat_objects_list_RNAonly_filtered_merged.RData"))

# 1. Normalize per-sample (operates on each counts.SAMPLE layer separately)
seurat_obj_merged <- NormalizeData(seurat_obj_merged)

# 2. Find variable features per-sample
seurat_obj_merged <- FindVariableFeatures(seurat_obj_merged)

# 3. Join layers (collapses into single counts + data matrix)
# seurat_obj_merged_joined <- JoinLayers(seurat_obj_merged)

# 4. Scale and continue
seurat_obj_merged <- ScaleData(seurat_obj_merged)

seurat_obj_merged[["percent.mt"]] <- PercentageFeatureSet(seurat_obj_merged, pattern = "^MT-")
seurat_obj_merged <- subset(
  seurat_obj_merged,
  subset = percent.mt < 10
)

# run pca
seurat_obj_merged <- RunPCA(seurat_obj_merged)

save(seurat_obj_merged, file= paste0(base_dir,"/data/seurat_objects_list_RNA_pca_nojoin.RData"))
load(paste0(base_dir,"/data/seurat_objects_list_RNA_pca_nojoin.RData"))
# seurat_obj_merged@meta.data$patient <- sub(".*-(Patient[0-9]+).*", "\\1", seurat_obj_merged@meta.data$orig.ident)


#######################################
# Integration
#######################################

seurat_obj_merged <- IntegrateLayers(
  object = seurat_obj_merged,
  method = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction = "harmony",
  group.by.vars = "orig.ident",
  theta = 2,
  verbose = TRUE
)

# 3. Join layers (collapses into single counts + data matrix)
seurat_obj_merged_joined <- JoinLayers(seurat_obj_merged)

# check the different layers 
head(Layers(seurat_obj_merged)) #before
head(Layers(seurat_obj_merged_joined)) #after joining
# colnames(seurat_obj_merged@meta.data)

# save integrated data (includes integration method)
save(seurat_obj_merged, paste0(base_dir,"/data/seurat_objects_list_RNA_pca_integrated_harmony_nojoin.RData"))