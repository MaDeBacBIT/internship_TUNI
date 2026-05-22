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

# load in preprocessed and PCA data from analysis_base_RNA_MDB.R
pca_file <- paste0(base_dir,"/data/seurat_objects_list_RNA_pca_nojoin.RData")
load(pca_file)

#######################################
# Integration
#######################################

# run harmony integration
seurat_obj_merged <- IntegrateLayers(
  object = seurat_obj_merged,
  method = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction = "harmony",
  group.by.vars = "orig.ident",
  # theta = 4,
  # lambda = NULL,
  verbose = TRUE
)

# save integrated data (includes integration method)
save(seurat_obj_merged, file = paste0(base_dir,"/data/seurat_objects_list_RNA_pca_integrated_harmony_unjoin.RData"))