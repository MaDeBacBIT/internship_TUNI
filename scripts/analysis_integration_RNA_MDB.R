# scRNA-seq - Integration
# written by Maarten De Backere

#######################################
# Packages
#######################################
# Check if the required packages are installed - if not, install them

req_packages <- c("Seurat","ggplot2","harmony","reticulate","devtools")

my_packages <- installed.packages()[, "Package"] #list installed packages
for (pkg in req_packages) {
  if(!(pkg %in% my_packages)) {
    install.packages(pkg)
  } 
}
# load in reticulate for connection with commandline
library (reticulate)

# select python version for reticulate
use_python("/home/pdr436/miniforge3/envs/scvi-env-310/bin/python", required = TRUE)

# show information of python version (used by reticulate)
# py_config()

# check if seuratwrappers is installed, install if needed
if (!("SeuratWrappers" %in% my_packages)) {
  library(devtools)
  install_github("satijalab/seurat-wrappers")
}

# load the different packages
library(Seurat)
library(ggplot2)
library(harmony)
library(SeuratWrappers)

# set a seed
set.seed(1234)

print("All packages loaded and installed.")

#######################################
# Load in the data
#######################################

# store location where the data, figures and scripts are located, without / on the end
base_dir <- "/scratch/svc_td_compbio/users/MaDeBa"

# load in preprocessed and PCA data from analysis_base_RNA_MDB.R
pca_file <- paste0(base_dir,"/data/seurat_objects_list_RNA_pca_nojoin.RData")
load(pca_file)
print("Seurat object with PCA reduction loaded.")

#######################################
# Integration
#######################################
# choose integration method: "harmony" or "scVI" => also used for object file name
integr_method <- "scVI"
print(paste0("Integration method: ", integr_method))

# check what integration method should be ran
if (integr_method == "harmony") {
  # run harmony integration
  seurat_obj_merged <- IntegrateLayers(
    object = seurat_obj_merged,
    method = HarmonyIntegration,
    orig.reduction = "pca",
    new.reduction = "harmony",
    group.by.vars = "orig.ident",
    # theta = 1,
    # lambda = NULL,
    verbose = TRUE
  )
}
# ?scVIIntegration
if (integr_method == "scVI") {
  # run scVI integration
  seurat_obj_merged <- IntegrateLayers(
    object = seurat_obj_merged,
    method = scVIIntegration,
    orig.reduction = "pca",
    new.reduction = "integrated.scvi",
    # specify conda environment location
    conda_env = "/home/pdr436/miniforge3/envs/scvi-env-310",
    verbose = TRUE
  )
}
print(paste0("Integration (",integr_method,") finished"))
#######################################
# Save data
#######################################
# save integrated data (includes integration method)
integr_file <- paste0(base_dir,"/data/seurat_objects_list_RNA_pca_integrated_",integr_method,"_unjoin.RData")
save(seurat_obj_merged, file = integr_file)

print(paste0("Saved Seurat object to: ",integr_file))
