# scRNA-seq QC generator
# written by Maarten De Backere

#######################################
# Packages
#######################################
# Check if the required packages are installed - if not, install them

req_packages <- c("Seurat","ggplot2","openxlsx")

my_packages <- installed.packages()[, "Package"] #list installed packages
for (pkg in req_packages) {
  if(!(pkg %in% my_packages)) {
    install.packages(pkg)
  } 
}
# print("All packages installed.")

library(Seurat)
library(ggplot2)
library(openxlsx)
set.seed(1234)

#######################################
# QC
#######################################

# load("/scratch/svc_td_compbio/users/MaDeBa/data/seurat_objects_list_RNA_umap.RData")

# only if not done in advance
# create patient column with patient id, extracted from orig.ident
# seurat_obj_merged_joined@meta.data$patient <- sub(".*-(Patient[0-9]+).*", "\\1", seurat_obj_merged_joined@meta.data$orig.ident)

# explore metadata
# dim(seurat_obj_merged_joined@meta.data)
# colnames(metadata)
# head(metadata)
# 
# str(metadata)
# summary(metadata)

# copy metadata to create new dataframe
metadata <- seurat_obj_merged_joined@meta.data

# select all unique patients and sort them from 0 to 9
patient_un <- sort(unique(metadata$patient))

# prepare empty dataframe, with rows = patients and columns = different qc metrics
qc_df_patient <- data.frame(
  patient = patient_un,
  nCells_RNA = NA
)

#prepare counter for to be able to store the data in the right row (patient)
qc_row<- 1

qc_metrics <- c("nCount_RNA","nFeature_RNA","percent.mt")

# loop that goes over each patient (sorted before) and calculates the different qc metrics based on metadata that is subset for each metric and patient
for (pt in patient_un) {
  subset_mtd <- metadata[metadata$patient == pt,]
  qc_df_patient$nCells_RNA[qc_row]<- nrow(subset_mtd)
  
  for (metric in qc_metrics){ 
    sub_mtd_metric <- subset_mtd[[metric]]
    qc_df_patient[[paste0("avg_",metric)]][qc_row] <- mean(sub_mtd_metric)
    qc_df_patient[[paste0("med_",metric)]][qc_row] <- median(sub_mtd_metric)
    qc_df_patient[[paste0("min_",metric)]][qc_row] <- min(sub_mtd_metric)
    qc_df_patient[[paste0("max_",metric)]][qc_row] <- max(sub_mtd_metric)
  }
  
  qc_row <- qc_row+1
}

# View(qc_df_patient)


#######################################
# Export
#######################################
write.xlsx(
  x = qc_df_patient,
  file = "/scratch/svc_td_compbio/users/MaDeBa/figures/qc_df.xlsx",
  sheetName = "patient",
  rownames=FALSE
)

