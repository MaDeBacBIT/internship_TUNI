# scRNA-seq - clustering and UMAP
# written by Maarten De Backere

#######################################
# Packages
#######################################
# Check if the required packages are installed - if not, install them
req_packages <- c("Seurat","ggplot2","viridis")

my_packages <- installed.packages()[, "Package"] #list installed packages
for (pkg in req_packages) {
  if(!(pkg %in% my_packages)) {
    install.packages(pkg)
  } 
}

# load the different packages
library(Seurat)
library(ggplot2)
library(viridis)
set.seed(1234)

print("All required packages loaded and installed.")

#######################################
# Load in the data
#######################################

# store location where the data, figures and scripts are located, without / on the end
base_dir <- "/scratch/svc_td_compbio/users/MaDeBa"
# dir.exists(base_dir) #check if file location exists

# choose integration method: "harmony", "scVI" or "none"
integr_method <- "scVI"

reduct_integr <- switch(integr_method,
                        "none" = "pca",
                        "harmony" = "harmony",
                        "scVI" = "integrated.scvi"
                        )
# load in data
if (integr_method == "none"){
  # PCA only (no integration)
  pca_file <- paste0(base_dir,"/data/seurat_objects_list_RNA_pca_nojoin.RData")
  load(pca_file)
} else {
  # integrated (depends on integration method)
  integr_file <- paste0(base_dir,"/data/seurat_objects_list_RNA_pca_integrated_",integr_method,"_unjoin.RData")
  load(integr_file)
}
print("Integrated/PCA Seurat object loaded.") # to confirm loading is done

#######################################
# Join the different layers
#######################################
# Join layers (collapses into single counts + data matrix)
seurat_obj_merged_joined <- JoinLayers(seurat_obj_merged)
print("Layers joined.")

# check the different layers 
print("Layers before joining:")
head(Layers(seurat_obj_merged)) #before
print("Layers after joining:")
head(Layers(seurat_obj_merged_joined)) #after joining

#######################################
# Specify parameters
#######################################

# choose pca dimensions and resolution for clustering
dims_pca <- 1:27
res_clust <- 0.6

# With clustree, multiple resolutions can be tested, set to True if this is wanted and a clustree will be created
clustree_check <- FALSE

# only specified when harmony integration and if needed. // if wanted, it can also be added to the title of the plots by uncommenting the #labs parts
# if (integr_method == "harmony"){
#   theta_harmony <- 1
#   }

#######################################
# Find neighbors + clusters
#######################################
# Find nearest neighbors, within pca dims
seurat_obj_merged_joined <- FindNeighbors(seurat_obj_merged_joined, reduction = reduct_integr, dims = dims_pca)

# Find clusters # ?FindClusters, if a clustree is wanted, clustreecheck should be set to TRUE
if (clustree_check){
  
  seurat_obj_merged_joined <- FindClusters(seurat_obj_merged_joined,
                                           # use multiple resolutions
                                           resolution = c(0.1, 0.2,0.3,0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.2, 1.4),
                                           # graph.name = "SCT_nn"
                                           )
} else {
  # find clusters with set resolution, not for clustree
  seurat_obj_merged_joined <- FindClusters(seurat_obj_merged_joined,
                                           resolution = res_clust)
}

# Check what cluster the first 5 cells belong to
# head(Idents(seurat_obj_merged_joined))

if (clustree_check) {
  # install package if needed and load package and seuratobject with multiple resolution graphs
  # install.packages("clustree")
  library(clustree)
  # load(paste0(base_dir,"/data/seurat_objects_list_RNA_clusters_louvain_clustree.RData"))
  
  # create clustree visualization plot in a pdf
  pdf(file = paste0(base_dir,"/figures/explore_seurat/clustree_louvain_rna_",integr_method,".pdf"), height = 15, width = 20)
  clustree(seurat_obj_merged_joined, prefix = "RNA_snn_res.")
  dev.off()
  
  # check the different clustering results
  # grep("res", colnames(seurat_obj_merged_joined@meta.data), value = TRUE)
}

#######################################
# Create UMAP Reduction + CellCycleScoring
#######################################

# Run umap reduction method, with PCA dimensions chosen
seurat_obj_merged_joined <- RunUMAP(seurat_obj_merged_joined, reduction = reduct_integr, dims = dims_pca)

# Run cellcyclescoring
seurat_obj_merged_joined <- CellCycleScoring(seurat_obj_merged_joined, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes)

#######################################
# Save object
#######################################

# save data
save(seurat_obj_merged_joined, file = "/scratch/svc_td_compbio/users/MaDeBa/data/seurat_objects_list_RNA_clusters_louvain_clustree.RData")

print("Saved to: /scratch/svc_td_compbio/users/MaDeBa/data/seurat_objects_list_RNA_clusters_louvain_clustree.RData")
#######################################
# UMAP Visualizations
#######################################

# check possible groupings on according to metadata
# colnames(seurat_obj_merged_joined[[]])

# create subdirectory based on integration(not ending on /)
subdir_umap <- switch(integr_method,
  "none" = "/figures/rna_umap_pca",    
  "harmony" = paste0("/figures/rna_umap_harmony_gb_ident_base"),
  "scVI" = "/figures/rna_umap_scVI_gb_ident_dim27_res07"
) 

if (!dir.exists(file.path(base_dir,subdir_umap))) {
  dir.create(file.path(base_dir,subdir_umap))
  print(paste0("Subdirectory created: ",subdir_umap))
} else {print(paste0("Umap subdirectory already exists: ", subdir_umap))}

# list different group by options (and easier name)
umap_groups <- list(
  "cluster" = "seurat_clusters",
  "patient" = "patient",
  "sample"= "orig.ident",
  "nCount" = "nCount_RNA",
  "nFeature" = "nFeature_RNA",
  "percent_mt" = "percent.mt"
)

# create part of plot title with different parameters
title_end <- paste0("(Reduction = ", integr_method,") (PCA dims = ", min(dims_pca),":",max(dims_pca),", res = ", res_clust, ")")

# visualize marker expression on umap, use color viridis and improve for all celltypes
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
    labs(title = paste("RNA -", grp, title_end)) +
    # labs(title = paste0("RNA - ",grp," ","(theta = ", theta_harmony,")",title_end)) +
    theme(
      plot.title = element_text(hjust=0.5)
      # ,legend.position = "bottom"
    )
  
  # print plot and close pdf
  print(umap_p)
  dev.off()
  # to check what has already been done: (prints in console)
  print(paste0("Umap plot created: ", grp))
}

# Cluster composition per sample
# table(seurat_obj_merged_joined$seurat_clusters, seurat_obj_merged_joined$orig.ident)

# CCS 

# UMAP visualization colored by Cellcyclephases
pdf(file = paste0(base_dir,subdir_umap,"/rna_CCS_umap_phases.pdf"),
    width = 12, height = 9)

DimPlot(seurat_obj_merged_joined,
        reduction = "umap",
        group.by = "Phase") +
  # labs(title = paste0("RNA - ",grp," ","(theta = ", theta_harmony,")",title_end)) +
  labs(title = paste("RNA - phases", title_end)) +
  theme(
    plot.title = element_text(hjust=0.5),
    legend.position = "bottom"
  )
dev.off()

# UMAP visualization split and colored by Cellcyclephases
pdf(file = paste0(base_dir,subdir_umap,"/rna_CCS_umap_split_phases.pdf"),
    width = 15, height = 9)
DimPlot(seurat_obj_merged_joined,
        reduction = "umap",
        split.by = "Phase",
        group.by = "Phase") +
  labs(title = paste("RNA - phases", title_end)) +
  
  theme(
    plot.title = element_text(hjust=0.5),
    legend.position = "bottom"
  )
dev.off()

# UMAP visualization split by Cellcyclephases and colored by seurat clusters
pdf(file = paste0(base_dir,subdir_umap,"/rna_CCS_umap_phase-cluster.pdf"),
    width = 15, height = 9)
DimPlot(seurat_obj_merged_joined,
        reduction = "umap",
        split.by = "Phase",
        group.by = "seurat_clusters",
        label = TRUE) +
  labs(title = paste("RNA - phases + clusters",title_end)) +
  
  theme(
    plot.title = element_text(hjust=0.5),
    legend.position = "none"
  )
dev.off()

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
dir_ftplt <- paste0(base_dir,subdir_umap,"/rna_marker_FeatPlot")
if (!dir.exists(dir_ftplt)) {
  dir.create(dir_ftplt)
} else (print("Umap-featureplot subdirectory already exists."))

for (ct in names(markers)) {
  # select celltype markers for each celltype
  ct_markers <- markers[[ct]]

  # create pdf with celltype in name
  pdf(file = paste0(dir_ftplt,"/rna_umap_ftplt_",ct,".pdf"),
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
pdf(file = paste0(base_dir,subdir_umap,"/rna_markers_dotplot_cluster.pdf"), width = 22, height = 10)
DotPlot(seurat_obj_merged_joined,
        features = markers,
        group.by = "seurat_clusters") + 
  scale_color_viridis_b() +
  labs(title = paste0("Marker expression - clusters - ", integr_method),
       x = "Cell type - gene markers",
       y = "Clusters") +
  theme(plot.title = element_text(hjust=0.5))+
  RotatedAxis()
dev.off()

# print confirmation
print("All visualizations created.")
