# scRNA-seq - clustering and UMAP
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
# library(readxl)
library(openxlsx)
library(viridis)
set.seed(1234)

#######################################
# Load in the data
#######################################

# store location where the data, figures and scripts are located, without / on the end
base_dir <- "/scratch/svc_td_compbio/users/MaDeBa"
# dir.exists(base_dir) #check if file location exists

# load in PCA data or integrated data
# PCA
# load(paste0(base_dir,"/data/seurat_objects_list_RNAonly_filtered_merged_normalized_HVG_joined_pca.RData"))
# integrated
# load(paste0(base_dir,"/data/seurat_objects_list_RNA_pca_integrated.RData"))
load(paste0(base_dir,"/data/seurat_objects_list_RNA_pca_integrated_harmony_nojoin.RData"))

#######################################
# Neighbors/clustering
#######################################

# pdf(file = "/scratch/svc_td_compbio/users/MaDeBa/figures/explore_seurat/rna_PCA_elbowplot_40dim_integr.pdf", width = 10, height = 6)
# ElbowPlot(seurat_obj_merged, ndims = 40)+
#   geom_line(aes(x=dims,y=y_data)) +
#   labs(title = "Elbow plot of PCA variance") +
#   theme(plot.title = element_text(hjust=0.5))
# dev.off()

# Find nearest neighbors, within pca dims 1:30 
# seurat_obj_merged_joined <- FindNeighbors(seurat_obj_merged_joined, dims = 1:30)
seurat_obj_merged_joined <- FindNeighbors(seurat_obj_merged_joined, reduction = "harmony", dims = 1:27)

# Find clusters
# seurat_obj_merged_joined <- FindClusters(seurat_obj_merged_joined, resolution = 0.5)
seurat_obj_merged_joined <- FindClusters(seurat_obj_merged_joined, resolution = 0.6)


# Check what cluster the first 5 cells belong to
# head(Idents(seurat_obj_merged_joined), 5)

#######################################
# UMAP Visualization
#######################################

# Run umap reduction method, with PCA dimensions chosen
seurat_obj_merged_joined <- RunUMAP(seurat_obj_merged_joined, reduction = "harmony", dims = 1:27)

# check possible groupings on according to metadata
# colnames(seurat_obj_merged_joined[[]])
# seurat_obj_merged_joined@meta.data$patient <- sub(".*-(Patient[0-9]+).*", "\\1", seurat_obj_merged_joined@meta.data$orig.ident)

# create subdirectory (not ending on /), chooses what to group by (also renamed) , visualize marker expression on umap, use color viridis and improve for all celltypes
subdir_umap <- "/figures/rna_umap_harmony_gb_ident"

if (!dir.exists(file.path(base_dir,subdir_umap))) {
  dir.create(file.path(base_dir,subdir_umap))
} else (print("Umap subdirectory already exists."))

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

pdf(file = paste0(base_dir,subdir_umap,"/rna_CCS_umap_phases.pdf"),
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

pdf(file = paste0(base_dir,subdir_umap,"/rna_CCS_umap_split_phases.pdf"),
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

pdf(file = paste0(base_dir,subdir_umap,"/rna_CCS_umap_phase-cluster.pdf"),
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
pdf(file = paste0(base_dir,subdir_umap,"/rna_markers_dotplot_cluster.pdf"), width = 22, height = 10)
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

