# Stacked barplot generator
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
# print("All packages installed.")

library(Seurat)
library(ggplot2)
set.seed(1234)

# required: clustering already done

#######################################
# Stacked barplot
#######################################
# store folder location
base_dir <- "/scratch/svc_td_compbio/users/MaDeBa"
# subdir umap
# subdir_umap <- "/figures/rna_umap_harmony_gb_pat"

# load(paste0(base_dir,"/data/seurat_objects_list_RNA_umap.RData"))

# only if not done in advance (create patient column with patient id, extracted from orig.ident)
# seurat_obj_merged_joined@meta.data$patient <- sub(".*-(Patient[0-9]+).*", "\\1", seurat_obj_merged_joined@meta.data$orig.ident)

# copy metadata to create new dataframe
metadata <- seurat_obj_merged_joined@meta.data

# new column, each cell counts as 1
metadata$nCells_RNA <- 1

fill_var_all <- c("patient", "sample")
for (fill_var in fill_var_all) {
  
  # store the pdf in the same location as umap visualizations and auto include var name
  pdf_fileloc <- paste0(base_dir,subdir_umap,"/rna_stackbp_",fill_var,"_harmony.pdf")
  pdf(file = pdf_fileloc,
      width = 10, height = 12)
  
  # make stacked barplot for the different clusters, with bars filled with the amount of cells and colored by patient
  if (fill_var == "sample") {
    p1 <- ggplot(metadata,
                 aes(fill=orig.ident,y=nCells_RNA,x=seurat_clusters))
  }
  if (fill_var == "patient") {
    p1 <- ggplot(metadata,
                 aes(fill=patient,y=nCells_RNA,x=seurat_clusters))
  }
  p1 <- p1 +
      # add stacked bars
      geom_bar(position="fill", stat="identity") +  
      # add theme
      theme_classic() + 
      # add custom labels for x, y and legend
      labs(title = paste0("Stacked barplot - ",fill_var," - harmony - th4"),
         x = "Seurat_clusters",
         y = "Proportion of cells",
         fill = fill_var) +
      # change legend and title position
      theme(plot.title = element_text(hjust=0.5),
            legend.position = "bottom")
  # add guide for large amount of samples (displayed in legend)
  if (fill_var == "sample")
    p1 <- p1 + guides(fill = guide_legend(nrow = 12, byrow = FALSE))
  # print plot
  print(p1)
  # close pdf
  dev.off()
}