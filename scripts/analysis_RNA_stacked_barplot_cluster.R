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

# load(paste0(base_dir,"/data/seurat_objects_list_RNA_umap.RData"))

# only if not done in advance
# create patient column with patient id, extracted from orig.ident
# seurat_obj_merged_joined@meta.data$patient <- sub(".*-(Patient[0-9]+).*", "\\1", seurat_obj_merged_joined@meta.data$orig.ident)

# copy metadata to create new dataframe
metadata <- seurat_obj_merged_joined@meta.data

# new column, each cell counts as 1
metadata$nCells_RNA <- 1

# create pdf
pdf(file = paste0(base_dir,"/figures/explore_seurat/rna_stackbp_patient.pdf"), width = 10, height = 12)

# make stacked barplot for the different clusters, with bars filled with the amount of cells and colored by patient
ggplot(metadata,
       aes(fill=patient,y=nCells_RNA,x=seurat_clusters)) +
  # add stacked bars
  geom_bar(position="fill", stat="identity") + 
  # add theme
  theme_classic() + 
  # add custom labels for x, y and legend
  labs(title = "Stacked barplot - patient",
       x = "Seurat_clusters",
       y = "Proportion of cells",
       fill = "Patient") +
  # change legend and title position
  theme(plot.title = element_text(hjust=0.5),
        legend.position = "bottom")
# close pdf
dev.off()

# same as above, but for sample
pdf(file = paste0(base_dir,"/figures/explore_seurat/rna_stackbp_sample.pdf"), width = 10, height = 12)
ggplot(metadata,
       aes(fill=orig.ident,y=nCells_RNA,x=seurat_clusters)) +
  geom_bar(position="fill", stat="identity") + theme_classic() +  
  labs(title = "Stacked barplot - sample",
       x = "Seurat_clusters",
       y = "Proportion of cells",
       fill = "Sample") +
  theme(plot.title = element_text(hjust=0.5),
        legend.position = "bottom")+
  guides(fill = guide_legend(nrow = 12, byrow = FALSE))
dev.off()
