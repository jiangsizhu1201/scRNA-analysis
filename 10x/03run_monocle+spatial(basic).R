library(Seurat)
library(ggplot2)
library(monocle3)
library(SeuratWrappers)

# Trajectory analysis

# Reference
# Pellin, et al. 2019](https://www.nature.com/articles/s41467-019-10291-0). This is a hematopoiesis dataset, i.e., blood development.

##### step1: Preprocessing #####

### After QC and normalization

sc_obj <- FindVariableFeatures(sc_obj, selection.method = "vst", nfeatures = 500)
sc_obj <- ScaleData(sc_obj, do.scale  = TRUE, do.center = TRUE)
sc_obj <- RunPCA(sc_obj)
ElbowPlot(sc_obj, ndims = 20, reduction = "pca")

### Make a monocle3 object
sc_cds <- as.cell_data_set(sc_obj)

### Pre-process the data for monocle3
sc_cds <- preprocess_cds(sc_cds, method = "PCA", num_dim = 15, 
                         norm_method = "log", use_genes = VariableFeatures(sc_obj))

##### step2: Batch correction #####
sc_cds <- align_cds(sc_cds, alignment_group = "batch")

### UMAP
sc_cds <- reduce_dimension(sc_cds, reduction_method = "UMAP")
plot_cells(sc_cds, 
           label_groups_by_cluster=FALSE,  
           color_cells_by = "cell_type", 
           reduction_method = "UMAP", 
           group_label_size = 6)

##### step3: Calculate the trajectory #####
sc_cds <- cluster_cells(sc_cds)
sc_cds <- learn_graph(sc_cds, use_partition = F)

### Compare visualization
plot_cells(sc_cds, 
           label_groups_by_cluster=FALSE,  
           color_cells_by = "cell_type", 
           reduction_method = "UMAP", 
           group_label_size = 6)


### Calculate pseudotime

sc_cds <- order_cells(sc_cds)
plot_cells(sc_cds, 
           label_groups_by_cluster=FALSE,  
           color_cells_by = "pseudotime", 
           reduction_method = "UMAP", 
           group_label_size = 6)

### Compare gene expression and pseudotime

rowData(sc_cds)$gene_name <- rownames(sc_cds)
rowData(sc_cds)$gene_short_name <- rowData(sc_cds)$gene_name
plot_cells(sc_cds, 
           genes = c("...."), # gene names
           label_groups_by_cluster=FALSE,
           reduction_method = "UMAP",
           show_trajectory_graph = FALSE)

# Spatial analysis
library(Rfast2)
library(glmGamPoi)

# Reference
# 10x's Visium platform from prostate cancer

##### step1: Preprocessing #####

### normalization, find variable features, and scale
obj <- SCTransform(obj, assay = "Spatial", verbose = F)

##### step2: Visualizing genes in spatial context #####
### example: KLK3(found in epithelial cells) and ACTA2(found in smooth muscle cells)
SpatialFeaturePlot(obj, features = c("KLK3", "ACTA2"))

##### step3: Clustering the data #####
obj <- RunPCA(obj, assay = "SCT", verbose = F)
obj <- RunUMAP(obj, reduction = "pca", dims = 1:30)
obj <- FindNeighbors(obj, reduction = "pca", dims = 1:30)
obj <- FindClusters(obj, verbose = FALSE, resolution = 0.5)

SpatialDimPlot(obj, label = T, label.size = 4)

##### step4: Find spatially variable genes #####

obj <- FindSpatiallyVariableFeatures(obj, assay = "SCT", 
                                     features = VariableFeatures(obj)[1:1000], 
                                     selection.method = "moransi")

###  top spatially differentially expressed genes

table <- obj@assays$SCT@meta.features
table %>% filter(moransi.spatially.variable == TRUE) %>% 
  arrange(moransi.spatially.variable.rank) %>% head() # tail()

