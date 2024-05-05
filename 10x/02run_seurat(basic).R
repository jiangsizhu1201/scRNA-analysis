library(readr)
library(Seurat)
library(ggplot2)
library(Matrix)
library(data.table)
library(dplyr)
library(stringr)
library(tidyverse)
library(harmony)
library(patchwork)
library(tidydr)
library(ggpubr)
library(fgsea)
library(msigdbr)
library(clusterProfiler)

# Reference:
# Azizi, et al. Single-Cell Map of Diverse Immune Phenotypes in the Breast Tumor Microenvironment. (2018) Cell. doi: 10.1016/j.cell.2018.05.060

###### step1: Read Data ######

raw_data <- fread("Azizi2018_BreastCancer/GSE114725_rna_raw.csv.gz")
data.matrix = Matrix(as.matrix(t(raw_data)), sparse = T) # dgCMatrix 
sc_obj <- CreateSeuratObject(counts = data.matrix, 
                             project = "BreastCancer",
                             meta.data = meta_data)

###### step2: Quality Check ######

### Mitochondria %

sc_obj[["percent_mt"]] <- PercentageFeatureSet(sc_obj, 
                                               pattern = "^MT-")

VlnPlot(sc_obj, 
        features = c("nFeature_RNA", "nCount_RNA", "percent_mt"), 
        ncol = 3, alpha = 0.1)

### Low-complexity cell filtration using linear regression

umis <- sc_obj$nCount_RNA
genes <- sc_obj$nFeature_RNA
residuals <- resid(lm(genes ~ umis))
cutoff <- mean(residuals) - 3 * sd(residuals)
cells_to_keep <- names(residuals)[residuals >= cutoff]

sc_obj$keep_qc <- ifelse(rownames(sc_obj[[]]) %in% cells_to_keep, "keep", "remove")

### Final QC
sc_obj <- subset(sc_obj, subset = (nFeature_RNA > 200) & 
                            (percent_mt < 20) &
                            (nCount_RNA > 500) &
                            (keep_qc == "keep"))

###### step3: Preprocessing + Batch Effect Correction ######

sc_obj <- NormalizeData(sc_obj, 
                        normalization.method = "LogNormalize", 
                        scale.factor = 10000, 
                        verbose = FALSE)

sc_obj <- FindVariableFeatures(sc_obj, 
                               selection.method = "vst", 
                               nfeatures = 2000,
                               verbose = FALSE) 

sc_obj <- ScaleData(sc_obj, 
                    do.scale  = TRUE, 
                    do.center = TRUE,
                    verbose = FALSE) 
## dimensino reduction
### PCA
sc_obj <- RunPCA(sc_obj, verbose = FALSE) 
ElbowPlot(sc_obj) # determine inflection point

### UMAP/TSNE
sc_obj <- RunUMAP(sc_obj, dims = 1:20, check_duplicates = FALSE) 

## perform harmony to correct batch effect
sc_obj = RunHarmony(sc_obj, c("batch"), verbose = FALSE)

sc_obj <- RunUMAP(sc_obj, 
                  dims = 1:20, 
                  reduction = "harmony", 
                  reduction.name = "umap_harmony")
#### visulization
DimPlot(sc_obj, reduction = "umap_harmony", group.by = "batch") +
  theme_dr(xlength = 0.3,
           ylength = 0.3,
           arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(panel.grid = element_blank())

###### step4: Clustering ######

### build the graph
sc_obj <- FindNeighbors(sc_obj, 
                                 dims = 1:20, 
                                 reduction = "harmony",
                                 verbose = FALSE)
### cluster on the graph (Louvain)
sc_obj <- FindClusters(sc_obj, 
                       resolution = seq(0.5, 3, by=0.5), 
                       algorithm = 1,
                       verbose = FALSE) 

### silhouett analysis

cell_dists <- dist(sc_obj@reductions$harmony@cell.embeddings,
                   method = "euclidean")

cluster_info <- sc_obj@meta.data[,grepl(paste0(DefaultAssay(sc_obj),
                                "_snn_res"),colnames(sc_obj@meta.data))] %>%
  dplyr::mutate_all(as.character) %>%
  dplyr::mutate_all(as.numeric)

###  test example
si_0.5= silhouette(cluster_info[, "RNA_snn_res.0.5"], cell_dists) # res = 0.5
mean(si_0.5[,'sil_width']) 

###### step5: Find Marker Genes ######

Idents(sc_obj) = "RNA_snn_res.0.5"

top_markers = FindAllMarkers(sc_obj, 
                             logfc.threshold = 0.5, 
                             test.use = "wilcox", 
                             min.pct = 0.25, 
                             only.pos = TRUE, 
                             verbose = FALSE)

top_10 = top_markers %>% 
  group_by(cluster) %>% 
  slice_head(n=10) %>% 
  ungroup()

# visualization: heatmap of the top 10 genes (based on log2FC) from each cluster
my_col <- viridis(15)
p1 <- DoHeatmap(sc_obj, features = top_10$gene, size = 3, 
                angle = 0, group.colors = my_col, hjust = 1) +
  theme(text = element_text(size = 5))+
  scale_colour_npg() +
  scale_fill_gradient2(low = '#0099CC', mid = 'white', high = '#CC0033') + 
  NoLegend()

# visualization: dotplot of marker genes within each cluster (manual annotation/ customize)
list_genes <- list(
  CD4_T_cell = c("CD3D", "IL7R", "TRAC", "TCF7", "CCR7", "NOSIP"),
  CD8_T_cell = c("CCL5", "TNFAIP3", "CD69"),
  NK_cell = c("PRF1", "FGFBP2"),
  B_cell = c("MS4A1", "IGHM", "CD74", "CD79A", "CD37"),
  Erythroid = c("HBB", "IGFBP7"),
  Monocyte = c("FCN1", "S100A9", "LYZ", "SERPINA1", "CTSS"),
  Macrophage_cell = c("C1QA", "CD68", "CTSB", "GRN", "NPC2"),
  Mast_cell = c("TPSAB1", "CPA3", "HDC", "MS4A2", "GATA2"),
  CD14_Monocyte = c("FCGR3B", "CSF3R", "S100A8"),
  Plasma_cell = c("MZB1", "IGHA2", "TNFRSF17", "DERL3"),
  Dendritic_cell = c("HLA-DPB1", "HLA-DPA1", "HLA-DQA1"),
  Cycling = c("TOP2A", "MKI67", "STMN1", "TPX2"),
  Plasmacytoid_dendritric_cell = c("JCHAIN", "IRF7", "TCF4", "LILRA4", "GZMB"),
  non_immune = c("MGP", "DCN", "APOD", "COL3A1", "COL6A", "CXCL12")
)

p1 <- DotPlot(sc_obj, features = list_genes) + RotatedAxis() +
  theme(panel.border = element_rect(color = "black"),
        panel.spacing = unit(1, "mm"),
        strip.text.x = element_text(size = 7),
        strip.text = element_text(margin = margin(b = 3, unit = "mm")),
        strip.placement = "outlet",
        axis.text.x = element_text(size = 4), 
        axis.text.y = element_text(size = 8),
        legend.title = element_text(size = 10), 
        legend.text = element_text(size = 8), 
        axis.line = element_blank()) +
  labs(x="", y="")

p2 <- ggplotGrob(p1)
lg <- linesGrob(x=unit(c(0,1),"npc"), y=unit(c(0,0)+0.2,"npc"),
                gp=gpar(col="black", lwd=4))

grid.newpage(); for (k in grep("strip-t",p2$layout$name)) {
  p2$grobs[[k]]$grobs[[1]]$children[[1]] <- lg}

grid.draw(p2)

###### step6: Differential Abundance ######

cell_types <- c("CD4_T_cell", "CD8_T_cell", "NK_cell", "B_cell", 
                "Macrophage", "Monocyte", "Erythroid", "Dendritic_cell", "Mast_cells",
                "CD14+_Monocyte", "Plasmacytoid_dendritri_cell", "Cycling", "Non-immune")

## calculate proportion across all cell types
metadata_filtered <- sc_obj[[]] %>%
  filter(tissue %in% c("NORMAL", "TUMOR"))
proportion_all <- data.frame()

for (cell in cell_types2) {
  temp <- metadata_filtered %>%
    group_by(patient, tissue) %>%
    summarise(
      total = n(),
      count = sum(cell_type == cell),
      proportion = count / total,
      .groups = 'drop'
    ) %>%
    mutate(cell_type = cell) 
  proportion_all <- bind_rows(proportion_all, temp)
}

ggplot(proportion_all, aes(x = tissue, y = proportion, fill = tissue)) +
  geom_boxplot() +
  facet_wrap(~ cell_type, scales = "free_y") +
  stat_compare_means(aes(group=tissue), label="p.signif") +
  labs(title = "",
       x = "Tissue",
       y = "Cell Type Proportion") +
  theme_classic() +
  scale_fill_brewer(palette = "Set2") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title.x = element_text(hjust = 0.5, face = "bold"), 
    axis.title.y = element_text(hjust = 0.5, face = "bold"), 
  )

###### step7: Differential Expression Analysis ######

### pseudobulk
pb_sc_obj <- AggregateExpression(sc_obj, assays = "RNA", 
                                 return.seurat = T, group.by = c("patient","cell_type", "tissue"))

pb_sc_obj$celltype_tissue <- paste(pb_sc_obj$cell_type, pb_sc_obj$tissue, sep = "_")
Idents(pb_sc_obj) <- "celltype_tissue"

###  example: DE testing on the pseudobulk level for CD4+ T cell
bulk.cd4_t <- FindMarkers(object = pb_sc_obj, 
                          ident.1 = "CD4-T-cell_TUMOR", 
                          ident.2 = "CD4-T-cell_NORMAL",
                          logfc.threshold = 1,
                          test.use = "DESeq2")

bulk.cd4_t %>% arrange(p_val, -avg_log2FC) %>% head(10)

genes.to.label <- bulk.cd4_t %>% 
  arrange(p_val, desc(avg_log2FC)) %>%  
  slice_head(n = 10) %>%               
  row.names()                         

VlnPlot(pb_sc_obj, features = genes.to.label[1:4], 
        idents = c("CD4-T-cell_NORMAL", "CD4-T-cell_TUMOR"), 
        ncol = 2,
        group.by = "tissue") 

### functinal enrichment analysis using fgsea

genesets <- msigdbr(species = "Homo sapiens", 
                    category = "C5", 
                    subcategory = "BP") 

genesets <- genesets %>% 
  split(x = .$gene_symbol, f = .$gs_name)

de.genes <- bulk.cd4_t %>% 
  rownames_to_column(var = "gene") %>% 
  arrange(-avg_log2FC) %>%
  select(gene, avg_log2FC) %>% 
  deframe()

gsea_out <- fgsea(genesets, stats = de.genes)
gsea_df <- gsea_out %>% arrange(pval) %>% filter(NES>0)

gsea_df$pathway_w <- tolower(gsub("_", " ", sub("GOBP_", "", gsea_df$pathway)))

ggdotchart(gsea_df[0:10,], x = "pathway_w", y = "NES", 
           color = "pval",
           sorting = "descending",
           add = "segments",
           rotate = TRUE,
           dot.size = 8,
           gradient.cols = c("darkgreen", "grey"),
           font.label = list(color = "white", size = 4, vjust = 0.5), 
           ggtheme = theme_pubr(),
           ylab = F)+
  geom_hline(yintercept = 0, linetype=2, linewidth=0.5)+
  theme(legend.title = element_text(size = 10), 
        legend.text = element_text(size = 5)) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 30))
