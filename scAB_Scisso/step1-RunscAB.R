rm(list = ls())

library(Seurat)
library(preprocessCore)
# devtools::install_github("jinworks/scAB")
# devtools::install_local('~/Downloads/jinworks-scAB-d142a20.tar.gz')
library(scAB)
# example data: phenotype, bulk_dataset

head(phenotype)
bulk_dataset[1:4,1:4]

load('scRNA_for_scAB_Scissor.Rdata')
# scRNA
# table(scRNA$celltype)
p1=DimPlot(scRNA, reduction ="umap", group.by="celltype",label = T)
p1
# lung cancer: post-harmony

load("../01-tcga_luad_from_xena/tcga-luad.for_survival.rdata")
head(meta)
exprSet[1:4,1:4]

bulk_dataset = exprSet
head(bulk_dataset[,1:10])
dim(bulk_dataset)
phenotype = meta[,c(3,2)]
colnames(phenotype)=c("time","status")
head(phenotype)
table(phenotype$status)
identical(colnames(bulk_dataset) ,row.names(phenotype))


# dimension reduction w/o harmony
sc_dataset <- scAB::run_seurat(scRNA,verbose = FALSE) 
p2=DimPlot(sc_dataset, reduction ="umap", group.by="orig.ident",label = T)
UMAP_celltype <- DimPlot(sc_dataset, reduction ="umap",
                         group.by="celltype",label = T)
# UMAP_celltype
# p1+p2
# ggplot2::ggsave('UMAP_harmony_or_not.pdf',width = 10)

# w/o harmony: heterogeneity across different patients

pdf('orig.ident-vs-phenotype.pdf')
gplots::balloonplot(
  table(sc_dataset$orig.ident,sc_dataset$celltype)
)
dev.off()

### cluster needed; computation-intense
if(T){
  start_time1 <- Sys.time()
  scAB_data <- create_scAB(sc_dataset,bulk_dataset,phenotype)
  
  K <- select_K(scAB_data)
  K
  scAB_result <- scAB(Object=scAB_data, K=K)
  sc_dataset <- findSubset(sc_dataset, 
                           scAB_Object = scAB_result, 
                           tred = 2)
  table(sc_dataset$scAB_select)
  end_time1 <- Sys.time()
  execution_time1 <- end_time1 - start_time1
  print(execution_time1)
}

save(sc_dataset,file = 'scAB_results.Rdata')

load(file = 'scAB_results.Rdata')
table(Idents(sc_dataset))
table(sc_dataset$scAB_select)
markers <- FindMarkers(sc_dataset, ident.1 = "scAB+ cells", 
                       group.by = 'scAB_select', 
                       logfc.threshold = 0.2)
markers <- markers[which(markers$p_val_adj<0.05),]
head(markers) 
UMAP_scAB <- DimPlot(sc_dataset,group.by="scAB_select",
                     cols=c("#80b1d3","red"),
                     pt.size=0.001,
                     order=c("scAB+ cells","Other cells"))
# UMAP_scAB
patchwork::wrap_plots(plots = list(UMAP_celltype,UMAP_scAB), ncol = 2)
ggplot2::ggsave('UMAP_scAB.pdf',width = 10)
gplots::balloonplot(
  table(sc_dataset$scAB_select,sc_dataset$celltype)
)





