rm(list = ls())

library(Seurat)
library(preprocessCore)
# devtools::install_github('sunduanchen/Scissor')
# https://api.github.com/repos/sunduanchen/Scissor/tarball/HEAD
# devtools::install_local('~/Downloads/sunduanchen-Scissor-311560a.tar.gz')
library(Scissor)

load('scRNA_for_scAB_Scissor.Rdata')
scRNA
table(scRNA$celltype)
sce=CreateSeuratObject(
  counts = scRNA@assays$RNA@counts,
  meta.data = scRNA@meta.data
)
# Scissor::Seurat_preprocessing 
sc_dataset <- Scissor::Seurat_preprocessing(
  scRNA@assays$RNA@counts , verbose = T)
sc_dataset$celltype =  scRNA$celltype
DimPlot(sc_dataset, reduction = 'umap', 
        label = T, label.size = 3) +
DimPlot(sc_dataset, reduction = 'umap',
        group.by = 'celltype', label = T, label.size = 3)
ggplot2::ggsave('umap_by_celltype_and_Scissor.pdf',width =8)

DimPlot(scRNA, reduction = 'umap',
        group.by = 'celltype', label = T, label.size = 3)+
  DimPlot(sc_dataset, reduction = 'umap',
          group.by = 'orig.ident', label = T, label.size = 3)


# sc_dataset
# scRNA

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

if(T){
  start_time1 <- Sys.time()
 
  # devtools::install_github('sunduanchen/Scissor')
  library(Seurat)
  library(Scissor)
  infos1 <- Scissor(as.matrix(bulk_dataset), 
                    sc_dataset, phenotype, alpha = 0.05, 
                    family = "cox", 
                    Save_file = './Scissor_survival.RData')
  
  end_time1 <- Sys.time()
  execution_time1 <- end_time1 - start_time1
  print(execution_time1)
}

names(infos1) 
length(infos1$Scissor_pos) 
infos1$Scissor_pos[1:4]
length(infos1$Scissor_neg) 
infos1$Scissor_neg

Scissor_select <- rep("Background", ncol(sc_dataset))
names(Scissor_select) <- colnames(sc_dataset)
Scissor_select[infos1$Scissor_pos] <- "Scissor+"
Scissor_select[infos1$Scissor_neg] <- "Scissor-"
sc_dataset <- AddMetaData(sc_dataset,
                          metadata = Scissor_select, 
                          col.name = "scissor")
table(sc_dataset$scissor)

p1 = DimPlot(sc_dataset, reduction = 'umap', group.by = 'scissor', 
        cols = c('grey','royalblue','indianred1'), pt.size = 1.2, order = T)
p1
p2 = DimPlot(sc_dataset, reduction = 'umap', group='celltype',
            label = T, label.size = 3)
p2
p1+p2

gplots::balloonplot(
  table(sc_dataset$celltype,sc_dataset$scissor)
)
save(infos1,sc_dataset,file = 'output_of_Scissor.Rdata')


if(F){
   
  library(scRNAstat)
  markers_figures <- basic_markers(sc_dataset,
                                   org='human',
                                   group='celltype',
                                   dir = './')
  markers_figures$all_markers + p
  sc_dataset
}

phe_scissor = sc_dataset@meta.data
load('scAB_results.Rdata')
sc_dataset
phe_scAB = sc_dataset@meta.data

gplots::balloonplot(table(phe_scAB$scAB_select,phe_scissor$scissor))
gplots::balloonplot(table(phe_scAB$celltype,phe_scAB$scAB_select))
gplots::balloonplot(table(phe_scissor$celltype,phe_scissor$scissor))




