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
# 下面的  Scissor::Seurat_preprocessing 函数需要的仅仅是表达量矩阵
sc_dataset <- Scissor::Seurat_preprocessing(
  scRNA@assays$RNA@counts , verbose = T)
sc_dataset$celltype =  scRNA$celltype
DimPlot(sc_dataset, reduction = 'umap', 
        label = T, label.size = 3) +
DimPlot(sc_dataset, reduction = 'umap',
        group.by = 'celltype', label = T, label.size = 3)
ggplot2::ggsave('umap_by_celltype_and_Scissor.pdf',width =8)
# 可以看到9个肺腺癌病人的上皮细胞被harmony整合后
# 上皮细胞可以区分出来正常上皮细胞亚群
# 以及一些未知的恶性肿瘤细胞亚群
DimPlot(scRNA, reduction = 'umap',
        group.by = 'celltype', label = T, label.size = 3)+
  DimPlot(sc_dataset, reduction = 'umap',
          group.by = 'orig.ident', label = T, label.size = 3)
# 可以看到不同肿瘤病人上皮细胞harmony处理与否
# 会导致病人个体差异被抹平
# 如果不走harmony整合，我们就需要针对每个病人研究肿瘤内部异质性
# 如果走harmony流程，就可以研究整体异质性

sc_dataset
scRNA

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

## 下面是主体代码，会很耗费时间和计算机资源
# 我的电脑超过了半个小时
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
# 是可以看到，不同算法，结果很难有很好的一致性
gplots::balloonplot(table(phe_scAB$scAB_select,phe_scissor$scissor))
gplots::balloonplot(table(phe_scAB$celltype,phe_scAB$scAB_select))
gplots::balloonplot(table(phe_scissor$celltype,phe_scissor$scissor))




