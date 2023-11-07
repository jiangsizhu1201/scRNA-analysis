rm(list = ls())

library(Seurat)
library(preprocessCore)
# devtools::install_github("jinworks/scAB")
# https://api.github.com/repos/jinworks/scAB/tarball/HEAD
# devtools::install_local('~/Downloads/jinworks-scAB-d142a20.tar.gz')
library(scAB)
# 这个scAB包自带了数据集 phenotype和bulk_dataset
# 加载scAB包即可查看它，有助于我们理解其输入文件的要求
head(phenotype)
bulk_dataset[1:4,1:4]

load('scRNA_for_scAB_Scissor.Rdata')
scRNA
table(scRNA$celltype)
p1=DimPlot(scRNA, reduction ="umap", group.by="celltype",label = T)
p1
# 可以看到9个肺腺癌病人的上皮细胞被harmony整合后
# 上皮细胞可以区分出来正常上皮细胞亚群
# 以及一些未知的恶性肿瘤细胞亚群


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


# 单细胞降维聚类分群(默认流程，无harmony)
sc_dataset <- scAB::run_seurat(scRNA,verbose = FALSE) 
p2=DimPlot(sc_dataset, reduction ="umap", group.by="orig.ident",label = T)
p2
UMAP_celltype <- DimPlot(sc_dataset, reduction ="umap",
                         group.by="celltype",label = T)
UMAP_celltype
p1+p2
ggplot2::ggsave('UMAP_harmony_or_not.pdf',width = 10)
# 可以看到不同肿瘤病人上皮细胞harmony处理与否
# 会导致病人个体差异被抹平
# 如果不走harmony整合，我们就需要针对每个病人研究肿瘤内部异质性
# 如果走harmony流程，就可以研究整体异质性


pdf('orig.ident-vs-phenotype.pdf')
gplots::balloonplot(
  table(sc_dataset$orig.ident,sc_dataset$celltype)
)
dev.off()

## 下面是主体代码，会很耗费时间和计算机资源
# 我自己的电脑耗时31分钟
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
UMAP_scAB
patchwork::wrap_plots(plots = list(UMAP_celltype,UMAP_scAB), ncol = 2)
ggplot2::ggsave('UMAP_scAB.pdf',width = 10)
gplots::balloonplot(
  table(sc_dataset$scAB_select,sc_dataset$celltype)
)
# 还是可以看到cycle以及1和7这3个亚群是有特殊性




