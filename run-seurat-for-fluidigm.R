###### step1: 获取单细胞转录组数据集 ######

suppressMessages(library(scater))
suppressMessages(library(scRNAseq))
library(ggplot2)
library(tidyr)
library(cowplot)
library("FactoMineR")
library("factoextra")
library("ROCR")

## ----- Load Example Data -----
fluidigm <- ReprocessedFluidigmData()
fluidigm

ct <- floor(assays(fluidigm)$rsem_counts)
ct[1:4,1:4] 
sample_ann <- as.data.frame(colData(fluidigm))
DT::datatable(sample_ann)

box <- lapply(colnames(sample_ann[,1:19]),function(i) {
  dat <-  sample_ann[,i,drop=F] 
  dat$sample=rownames(dat)
  ## 画boxplot 
  ggplot(dat, aes('all cells', get(i))) +
    geom_boxplot() +
    xlab(NULL)+ylab(i)
})
plot_grid(plotlist=box, ncol=5 )

library(Seurat)
sce.all=CreateSeuratObject(counts = ct)
sce.all 
as.data.frame(sce.all@assays$RNA@counts[1:10, 1:2])
head(sce.all@meta.data, 10)
table(sce.all@meta.data$orig.ident) 



###### step2:QC质控 ######

dir.create("./1-QC")
setwd("./1-QC") 
#计算线粒体基因比例
# 人和鼠的基因名字稍微不一样 
mito_genes=rownames(sce.all)[grep("^MT-", rownames(sce.all))] 
mito_genes #13个线粒体基因
sce.all=PercentageFeatureSet(sce.all, "^MT-", col.name = "percent_mito")
fivenum(sce.all@meta.data$percent_mito)
#计算核糖体基因比例
ribo_genes=rownames(sce.all)[grep("^Rp[sl]", rownames(sce.all),ignore.case = T)]
ribo_genes
sce.all=PercentageFeatureSet(sce.all, "^RP[SL]", col.name = "percent_ribo")
fivenum(sce.all@meta.data$percent_ribo)
#计算红血细胞基因比例
rownames(sce.all)[grep("^Hb[^(p)]", rownames(sce.all),ignore.case = T)]
sce.all=PercentageFeatureSet(sce.all, "^HB[^(P)]", col.name = "percent_hb")
fivenum(sce.all@meta.data$percent_hb)
#可视化细胞的上述比例情况
feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "percent_hb")
feats <- c("nFeature_RNA", "nCount_RNA")
p1=VlnPlot(sce.all, group.by = "orig.ident", features = feats, pt.size = 0.01, ncol = 2) + 
  NoLegend()
p1
library(ggplot2) 
ggsave(filename="Vlnplot1.png",plot=p1)
feats <- c("percent_mito", "percent_ribo", "percent_hb")
p2=VlnPlot(sce.all, group.by = "orig.ident", features = feats, pt.size = 0.01, ncol = 3, same.y.lims=T) + 
  scale_y_continuous(breaks=seq(0, 100, 5)) +
  NoLegend()
p2	
ggsave(filename="Vlnplot2.png",plot=p2)

p3=FeatureScatter(sce.all, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident", pt.size = 0.5)
ggsave(filename="Scatterplot.png",plot=p3)
#根据上述指标，过滤低质量细胞/基因
#过滤指标1:最少表达基因数的细胞&最少表达细胞数的基因
selected_c <- WhichCells(sce.all, expression = nFeature_RNA > 300)
selected_f <- rownames(sce.all)[Matrix::rowSums(sce.all@assays$RNA@counts > 0 ) > 3]

sce.all.filt <- subset(sce.all, features = selected_f, cells = selected_c)
dim(sce.all) 
dim(sce.all.filt) 
#  可以看到，主要是过滤了基因，其次才是细胞

# par(mar = c(4, 8, 2, 1))
C=sce.all.filt@assays$RNA@counts
dim(C)
C=Matrix::t(Matrix::t(C)/Matrix::colSums(C)) * 100
# 这里的C 这个矩阵，有一点大，可以考虑随抽样
C=C[,sample(1:ncol(C),100)]
most_expressed <- order(apply(C, 1, median), decreasing = T)[50:1]
pdf("TOP50_most_expressed_gene.pdf",width=14)
boxplot(as.matrix(Matrix::t(C[most_expressed, ])),
        cex = 0.1, las = 1, 
        xlab = "% total count per cell", 
        col = (scales::hue_pal())(50)[50:1], 
        horizontal = TRUE)
dev.off()
rm(C)

#细胞周期评分
sce.all.filt = NormalizeData(sce.all.filt)
s.genes=Seurat::cc.genes.updated.2019$s.genes
g2m.genes=Seurat::cc.genes.updated.2019$g2m.genes
sce.all.filt=CellCycleScoring(object = sce.all.filt, 
                              s.features = s.genes, 
                              g2m.features = g2m.genes, 
                              set.ident = TRUE)
p4=VlnPlot(sce.all.filt, features = c("S.Score", "G2M.Score"), group.by = "orig.ident", 
           ncol = 2, pt.size = 0.1)
p4

ggsave(filename="Vlnplot4_cycle.png",plot=p4)
sce.all.filt@meta.data  %>% ggplot(aes(S.Score,G2M.Score))+geom_point(aes(color=Phase))+
  theme_minimal()
ggsave(filename="cycle_details.pdf" )
# S.Score较高的为S期，G2M.Score较高的为G2M期，都比较低的为G1期

###### step3: seurat 标准流程（降维聚类分群） ######

sce.all.filt = FindVariableFeatures(sce.all.filt)
sce.all.filt = ScaleData(sce.all.filt, 
                         vars.to.regress = c("nFeature_RNA", "percent_mito"))
sce.all.filt = RunPCA(sce.all.filt, npcs = 20)
sce.all.filt = RunTSNE(sce.all.filt, npcs = 20)
sce.all.filt = RunUMAP(sce.all.filt, dims = 1:10)

sce.all.filt
sce.all.filt <- FindNeighbors(sce.all.filt, dims = 1:15)
sce.all.filt <- FindClusters(sce.all.filt, resolution = 0.8)
table(sce.all.filt@meta.data$RNA_snn_res.0.8)  
 
DimPlot(sce.all.filt, 
        reduction = 'umap')
ggsave('first_umap_by_seurat_cluster.pdf')

library(ggplot2) 
genes_to_check = c('PTPRC', 'CD3D', 'CD3E', 'CD4','CD8A','CD19', 'CD79A', 'MS4A1' ,
                   'IGHG1', 'MZB1', 'SDC1',
                   'CD68', 'CD163', 'CD14', 
                   'TPSAB1' , 'TPSB2',  # mast cells,
                   'RCVRN','FPR1' , 'ITGAM' ,
                   'FGF7','MME', 'ACTA2',
                   'PECAM1', 'VWF', 
                   'EPCAM' , 'KRT19', 'PROM1', 'ALDH1A1' )
library(stringr)  
p_all_markers <- DotPlot(sce.all.filt, features = genes_to_check,
                         assay='RNA'  )  + coord_flip()

p_all_markers
ggsave(plot=p_all_markers,
       filename="first_check_all_marker_by_seurat_cluster.pdf",width = 12)


###### step4: 注释（手动） ######


sce.all=sce.all.filt
colnames(sce.all@meta.data)
# 需要自行看图，定细胞亚群：  
celltype=data.frame(ClusterID=0:2,
                    celltype=0:2) 
celltype[celltype$ClusterID %in% c(2),2]='myeloid'
celltype[celltype$ClusterID %in% c(0 ),2]='epithelial' 
celltype[celltype$ClusterID %in% c( 1),2]='fibo'   

head(celltype)
celltype 
table(celltype$celltype)
sce.all@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  sce.all@meta.data[which(sce.all@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(sce.all@meta.data$celltype)

genes_to_check = c('PTPRC', 'CD3D', 'CD3E', 'CD4','CD8A','CD19', 'CD79A', 'MS4A1' ,
                   'IGHG1', 'MZB1', 'SDC1',
                   'CD68', 'CD163', 'CD14', 
                   'TPSAB1' , 'TPSB2',  # mast cells,
                   'RCVRN','FPR1' , 'ITGAM' ,
                   'FGF7','MME', 'ACTA2',
                   'PECAM1', 'VWF', 
                   'EPCAM' , 'KRT19', 'PROM1', 'ALDH1A1' )
p <- DotPlot(sce.all, features = genes_to_check,
             assay='RNA' ,group.by = 'celltype' )  + coord_flip()

p
ggsave(plot=p, filename="check_marker_by_celltype.pdf")
table(sce.all@meta.data$celltype,sce.all@meta.data$seurat_clusters)

DimPlot(sce.all, reduction = "umap", group.by = "celltype",label = T)  
ggsave('umap_by_celltype.pdf')

library(patchwork)
p_all_markers=DotPlot(sce.all, features = genes_to_check,
                      assay='RNA' ,group.by = 'celltype' )  + coord_flip()
p_umap=DimPlot(sce.all, reduction = "umap", group.by = "celltype",label = T) 
p_all_markers+p_umap
ggsave('markers_umap_by_celltype.pdf',units = 'cm',
       width = 22,height = 11)


###### step5: 注释（自动） ######


###### step6: 差异基因 ######

sce=sce.all
sce
Idents(sce)=sce@meta.data$seurat_clusters
table(Idents(sce))  
sce.markers <- FindAllMarkers(object = sce, only.pos = TRUE, min.pct = 0.25, 
                              thresh.use = 0.25)

DT::datatable(sce.markers)
pro='fluidigm'
write.csv(sce.markers,file=paste0(pro,'_sce.markers.csv'))
library(dplyr) 
top10 <- sce.markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
DoHeatmap(sce,top10$gene,size=3)
ggsave(filename=paste0(pro,'_sce.markers_heatmap.pdf'))

# 数据集验证，取决于单细胞是否提供
table(sce$celltype,sample_ann$Biological_Condition)












