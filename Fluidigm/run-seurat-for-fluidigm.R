suppressMessages(library(scater))
suppressMessages(library(scRNAseq))
library(ggplot2)
library(tidyr)
library(cowplot)
library("FactoMineR")
library("factoextra")
library("ROCR")

## ----- Load Example Data -----

###### step1 ######

fluidigm <- ReprocessedFluidigmData()
ct <- floor(assays(fluidigm)$rsem_counts)
sample_ann <- as.data.frame(colData(fluidigm))
DT::datatable(sample_ann)

box <- lapply(colnames(sample_ann[,1:19]),function(i) {
  dat <-  sample_ann[,i,drop=F] 
  dat$sample=rownames(dat)
  ## boxplot 
  ggplot(dat, aes('all cells', get(i))) +
    geom_boxplot() +
    xlab(NULL)+ylab(i)
})
plot_grid(plotlist=box, ncol=5 )

library(Seurat)
sce.all=CreateSeuratObject(counts = ct)

###### step2:QC ######

## calculate mitochondrial 
mito_genes=rownames(sce.all)[grep("^MT-", rownames(sce.all))] 

sce.all=PercentageFeatureSet(sce.all, "^MT-", col.name = "percent_mito")
fivenum(sce.all@meta.data$percent_mito)

## calculate ribosome 
ribo_genes=rownames(sce.all)[grep("^Rp[sl]", rownames(sce.all),ignore.case = T)]
sce.all=PercentageFeatureSet(sce.all, "^RP[SL]", col.name = "percent_ribo")
fivenum(sce.all@meta.data$percent_ribo)

# calculate hemoglobin
rownames(sce.all)[grep("^Hb[^(p)]", rownames(sce.all),ignore.case = T)]
sce.all=PercentageFeatureSet(sce.all, "^HB[^(P)]", col.name = "percent_hb")
fivenum(sce.all@meta.data$percent_hb)

# visualization 
feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "percent_hb")
feats <- c("nFeature_RNA", "nCount_RNA")
p1=VlnPlot(sce.all, group.by = "orig.ident", features = feats, pt.size = 0.01, ncol = 2) + 
  NoLegend()

feats <- c("percent_mito", "percent_ribo", "percent_hb")

p2=VlnPlot(sce.all, group.by = "orig.ident", features = feats, pt.size = 0.01, ncol = 3, same.y.lims=T) + 
  scale_y_continuous(breaks=seq(0, 100, 5)) +
  NoLegend()

p3=FeatureScatter(sce.all, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident", pt.size = 0.5)

# filter low-quality cells

selected_c <- WhichCells(sce.all, expression = nFeature_RNA > 300)
selected_f <- rownames(sce.all)[Matrix::rowSums(sce.all@assays$RNA@counts > 0 ) > 3]

sce.all.filt <- subset(sce.all, features = selected_f, cells = selected_c)

# par(mar = c(4, 8, 2, 1))
C=sce.all.filt@assays$RNA@counts
dim(C)
C=Matrix::t(Matrix::t(C)/Matrix::colSums(C)) * 100
#### random sampling (C too large) 
C=C[,sample(1:ncol(C),100)]
most_expressed <- order(apply(C, 1, median), decreasing = T)[50:1]
boxplot(as.matrix(Matrix::t(C[most_expressed, ])),
        cex = 0.1, las = 1, 
        xlab = "% total count per cell", 
        col = (scales::hue_pal())(50)[50:1], 
        horizontal = TRUE)


# scoring by cell cycle 
sce.all.filt = NormalizeData(sce.all.filt)
s.genes=Seurat::cc.genes.updated.2019$s.genes
g2m.genes=Seurat::cc.genes.updated.2019$g2m.genes
sce.all.filt=CellCycleScoring(object = sce.all.filt, 
                              s.features = s.genes, 
                              g2m.features = g2m.genes, 
                              set.ident = TRUE)
p4=VlnPlot(sce.all.filt, features = c("S.Score", "G2M.Score"), group.by = "orig.ident", 
           ncol = 2, pt.size = 0.1)

sce.all.filt@meta.data  %>% ggplot(aes(S.Score,G2M.Score))+geom_point(aes(color=Phase))+
  theme_minimal()

# high S.Score: S，high G2M.Score: G2M，both low: G1

###### step3: seurat ######

sce.all.filt = FindVariableFeatures(sce.all.filt)
sce.all.filt = ScaleData(sce.all.filt, 
                         vars.to.regress = c("nFeature_RNA", "percent_mito"))
sce.all.filt = RunPCA(sce.all.filt, npcs = 20)
sce.all.filt = RunTSNE(sce.all.filt, npcs = 20)
sce.all.filt = RunUMAP(sce.all.filt, dims = 1:10)

sce.all.filt <- FindNeighbors(sce.all.filt, dims = 1:15)
sce.all.filt <- FindClusters(sce.all.filt, resolution = 0.8)
# table(sce.all.filt@meta.data$RNA_snn_res.0.8)  
 
DimPlot(sce.all.filt, 
        reduction = 'umap')

genes_to_check = c('PTPRC', 'CD3D', 'CD3E', 'CD4','CD8A','CD19', 'CD79A', 'MS4A1' ,
                   'IGHG1', 'MZB1', 'SDC1',
                   'CD68', 'CD163', 'CD14', 
                   'TPSAB1' , 'TPSB2',  # mast cells,
                   'RCVRN','FPR1' , 'ITGAM' ,
                   'FGF7','MME', 'ACTA2',
                   'PECAM1', 'VWF', 
                   'EPCAM' , 'KRT19', 'PROM1', 'ALDH1A1' )

p_all_markers <- DotPlot(sce.all.filt, features = genes_to_check,
                         assay='RNA'  )  + coord_flip()


###### step4 ######

sce.all=sce.all.filt

celltype=data.frame(ClusterID=0:2,
                    celltype=0:2) 
celltype[celltype$ClusterID %in% c(2),2]='myeloid'
celltype[celltype$ClusterID %in% c(0 ),2]='epithelial' 
celltype[celltype$ClusterID %in% c( 1),2]='fibo'   

sce.all@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  sce.all@meta.data[which(sce.all@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}

# table(sce.all@meta.data$celltype)

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

# table(sce.all@meta.data$celltype,sce.all@meta.data$seurat_clusters)

DimPlot(sce.all, reduction = "umap", group.by = "celltype",label = T)  

p_all_markers=DotPlot(sce.all, features = genes_to_check,
                      assay='RNA' ,group.by = 'celltype' )  + coord_flip()
p_umap=DimPlot(sce.all, reduction = "umap", group.by = "celltype",label = T) 


###### step5 (auto cell annotation) ######


###### step6 ######

sce=sce.all
Idents(sce)=sce@meta.data$seurat_clusters
sce.markers <- FindAllMarkers(object = sce, only.pos = TRUE, min.pct = 0.25, 
                              thresh.use = 0.25)

DT::datatable(sce.markers)
pro='fluidigm'
write.csv(sce.markers,file=paste0(pro,'_sce.markers.csv'))

library(dplyr) 
top10 <- sce.markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
DoHeatmap(sce,top10$gene,size=3)

# table(sce$celltype,sample_ann$Biological_Condition)












