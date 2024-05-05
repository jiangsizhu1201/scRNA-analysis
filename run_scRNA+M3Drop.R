#load all the necessary libraries
options(warn=-1) 
suppressMessages(library(scater))
suppressMessages(library(scRNAseq))
library(ggplot2)
library(tidyr)
library(cowplot)
library("FactoMineR")
library("factoextra")
library("ROCR")


#######  scRNAseq   #######

Example data: Pollen et al. 2014
pluripotent stem cells --> neural progenitor cells (“NPC”) , “GW16” and “GW21” ，“GW21+3” 

Complete data: 23730 features，301 samples <https://hemberg-lab.github.io/scRNA.seq.datasets/human/tissues/> 

Expression matrix: RSEM (hg38 RefSeq transcriptome)

library(scRNAseq)

## ----- Load Example Data -----
data(fluidigm) 
ct <- floor(assays(fluidigm)$rsem_counts)
# ct[1:4,1:4] 
sample_ann <- as.data.frame(colData(fluidigm))
DT::datatable(sample_ann)

## phenotype/metadata
box <- lapply(colnames(sample_ann[,1:19]),function(i) {
    dat <-  sample_ann[,i,drop=F] 
    dat$sample=rownames(dat)
   ggplot(dat, aes('all cells', get(i))) +
          geom_boxplot() +
          xlab(NULL)+ylab(i)
})
plot_grid(plotlist=box, ncol=5 )

## group
pa <- colnames(sample_ann[,c(1:9,11:16,18,19)])
tf <- lapply(pa,function(i) {
 # i=pa[1]
  dat <-  sample_ann[,i]  
  dat <- abs(log10(dat))
  fivenum(dat)
  (up <- mean(dat)+2*sd(dat))
  (down <- mean(dat)- 2*sd(dat) ) 
  valid <- ifelse(dat > down & dat < up, 1,0 ) 
})

tf <- do.call(cbind,tf)
choosed_cells <- apply(tf,1,function(x) all(x==1))
table(sample_ann$Biological_Condition)
sample_ann=sample_ann[choosed_cells,]
table(sample_ann$Biological_Condition)
ct <- ct[,choosed_cells]

# ct[1:4,1:4] 
#### EDA                      
counts <- ct
fivenum(apply(counts,1,function(x) sum(x>0) ))
boxplot(apply(counts,1,function(x) sum(x>0) ))
fivenum(apply(counts,2,function(x) sum(x>0) ))
hist(apply(counts,2,function(x) sum(x>0) ))
choosed_genes=apply(counts,1,function(x) sum(x>0) )>0
# table(choosed_genes)
counts <- counts[choosed_genes,]
                    
## ----- Gene Correlation Analysis -----
                    
dat <- log2(edgeR::cpm(counts) + 1)
dat[1:4, 1:4]
dat_back <- dat

exprSet <- dat_back
# colnames(exprSet)
pheatmap::pheatmap(cor(exprSet))
group_list <- sample_ann$Biological_Condition
tmp <- data.frame(g = group_list)
rownames(tmp) <-  colnames(exprSet)
# intra-group > inter-group
pheatmap::pheatmap(cor(exprSet), annotation_col = tmp)
# dim(exprSet)
exprSet = exprSet[apply(exprSet, 1, function(x)
# sum(x > 1) > 5), ]
# dim(exprSet)

exprSet <-  exprSet[names(sort(apply(exprSet, 1, mad), decreasing = T)[1:500]), ]
M <-cor(log2(exprSet + 1))
tmp <- data.frame(g = group_list)
rownames(tmp) <-  colnames(M)
pheatmap::pheatmap(M, annotation_col = tmp)

# table(sample_ann$LibraryName)

## ----- Clustering -----
                        
dat <- dat_back
hc <- hclust(dist(t(dat))) 
plot(hc,labels = FALSE)
clus <-  cutree(hc, 4) 
group_list <-  as.factor(clus) 
# table(group_list)
# table(group_list,sample_ann$Biological_Condition)   

## ----- Dimension Reduction -----
### PCA
                        
dat <- dat_back
dat <- t(dat)
dat <- as.data.frame(dat)
plate <- sample_ann$Biological_Condition # 
dat <-  cbind(dat, plate) #

# The variable plate (index = ) is removed
# before PCA analysis
dat.pca <- PCA(dat[, -ncol(dat)], graph = FALSE)

fviz_pca_ind(
      dat.pca,
      #repel =T,
      geom.ind = "point",
      # show points only (nbut not "text")
      col.ind = dat$plate,
      # color by groups
      #palette = c("#00AFBB", "#E7B800"),
      addEllipses = TRUE,
      # Concentration ellipses
      legend.title = "Groups"
) 
                        
### TSNE
dat_matrix <- dat.pca$ind$coord
library(Rtsne) 
# dat_matrix = dat_back
# if Remove duplicates before running TSNE, check_duplicated = FALSE

# tsne_out <- Rtsne(dat_matrix,pca=FALSE,perplexity=30,theta=0.0, check_duplicates = FALSE) # Run TSNE
tsne_out <- Rtsne(dat_matrix,perplexity=10)
plate <- sample_ann$Biological_Condition # 
plot(tsne_out$Y,col= rainbow(4)[as.numeric(as.factor(plate))], pch=19) 

# tsne_out$Y is used to save time
opt_tsne=tsne_out$Y
# table(kmeans(opt_tsne,centers = 4)$clust)
plot(opt_tsne,  col=kmeans(opt_tsne,centers = 4)$clust, pch=19, xlab="tSNE dim 1", ylab="tSNE dim 2")
library(dbscan)
                        
plot(opt_tsne,  col=dbscan(opt_tsne,eps=3.1)$cluster, pch=19, xlab="tSNE dim 1", ylab="tSNE dim 2")
# table(dbscan(opt_tsne,eps=3.1)$cluster)
# compare different algorithmn
# table(kmeans(opt_tsne,centers = 4)$clust,dbscan(opt_tsne,eps=3.1)$cluster)

#######  M3Drop  #######
## ----- Create M3  Object -----
library(M3Drop) 
Normalized_data <- M3DropCleanData(counts, 
                                   labels = sample_ann$Biological_Condition , 
                                   is.counts=TRUE, min_detected_genes=2000)

### Michaelis-Menten

fits <- M3DropDropoutModels(Normalized_data$data)

# Sum absolute residuals
data.frame(MM=fits$MMFit$SAr, Logistic=fits$LogiFit$SAr,
           DoubleExpo=fits$ExpoFit$SAr) 
# Sum squared residuals
data.frame(MM=fits$MMFit$SSr, Logistic=fits$LogiFit$SSr,
           DoubleExpo=fits$ExpoFit$SSr)

## ----- DE Analysis -----
                        
DE_genes <- M3DropDifferentialExpression(Normalized_data$data, 
                                         mt_method="fdr", mt_threshold=0.01)

par(mar=c(1,1,1,1)) 
heat_out <- M3DropExpressionHeatmap(DE_genes$Gene, Normalized_data$data, 
                                    cell_labels = Normalized_data$labels)


 ## ----- Clustering + Find Marker Genes -----

cell_populations <- M3DropGetHeatmapCellClusters(heat_out, k=4)

marker_genes <- M3DropGetMarkers(Normalized_data$data, cell_populations)
# table(cell_populations,Normalized_data$labels)

### marker genes
# head(marker_genes[marker_genes$Group==4,],20) 
marker_genes[rownames(marker_genes)=="FOS",] 

### visualization
par(mar=c(1,1,1,1)) 
choosed_marker_genes=as.character(unlist(lapply(split(marker_genes,marker_genes$Group), function(x) (rownames(head(x,20))))))
heat_out <- M3DropExpressionHeatmap(choosed_marker_genes, Normalized_data$data, cell_labels =  cell_populations)

