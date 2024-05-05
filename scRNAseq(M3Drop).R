
if (!requireNamespace("Rtsne"))
    install.packages("Rtsne")
if (!requireNamespace("FactoMineR"))
    install.packages("FactoMineR")
if (!requireNamespace("factoextra"))
    install.packages("factoextra")
if (!requireNamespace("BiocManager"))
    install.packages("BiocManager")
if (!requireNamespace("scater"))
    BiocManager::install("scater")
if (!requireNamespace("scRNAseq"))
    BiocManager::install("scRNAseq") 
if (!requireNamespace("M3Drop"))
    BiocManager::install("M3Drop") 
if (!requireNamespace("ROCR"))
    BiocManager::install("ROCR") 


rm(list = ls()) # clear the environment
#load all the necessary libraries
options(warn=-1) # turn off warning message globally
suppressMessages(library(scater))
suppressMessages(library(scRNAseq))
library(ggplot2)
library(tidyr)
library(cowplot)
library("FactoMineR")
library("factoextra")
library("ROCR")


##  scRNAseq

Pollen et al. 2014
pluripotent stem cells --> neural progenitor cells (“NPC”) , “GW16” and “GW21” ，“GW21+3” 


library(scRNAseq)
## ----- Load Example Data -----
data(fluidigm) 
ct <- floor(assays(fluidigm)$rsem_counts)
ct[1:4,1:4] 
sample_ann <- as.data.frame(colData(fluidigm))
DT::datatable(sample_ann)

## phenotype

box <- lapply(colnames(sample_ann[,1:19]),function(i) {
    dat <-  sample_ann[,i,drop=F] 
    dat$sample=rownames(dat)
   ggplot(dat, aes('all cells', get(i))) +
          geom_boxplot() +
          xlab(NULL)+ylab(i)
})
plot_grid(plotlist=box, ncol=5 )
# ggsave(file="stat_all_cells.pdf")


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

ct[1:4,1:4] 
counts <- ct
fivenum(apply(counts,1,function(x) sum(x>0) ))
boxplot(apply(counts,1,function(x) sum(x>0) ))
fivenum(apply(counts,2,function(x) sum(x>0) ))
hist(apply(counts,2,function(x) sum(x>0) ))
choosed_genes=apply(counts,1,function(x) sum(x>0) )>0
table(choosed_genes)
counts <- counts[choosed_genes,]

## gene corr


dat <- log2(edgeR::cpm(counts) + 1)
dat[1:4, 1:4]
dat_back <- dat

exprSet <- dat_back
colnames(exprSet)
pheatmap::pheatmap(cor(exprSet))
group_list <- sample_ann$Biological_Condition
tmp <- data.frame(g = group_list)
rownames(tmp) <-  colnames(exprSet)
# intra-group > inter-group
pheatmap::pheatmap(cor(exprSet), annotation_col = tmp)
dim(exprSet)
exprSet = exprSet[apply(exprSet, 1, function(x)
sum(x > 1) > 5), ]
dim(exprSet)
 
dim(exprSet)
exprSet <-  exprSet[names(sort(apply(exprSet, 1, mad), decreasing = T)[1:500]), ]
dim(exprSet)
M <-cor(log2(exprSet + 1))
tmp <- data.frame(g = group_list)
rownames(tmp) <-  colnames(M)
pheatmap::pheatmap(M, annotation_col = tmp)

table(sample_ann$LibraryName)



dat <- dat_back
hc <- hclust(dist(t(dat))) 
plot(hc,labels = FALSE)
clus <-  cutree(hc, 4) 
group_list <-  as.factor(clus) 
table(group_list)
table(group_list,sample_ann$Biological_Condition)   

dat <- dat_back
dat <- t(dat)
dat <- as.data.frame(dat)
plate <- sample_ann$Biological_Condition # 
dat <-  cbind(dat, plate) #
dat[1:4, 1:4]
table(dat$plate)

# The variable plate (index = ) is removed
# before PCA analysis
dat.pca <- PCA(dat[, -ncol(dat)], graph = FALSE)
head(dat.pca$var$coord) ## 
head(dat.pca$ind$coord) ## 
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


dat_matrix <- dat.pca$ind$coord
# Set a seed if you want reproducible results
set.seed(42)
library(Rtsne) 
# dat_matrix = dat_back
# if Remove duplicates before running TSNE, check_duplicated = FALSE

# tsne_out <- Rtsne(dat_matrix,pca=FALSE,perplexity=30,theta=0.0, check_duplicates = FALSE) # Run TSNE
tsne_out <- Rtsne(dat_matrix,perplexity=10)
plate <- sample_ann$Biological_Condition # 
plot(tsne_out$Y,col= rainbow(4)[as.numeric(as.factor(plate))], pch=19) 

# tsne_out$Y is used to save time
head(tsne_out$Y)
opt_tsne=tsne_out$Y
table(kmeans(opt_tsne,centers = 4)$clust)
plot(opt_tsne,  col=kmeans(opt_tsne,centers = 4)$clust, pch=19, xlab="tSNE dim 1", ylab="tSNE dim 2")
library(dbscan)
plot(opt_tsne,  col=dbscan(opt_tsne,eps=3.1)$cluster, pch=19, xlab="tSNE dim 1", ylab="tSNE dim 2")
table(dbscan(opt_tsne,eps=3.1)$cluster)
# compare different algorithmn
table(kmeans(opt_tsne,centers = 4)$clust,dbscan(opt_tsne,eps=3.1)$cluster)





## M3Drop


library(M3Drop) 
Normalized_data <- M3DropCleanData(counts, 
                                   labels = sample_ann$Biological_Condition , 
                                   is.counts=TRUE, min_detected_genes=2000)
dim(Normalized_data$data)
length(Normalized_data$labels)
class(Normalized_data)
str(Normalized_data)



### Michaelis-Menten

fits <- M3DropDropoutModels(Normalized_data$data)

# Sum absolute residuals
data.frame(MM=fits$MMFit$SAr, Logistic=fits$LogiFit$SAr,
           DoubleExpo=fits$ExpoFit$SAr) 
# Sum squared residuals
data.frame(MM=fits$MMFit$SSr, Logistic=fits$LogiFit$SSr,
           DoubleExpo=fits$ExpoFit$SSr)

DE_genes <- M3DropDifferentialExpression(Normalized_data$data, 
                                         mt_method="fdr", mt_threshold=0.01)
dim(DE_genes)
head(DE_genes)


par(mar=c(1,1,1,1)) 
heat_out <- M3DropExpressionHeatmap(DE_genes$Gene, Normalized_data$data, 
                                    cell_labels = Normalized_data$labels)


 

cell_populations <- M3DropGetHeatmapCellClusters(heat_out, k=4)
library("ROCR") 
marker_genes <- M3DropGetMarkers(Normalized_data$data, cell_populations)
table(cell_populations,Normalized_data$labels)


### marker genes
head(marker_genes[marker_genes$Group==4,],20) 
marker_genes[rownames(marker_genes)=="FOS",] 

par(mar=c(1,1,1,1)) 
choosed_marker_genes=as.character(unlist(lapply(split(marker_genes,marker_genes$Group), function(x) (rownames(head(x,20))))))
heat_out <- M3DropExpressionHeatmap(choosed_marker_genes, Normalized_data$data, cell_labels =  cell_populations)

