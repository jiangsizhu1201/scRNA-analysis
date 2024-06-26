suppressMessages(library(scater))
suppressMessages(library(scRNAseq))
library(ggplot2)
library(tidyr)
library(cowplot)
library(monocle)
library("FactoMineR")
library("factoextra")
library("ROCR")

## ----- Example Data -----

###### step1 ######

fluidigm <- ReprocessedFluidigmData()
ct <- floor(assays(fluidigm)$rsem_counts)
# ct[1:4,1:4] 
sample_ann <- as.data.frame(colData(fluidigm))
DT::datatable(sample_ann)

gene_ann <- data.frame(
  gene_short_name = rownames( ct ) , 
  row.names =  rownames( ct ) 
)

pd <- new("AnnotatedDataFrame",
          data=sample_ann)
fd <- new("AnnotatedDataFrame",
          data=gene_ann) 

###### step2: monocle ######

# create object + basic dimension reduction

sc_cds <- newCellDataSet(
  as.matrix(ct), 
  phenoData = pd,
  featureData =fd,
  expressionFamily = negbinomial.size(),
  lowerDetectionLimit=1)
sc_cds # CellDataSet 

# standard monocle processing

sc_cds <- detectGenes(sc_cds, min_expr = 1) 
sc_cds <- sc_cds[fData(sc_cds)$num_cells_expressed > 10, ]

cds <- sc_cds
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds) 

# select genes + find clusters

disp_table <- dispersionTable(cds)
unsup_clustering_genes <- subset(disp_table,
                                 mean_expression >= 0.1)

cds <- setOrderingFilter(cds, unsup_clustering_genes$gene_id)

plot_ordering_genes(cds) 
plot_pc_variance_explained(cds, 
                           return_all = F) # norm_method='log'

cds <- reduceDimension(cds, max_components = 2, num_dim = 6, # num_dim determined by PCA 
                       reduction_method = 'tSNE', verbose = T)
cds <- clusterCells(cds, num_clusters = 6) 
plot_cell_clusters(cds, 1, 2 )

###### step3: trajectory ######

###### step3.1: monocle ######

# check the trajectory for certain biological condition

# table(pData(cds)$Cluster)
# table(pData(cds)$Cluster,pData(cds)$Biological_Condition)
# plot_cell_clusters(cds, 1, 2 )

# Note: use the seurat_clusters, not monocle results

# Biological conditions

pData(cds)$Cluster=pData(cds)$Biological_Condition
# table(pData(cds)$Cluster)

Sys.time()
diff_test_res <- differentialGeneTest(cds,
                                      fullModelFormulaStr = "~Cluster")
Sys.time()

# Select genes that are significant at an FDR < 10%
sig_genes <- subset(diff_test_res, qval < 0.1)
sig_genes=sig_genes[order(sig_genes$pval),]
# head(sig_genes[,c("gene_short_name", "pval", "qval")] ) 
cg=as.character(head(sig_genes$gene_short_name)) 

plot_genes_jitter(cds[cg,],
                  grouping = "Cluster",
                  color_by = "Cluster",
                  nrow= 3,
                  ncol = NULL )

cg2=as.character(tail(sig_genes$gene_short_name)) 

# Find the most significant genes for visualization

plot_genes_jitter(cds[cg2,],
                  grouping = "Cluster",
                  color_by = "Cluster",
                  nrow= 3,
                  ncol = NULL )

###### step3.2: trejactory ######

# step 3.2.1: select DEGs
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
cds <- setOrderingFilter(cds, ordering_genes)
plot_ordering_genes(cds) 

# step 3.2.2: dimension reduction
cds <- reduceDimension(cds, max_components = 2,
                       method = 'DDRTree') # default by Monocle2

# step 3.2.3: order cells
cds <- orderCells(CDs)

# step 3.2.4: visualization
plot_cell_trajectory(cds, color_by = "Cluster")  

plot_genes_in_pseudotime(cds[cg,],
                         color_by = "Cluster") 

phe=pData(cds)
boxplot(phe$Pseudotime,phe$Cluster)

# https://davetang.org/muse/2017/10/01/getting-started-monocle/

# step 4: find genes that change along the trajectory

my_cds_subset=cds
# pseudotime is now a column in the phenotypic data as well as the cell state

# differentialGeneTest - time consuming
my_pseudotime_de <- differentialGeneTest(my_cds_subset,
                                         fullModelFormulaStr = "~sm.ns(Pseudotime)",
                                         cores = 1 )
# cluster needed : cores = 4
# if mac: cores = 1 

###### step4: Visualization of Monocle ######

library(Seurat)
library(gplots)
library(ggplot2)
library(monocle)
library(dplyr)
load(file = 'output_of_monocle.Rdata')

cds=my_cds_subset
phe=pData(cds)

library(ggsci)

p1=plot_cell_trajectory(cds, color_by = "Cluster")  + scale_color_nejm() 
plot_cell_trajectory(cds, color_by = "Biological_Condition")  
p2=plot_cell_trajectory(cds, color_by = "Pseudotime")  
p3=plot_cell_trajectory(cds, color_by = "State")  + scale_color_npg()

library(patchwork)
p1+p2/p3

phe=pData(cds)
# table(phe$State,phe$Cluster) 

library(dplyr)
my_pseudotime_de %>% arrange(qval) %>% head() 

# test the top genes
my_pseudotime_de %>% arrange(qval) %>% head() %>% select(gene_short_name) -> my_pseudotime_gene
my_pseudotime_gene=my_pseudotime_gene[,1]
plot_genes_in_pseudotime(my_cds_subset[my_pseudotime_gene,])+ scale_color_npg()


plot_genes_jitter(my_cds_subset[my_pseudotime_gene,],
                  grouping = "Cluster",
                  color_by = "Cluster",
                  nrow= 3,
                  ncol = NULL )+ scale_color_nejm()


# cluster the top 50 genes that vary as a function of pseudotime
my_pseudotime_de %>% arrange(qval) %>% head(50) %>% select(gene_short_name) -> gene_to_cluster
gene_to_cluster <- gene_to_cluster[,1]

# table(pData(my_cds_subset)$Cluster,pData(my_cds_subset)$State) 
ac=pData(my_cds_subset)[c('Biological_Condition','State','Pseudotime')]

# Pseudotime: 100col x 50 row

my_pseudotime_cluster <- plot_pseudotime_heatmap(my_cds_subset[gene_to_cluster,],
                                                 # num_clusters = 2, 
                                                 # add_annotation_col = ac,
                                                 show_rownames = TRUE,
                                                 return_heatmap = TRUE)




 

