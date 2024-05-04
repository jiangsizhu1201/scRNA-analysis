###### step1 ######

suppressMessages(library(scater))
suppressMessages(library(scRNAseq))
library(ggplot2)
library(tidyr)
library(cowplot)
library(monocle)
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

gene_ann <- data.frame(
  gene_short_name = rownames( ct ) , 
  row.names =  rownames( ct ) 
)
head(gene_ann)

pd <- new("AnnotatedDataFrame",
          data=sample_ann)
fd <- new("AnnotatedDataFrame",
          data=gene_ann) 
ct[1:4,1:4]


###### step2: monocle ######

sc_cds <- newCellDataSet(
  as.matrix(ct), 
  phenoData = pd,
  featureData =fd,
  expressionFamily = negbinomial.size(),
  lowerDetectionLimit=1)
sc_cds

# CellDataSet 

library(monocle)
sc_cds
sc_cds <- detectGenes(sc_cds, min_expr = 1) 
sc_cds <- sc_cds[fData(sc_cds)$num_cells_expressed > 10, ]
sc_cds

cds <- sc_cds
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds) 
# plyr 

# select genes
disp_table <- dispersionTable(cds)
unsup_clustering_genes <- subset(disp_table,
                                 mean_expression >= 0.1)
unsup_clustering_genes
cds <- setOrderingFilter(cds, 
                         unsup_clustering_genes$gene_id)
plot_ordering_genes(cds) 
plot_pc_variance_explained(cds, 
                           return_all = F) # norm_method='log'
# num_dim determined by PCA
cds <- reduceDimension(cds, max_components = 2, num_dim = 6,
                       reduction_method = 'tSNE', verbose = T)
cds <- clusterCells(cds, num_clusters = 6) 
plot_cell_clusters(cds, 1, 2 )
table(pData(cds)$Cluster) 
colnames(pData(cds)) 
table(pData(cds)$Biological_Condition)
table(pData(cds)$Cluster,pData(cds)$Biological_Condition) 
plot_cell_clusters(cds, 1, 2 )

save(cds,file = 'input_cds.Rdata')



###### step3: trajectory ######

### monocle ### 
cds 

table(pData(cds)$Cluster)
table(pData(cds)$Cluster,pData(cds)$Biological_Condition)
plot_cell_clusters(cds, 1, 2 )

# Biological conditions
pData(cds)$Cluster=pData(cds)$Biological_Condition
table(pData(cds)$Cluster)

Sys.time()
diff_test_res <- differentialGeneTest(cds,
                                      fullModelFormulaStr = "~Cluster")
Sys.time()

# Select genes that are significant at an FDR < 10%
sig_genes <- subset(diff_test_res, qval < 0.1)
sig_genes=sig_genes[order(sig_genes$pval),]
head(sig_genes[,c("gene_short_name", "pval", "qval")] ) 
cg=as.character(head(sig_genes$gene_short_name)) 

plot_genes_jitter(cds[cg,],
                  grouping = "Cluster",
                  color_by = "Cluster",
                  nrow= 3,
                  ncol = NULL )
cg2=as.character(tail(sig_genes$gene_short_name)) 
plot_genes_jitter(cds[cg2,],
                  grouping = "Cluster",
                  color_by = "Cluster",
                  nrow= 3,
                  ncol = NULL )



# select gene set (significant)

ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
ordering_genes
cds <- setOrderingFilter(cds, ordering_genes)
plot_ordering_genes(cds) 

cds <- reduceDimension(cds, max_components = 2,
                       method = 'DDRTree')

cds <- orderCells(cds)

plot_cell_trajectory(cds, color_by = "Cluster")  
ggsave('monocle_cell_trajectory_for_seurat.pdf')

length(cg)
plot_genes_in_pseudotime(cds[cg,],
                         color_by = "Cluster") 
ggsave('monocle_plot_genes_in_pseudotime_for_seurat.pdf')

phe=pData(cds)
boxplot(phe$Pseudotime,phe$Cluster)

# https://davetang.org/muse/2017/10/01/getting-started-monocle/

my_cds_subset=cds
# pseudotime is now a column in the phenotypic data as well as the cell state
head(pData(my_cds_subset))
# differentialGeneTest - time consuming
my_pseudotime_de <- differentialGeneTest(my_cds_subset,
                                         fullModelFormulaStr = "~sm.ns(Pseudotime)",
                                         cores = 1 )
# cluster needed : cores = 4
# if mac: cores = 1 
head(my_pseudotime_de)
save( my_cds_subset,my_pseudotime_de,
      file = 'output_of_monocle.Rdata')


###### step4 ######

library(Seurat)
library(gplots)
library(ggplot2)
library(monocle)
library(dplyr)
load(file = 'output_of_monocle.Rdata')

cds=my_cds_subset
phe=pData(cds)
colnames(phe)
library(ggsci)
p1=plot_cell_trajectory(cds, color_by = "Cluster")  + scale_color_nejm() 
p1
ggsave('trajectory_by_cluster.pdf')
plot_cell_trajectory(cds, color_by = "Biological_Condition")  

p2=plot_cell_trajectory(cds, color_by = "Pseudotime")  
p2
ggsave('trajectory_by_Pseudotime.pdf')

p3=plot_cell_trajectory(cds, color_by = "State")  + scale_color_npg()
p3
ggsave('trajectory_by_State.pdf')
library(patchwork)
p1+p2/p3

phe=pData(cds)
head(phe)
table(phe$State,phe$Cluster) 

library(dplyr)
my_pseudotime_de %>% arrange(qval) %>% head() 
# save the top 6 genes
my_pseudotime_de %>% arrange(qval) %>% head() %>% select(gene_short_name) -> my_pseudotime_gene
my_pseudotime_gene=my_pseudotime_gene[,1]
my_pseudotime_gene
plot_genes_in_pseudotime(my_cds_subset[my_pseudotime_gene,])+ scale_color_npg()
ggsave('monocle_top6_pseudotime_by_state.pdf')

plot_genes_jitter(my_cds_subset[my_pseudotime_gene,],
                  grouping = "Cluster",
                  color_by = "Cluster",
                  nrow= 3,
                  ncol = NULL )+ scale_color_nejm()
ggsave('monocle_top6_pseudotime_by_cluster.pdf')


# cluster the top 50 genes that vary as a function of pseudotime
my_pseudotime_de %>% arrange(qval) %>% head(50) %>% select(gene_short_name) -> gene_to_cluster
gene_to_cluster <- gene_to_cluster[,1]
gene_to_cluster
colnames(pData(my_cds_subset))
table(pData(my_cds_subset)$Cluster,pData(my_cds_subset)$State) 
ac=pData(my_cds_subset)[c('Biological_Condition','State','Pseudotime')]
head(ac)
# Pseudotime: 100col x 50 row

my_pseudotime_cluster <- plot_pseudotime_heatmap(my_cds_subset[gene_to_cluster,],
                                                 # num_clusters = 2, 
                                                 # add_annotation_col = ac,
                                                 show_rownames = TRUE,
                                                 return_heatmap = TRUE)
my_pseudotime_cluster

pdf('monocle_top50_heatmap.pdf')
print(my_pseudotime_cluster)
dev.off()



 

