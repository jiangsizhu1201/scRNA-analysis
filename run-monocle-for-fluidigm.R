###### step1: 获取单细胞转录组数据集 ######

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


###### step2: 构建monocle的对象，并且进行基础的降维聚类分群 ######

# 构建对象，seurat，monocle，scater 都类似
# monocle标准流程，降维聚类分群 

sc_cds <- newCellDataSet(
  as.matrix(ct), 
  phenoData = pd,
  featureData =fd,
  expressionFamily = negbinomial.size(),
  lowerDetectionLimit=1)
sc_cds

# 这个 CellDataSet 对象一定要认识清楚，务必花两个小时去摸索它。


# 接下来仅仅是  monocle的标准流程而已
library(monocle)
sc_cds
sc_cds <- detectGenes(sc_cds, min_expr = 1) 
# 数值可以自行摸索
sc_cds <- sc_cds[fData(sc_cds)$num_cells_expressed > 10, ]
sc_cds

cds <- sc_cds
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds) 
# plyr 

# 并不是所有的基因都有作用，所以先进行挑选，合适的基因用来进行聚类。
disp_table <- dispersionTable(cds)
unsup_clustering_genes <- subset(disp_table,
                                 mean_expression >= 0.1)
unsup_clustering_genes
cds <- setOrderingFilter(cds, 
                         unsup_clustering_genes$gene_id)
plot_ordering_genes(cds) 
plot_pc_variance_explained(cds, 
                           return_all = F) # norm_method='log'
# 其中 num_dim 参数选择基于上面的PCA图
cds <- reduceDimension(cds, max_components = 2, num_dim = 6,
                       reduction_method = 'tSNE', verbose = T)
cds <- clusterCells(cds, num_clusters = 6) 
plot_cell_clusters(cds, 1, 2 )
table(pData(cds)$Cluster) 
colnames(pData(cds)) 
table(pData(cds)$Biological_Condition)
table(pData(cds)$Cluster,pData(cds)$Biological_Condition) 
plot_cell_clusters(cds, 1, 2 )

# 可以看到 monocle 给细胞重新定义了亚群，亚群数量是自己选择的
# 整体来说，monocle和seurat 各自独立流程定义的亚群的一致性还不错。

# 只是跑流程而已
save(cds,file = 'input_cds.Rdata')



###### step3: 使用monocle的拟时序分析######

### 然后查看monocle ### 
cds 
# 接下来很重要，到底是看哪个性状的轨迹
table(pData(cds)$Cluster)
table(pData(cds)$Cluster,pData(cds)$Biological_Condition)
plot_cell_clusters(cds, 1, 2 )

## 我们这里并不能使用 monocle的分群
# 还是依据前面的 seurat分群, 也就是说前面的代码仅仅是流程而已，我们没有使用那些结果哦

# 其实取决于自己真实的生物学意图
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
#  挑选差异最显著的基因可视化
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

# 前面是找差异基因，后面是做拟时序分析

# 第一步: 挑选合适的基因. 有多个方法，例如提供已知的基因集，
# 这里选取统计学显著的差异基因列表
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
ordering_genes
cds <- setOrderingFilter(cds, ordering_genes)
plot_ordering_genes(cds) 
# 第二步: 降维。降维的目的是为了更好的展示数据。函数里提供了很多种方法,
# 不同方法的最后展示的图都不太一样, 其中“DDRTree”是Monocle2使用的默认方法
cds <- reduceDimension(cds, max_components = 2,
                       method = 'DDRTree')
# 第三步: 对细胞进行排序
cds <- orderCells(cds)
# 最后两个可视化函数，对结果进行可视化
plot_cell_trajectory(cds, color_by = "Cluster")  
ggsave('monocle_cell_trajectory_for_seurat.pdf')

length(cg)
plot_genes_in_pseudotime(cds[cg,],
                         color_by = "Cluster") 
ggsave('monocle_plot_genes_in_pseudotime_for_seurat.pdf')

phe=pData(cds)
boxplot(phe$Pseudotime,phe$Cluster)

# https://davetang.org/muse/2017/10/01/getting-started-monocle/
# 前面根据差异基因，推断好了拟时序，也就是说把差异基因动态化了

# 后面就可以具体推断哪些基因随着拟时序如何的变化
my_cds_subset=cds
# pseudotime is now a column in the phenotypic data as well as the cell state
head(pData(my_cds_subset))
# 这个differentialGeneTest会比较耗费时间
my_pseudotime_de <- differentialGeneTest(my_cds_subset,
                                         fullModelFormulaStr = "~sm.ns(Pseudotime)",
                                         cores = 1 )
# 不知道为什么在Mac电脑无法开启并行计算了 ，不过我测试了在Windows 电脑设置cores = 4是可以的
# 如果你是Mac电脑，自己修改 cores = 1 即可 
head(my_pseudotime_de)
save( my_cds_subset,my_pseudotime_de,
      file = 'output_of_monocle.Rdata')


###### step4: 可视化monocle的拟时序分析######

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
# 这个热图绘制的并不是纯粹的细胞基因表达量矩阵，而是被 Pseudotime 好了的100列，50行的矩阵

my_pseudotime_cluster <- plot_pseudotime_heatmap(my_cds_subset[gene_to_cluster,],
                                                 # num_clusters = 2, 
                                                 # add_annotation_col = ac,
                                                 show_rownames = TRUE,
                                                 return_heatmap = TRUE)
my_pseudotime_cluster

pdf('monocle_top50_heatmap.pdf')
print(my_pseudotime_cluster)
dev.off()



 

