rm(list=ls())
library(Seurat)
library(preprocessCore)
library(survival)
library(survminer) 
library(ggstatsplot) 
library(remotes)
library(UCell)

load( file = 'batch_cox_results.Rdata')
cox_results = as.data.frame(cox_results)
cox_results=cox_results[order(cox_results$HR,decreasing = T),]

load('scRNA_for_scAB_Scissor.Rdata')

p1=DimPlot(scRNA, reduction ="umap", group.by="celltype",label = T)


cox_results=cox_results[rownames(cox_results)%in% rownames(scRNA),]
cox_markers=list(
  pos = head(rownames(cox_results),100),
  neg = tail(rownames(cox_results),100)
)

sc_dataset <- AddModuleScore_UCell(scRNA, 
                                   features = cox_markers) 
signature.names <- paste0(names(cox_markers), "_UCell") 

options(repr.plot.width=6, repr.plot.height=4)

VlnPlot(sc_dataset, features = signature.names, 
        #group.by = "celltype", 
        stack=TRUE ) + NoLegend()

FeaturePlot(sc_dataset,'pos_UCell')

fivenum(sc_dataset$pos_UCell)
b1=table(sc_dataset$pos_UCell> fivenum(sc_dataset$pos_UCell)[4],
      sc_dataset$celltype)
gplots::balloonplot(b1)
phe=sc_dataset@meta.data

load(file = 'output_of_Scissor.Rdata')  
sc_dataset$pos_UCell = phe$pos_UCell
FeaturePlot(sc_dataset,'pos_UCell',split.by ='scissor' )
sc_dataset$neg_UCell = phe$neg_UCell
FeaturePlot(sc_dataset,'neg_UCell',split.by ='scissor' )

#Idents(sc_dataset)= sc_dataset$scissor 
deg_Scissor_markers <- FindMarkers(sc_dataset, ident.1 = "Scissor+", 
                           group.by = 'scissor', 
                           logfc.threshold = 0.25 )

deg_Scissor_markers <- deg_Scissor_markers[which(deg_Scissor_markers$p_val_adj<0.05),]

cox_markers$up = rownames(deg_Scissor_markers)[deg_Scissor_markers$avg_log2FC>0]
cox_markers$down = rownames(deg_Scissor_markers)[deg_Scissor_markers$avg_log2FC < 0]

require("VennDiagram")
grid.newpage()
venn.plot <- venn.diagram(cox_markers , NULL,
                          fill=c("red", "blue",'green','black'),
                          alpha=c(0.5,0.5,0.5,0.5), cex = 2, cat.fontface=4,
                          category.names= names(cox_markers),
                          main="venn.diagram")
grid.draw(venn.plot)


