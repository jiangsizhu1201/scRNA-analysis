# Read files
library(Seurat)
library(ggplot2)
library(data.table)
library(Matrix)

## Seurat V5
dir = 'dirname' # include samples(each with varcods, features, and matrix)
sc_data <- Read10x(dir)

sc_obj(CreateSeuratObject(counts = sc_data, 
                            project = "obj_name", #
                            min.cells = 3, # 
                            min.features = 200) # 
##### csv.gz file

sc_data <- fread('file.csv.gz') # data.table
sc_data <- as.data.frame(sc_data)
## set row names
rownames(sc_data) <- sc_data$V1
sc_data$V1 <- NULL

data.matrix = Matrix(as.matrix(sc_data), sparse = T) # ddiMatrix
data.matrix = as(as.matrix(sc_data),'dgCMatrix') # dgCMatrix

sc_obj <- CreateSeuratObject(counts = data.matrix, 
                             project = "obj_name") #
