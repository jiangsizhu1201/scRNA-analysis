cran_packages <- c('Seurat',
                   'tidyverse'
                   ) 
Biocductor_packages <- c('sva',
                         'monocle',
                         'GOplot',
                         'GSVA',
                         'plotmo',
                         'regplot',
                         'scRNAseq',
                         'BiocStyle',
                         'celldex',
                         'SingleR',
                         'BiocParallel',
                         'scater'
)

for (pkg in cran_packages){
  if (! require(pkg,character.only=T) ) {
    install.packages(pkg,ask = F,update = F)
    require(pkg,character.only=T) 
  }
}


for (pkg in Biocductor_packages){
  if (! require(pkg,character.only=T) ) {
    BiocManager::install(pkg,ask = F,update = F)
    require(pkg,character.only=T) 
  }
}

for (pkg in c(Biocductor_packages,cran_packages)){
  require(pkg,character.only=T) 
}


