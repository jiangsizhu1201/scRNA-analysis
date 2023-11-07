library(scRNAseq)

fluidigm <- ReprocessedFluidigmData()
fluidigm

out <- listDatasets()
fs=out$Call
fs

lapply( fs , function(x){
  pro= gsub("'",'-', gsub('[()]','', x ))
  print(pro)
  f=paste0( pro ,
            '.Rdata')
  if( ! file.exists(f)){
    sce=eval(parse(text=x))
    save(sce,file = f)
  }

})
