fs=list.files('./','features')
fs

samples1= gsub('.tsv.gz','',gsub('features.','',fs))
samples1

samples2 = samples1


lapply(1:length(samples2), function(i){
  x=samples2[i]
  y=fs[i]
  dir.create(x,recursive = T)
  file.copy(from= y ,
            to=file.path(x,  'features.tsv.gz' )) 
  file.copy(from= gsub('features','matrix',gsub('tsv','mtx',y)),
            to= file.path(x, 'matrix.mtx.gz' ) ) 
  file.copy(from= gsub('features','barcodes',y),
            to= file.path(x, 'barcodes.tsv.gz' )) 
  
})