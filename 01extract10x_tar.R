### 提供的是10X格式的标准三个文件，选择下载数据之后需要对数据进行整理，将三个文件分别整理到对应的文件夹中。

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
