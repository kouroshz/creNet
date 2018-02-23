## Given gene expression data, this function applies a MS5 filter to the probes
## and normalizes the data
##
## Arguments:
##
## celphath       -    path to the directory containing the .CEL files
##
## phenotypeFile  -    tab delimited file of phenotype of CEL files.
##
## eset           -    output normalized gene expression data
##
## Note: In this version probes are filtered with MS5 method.
normalizeChipMas5 = function(celpath, phenotypeFile){

  phenotypes = read.table(phenotypeFile, header=F, sep="\t", stringsAsFactors = F)

  if(length(grep(".CEL",  phenotypes$V1)) != nrow(phenotypes)){
    phenotypes$V1 = paste(phenotypes$V1, ".CEL", sep ='')
  }

  print('Reading CEL files')
  abatch = ReadAffy(celfile.path=celpath)
  ## Removing uninformative probes.
  x = mas5calls(abatch)
  y = exprs(x)
  absent = rowSums(y == 'A' | y == 'M')
  absent = which(absent == ncol(y))

  ## Normalizing the cel files
  if(length(grep(".CEL.gz",  row.names(pData(abatch)))) > 0){
    phenotypes$V1 = paste(gsub('.CEL', '', phenotypes$V1))
    phenotypes$V1 = paste(phenotypes$V1, ".CEL.gz", sep ='')
  }else if(length(grep(".CEL",  row.names(pData(abatch)))) > 0){
    if(length(grep(".CEL",  phenotypes$V1)) != nrow(phenotypes)){
      phenotypes$V1 = paste(phenotypes$V1, ".CEL", sep ='')
    }
  }

  tmp = merge(row.names(pData(abatch)), phenotypes, by.x = 1, by.y = 1)
  if(ncol(phenotypes) == 3){
    pData(abatch) = as.data.frame(tmp[,2:3])
    colnames(pData(abatch)) = c('Cohort', 'Response')
  }else{
    pData(abatch) = as.data.frame(tmp[,2])
    colnames(pData(abatch)) = c('Response')
  }
  row.names(pData(abatch)) = tmp[,1]


  ## Normalize the data
  print('Normalizing')
  eset = rma(abatch)


  eset = eset[-absent,]

  AFFX = grep('AFF*',rownames(exprs(eset)))
  if(length(AFFX[1]) > 0){
    eset = eset[-AFFX,]
  }

  eset
}
