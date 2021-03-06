processDataAndKB <- function(ents.file, rels.file, data.train.file, data.test.file = NULL, 
                             verbose = FALSE, uids = NULL, isLasso = FALSE){
  
  
  x.test <- NULL
  y.test <- NULL
  ## Rading the data
  if(verbose)
    cat('\n reading the training data\n')
  
  data.train <- read.table(data.train.file, header = T, sep = '\t', stringsAsFactors = F, 
                           check.names = F)
  
  if(!is.null(data.test.file)){
    if(verbose)
      cat('\n reading the testing data\n')
    
    data.test <- read.table(data.test.file, header = T, sep = '\t', stringsAsFactors = F, 
                            check.names = F)
    
    ## Match train and test files so they have the same columns
    L <- matchTrainTest(data.train, data.test)
    
    x.train = data.matrix(L$train.dat[,2:ncol(L$train.dat)])
    y.train = as.vector(L$train.dat[,1])
    x.test = data.matrix(L$test.dat[,2:ncol(L$test.dat)])
    y.test = as.vector(L$test.dat[,1])
  }else{
    x.train = data.matrix(data.train[,2:ncol(data.train)])
    y.train = as.vector(data.train[,1])
  }
  
  
  if(isLasso & is.null(uids)){
    if(!is.null(data.test.file)){
      D <- rbind(x.train, x.test)
    }else{
      D <- x.train
    }
    
    ## merge genes with same entrez ids.
    D <- mergeMultProbs(D, method = 'var') 
    
    x.train <- D[1:nrow(x.train), ]
    if(!is.null(data.test.file)){
      x.test  <- D[(nrow(x.train)+1):(nrow(x.train)+nrow(x.test)), ]
    }
    ents <- NULL
    rels <- NULL
  }else{
    if(verbose)
      cat('\n processing the network\n')
    L <- processKB(ents.file, rels.file)
    ents <- L$ents
    rels <- L$rels
    id.map <- L$id.map
    if(!is.null(uids)){
      ind <- match(uids, L$id.map$uid.orig)
      if(is.na(ind[1])){
        cat('\n Your specified uids do not exist in the network')
        stop()
      }
      uids <- L$id.map$uid.new[ind]
      L <- nodeNetList(uids, ents, rels, levels = F)
      ents <- L$ents
      rels <- L$rels
    }
    
    if(!is.null(data.test.file)){
      D <- rbind(x.train, x.test)
    }else{
      D <- x.train
    }
    
    ## merge genes with same entrez ids.
    D <- mergeMultProbs(D, method = 'var') 
    
    D <- matchDataToKb(ents, rels, D)
    ents <- D$ents
    rels <- D$rels
    x.train <- D$D.matched[1:nrow(x.train), ]
    if(!is.null(data.test.file)){
      x.test  <- D$D.matched[(nrow(x.train)+1):(nrow(x.train)+nrow(x.test)), ]
    }
  }
  
  L <- list(ents = ents, rels = rels, x.train = x.train, x.test = x.test, 
            y.train = y.train, y.test = y.test)
}
