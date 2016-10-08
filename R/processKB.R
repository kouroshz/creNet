processKB <- function(ents.file, rels.file, verbose = F)
{
  ## Generating the one-level network form the knowledge base
  ## Note: There are ' and # characters in KB that mess up the tables. 
  ## The following will clean them up.
  ents <- read.table(ents.file, header = TRUE, stringsAsFactors = FALSE, strip.white=TRUE, sep = '\t',
                     quote = NULL, comment.char = '')
  colnames(ents) = c('uid', 'name', 'id', 'type')
  ents <- na.omit(ents) # remove incomplete entries
  
  rels <- read.table(rels.file, header = TRUE, stringsAsFactors = FALSE, strip.white=TRUE, sep = '\t', 
                     quote = NULL, comment.char = '')
  colnames(rels) = c('uid', 'srcuid', 'trguid', 'type', 'pmids', 'nls')
  
  cleanup <- function(x) {
    x <- gsub('\"', '', x)
    x <- gsub('\'', 'p', x)
    x <- gsub('#', '_', x)
    x <- gsub(' ', '', x)
    return(x) 
  } 
  
  ents <- as.data.frame(lapply(ents,cleanup), stringsAsFactors = F)
  rownames(ents) <- 1:nrow(ents)
  rels <- as.data.frame(lapply(rels,cleanup), stringsAsFactors = F)
  rownames(rels) = 1:nrow(rels)
  
  ## Convert uids to integers
  uid.orig = ents$uid
  uid.new = seq(1, length(ents$uid))
  id.map = data.frame(uid.orig = uid.orig, uid.new = uid.new, stringsAsFactors = F)
  ents$uid = uid.new
  rels$srcuid = uid.new[match(rels$srcuid, uid.orig)]
  rels$trguid = uid.new[match(rels$trguid, uid.orig)]
  
  
  ## Protein, Compound or mRNA entries
  ents = unique(ents[which(ents$type %in% c('mRNA', 'Protein', 'Compound')),])
  
  ## (unique) Protein/Compound entries
  ents.pc = unique(ents[which(ents$type %in% c('Protein', 'Compound')),])
  ## (unique) mRNA entries
  ents.mRNA = unique(ents[which(ents$type == 'mRNA'),])
  
  ## src has to be protein or compound
  rels = rels[which(rels$srcuid %in% ents.pc$uid),]
  rownames(rels) = 1:nrow(rels)
  ## target has to be mRNA
  rels = rels[which(rels$trguid %in% ents.mRNA$uid),]
  rownames(rels) = 1:nrow(rels)
  
  ## remove -1's if any
  ents = ents[(ents$id != -1 | is.na(ents$uid)),]
  rels = rels[which(rels$trguid %in% ents$uid & rels$srcuid %in% ents$uid), ] ##update rels
  ents = ents[which(ents$uid %in% c(rels$srcuid, rels$trguid)), ]
  rownames(ents) = 1:nrow(ents)
  
  ## Remove duplicates
  if (anyDuplicated(ents)) ents <- ents[!duplicated(ents),]
  if (anyDuplicated(rels)) rels <- rels[!duplicated(rels),]
  
  ## Check that relation types are valid
  valid <- (rels$type %in% c("increase","decrease","conflict"))
  if (!all(valid)) rels <- rels[valid,]
  

  if (anyDuplicated(rels[,c("srcuid","trguid")])) {
    ## Identify duplicated rows
    dup.rels = rels[duplicated(rels[,c(2,3)]), ]
    ## Take one example from each
    dup.rels.uniq = dup.rels[!duplicated(dup.rels[,c(2,3)]), ]
    
    for(i in 1:nrow(dup.rels.uniq)){
      ind = which(rels$srcuid == dup.rels.uniq$srcuid[i] & rels$trguid == dup.rels.uniq$trguid[i])
      if(all(c("increase","decrease") %in% unique(rels$type[ind]))){
        rels$type[ind] = 'conflict'
      }else if("increase" %in% unique(rels$type[ind])){
        rels$type[ind] = 'increase'
      }else if("decrease" %in% unique(rels$type[ind])){
        rels$type[ind] = 'decrease'
      }else{
        rels$type[ind] = 'conflict'
      }
    }
    rels$type[which(!(rels$type %in% c('increase', 'decrease', 'conflict')))] = 'conflict'
  }
  
  ## remove duplicated rels
  if (any(duplicated(rels[,c("srcuid","trguid", "type")]))) 
    rels <- rels[!duplicated(rels[,c("srcuid","trguid", "type")]),]
  
  ## Check that sources in rels are protein/compound and targets are mRNA
  # Note: this also takes care of removing lines in rels with uids not matched in ents
  ind.pc <- which(ents$type %in% c("Protein","Compound"))
  inds <- which(!(rels$srcuid %in% ents$uid[ind.pc]))
  if (length(inds)) rels <- rels[-inds,]
  ind.mRNA <- which(ents$type=="mRNA")
  inds <- which(!(rels$trguid %in% ents$uid[ind.mRNA]))
  if (length(inds)) rels <- rels[-inds,]
  
  ## Consolidate uids in ents	(remove duplicates)
  ## KZ: preordering the factors. split does not preserve the order
  ##uid <- split(ents$uid,list(ents$id,ents$type)) # group uids by id and type (mRNA, protein, compound...)
  uid <- split(ents$uid,list(factor(ents$id, levels = unique(ents$id)),factor(ents$type, levels = unique(ents$type))))
  muid <- sapply(uid,length) # multiplicities
  uid <- uid[muid>1] # duplicates are present and must be removed only if multiplicity = 1,
  dups <- lapply(uid, function(x) matrix(c(rep(x[1],length(x)-1),x[-1]), ncol=2))
  dups <- do.call(rbind,dups) 
  colnames(dups) <- c("orig","dups")	
  # matrix with 2 columns containing "good" uids in 1st column 
  # and duplicated uids in 2nd
  ents <- ents[!(ents$uid %in% dups[,2]),] # remove duplicates
  
  ## Map srcuid & trguid in rels to consolidated uids in ents
  ind.loc <- which(rels$trguid %in% dups[,2])
  if (length(ind.loc)) {
    ind.val <- match(rels$trguid[ind.loc],dups[,2])
    rels$trguid[ind.loc] <- dups[ind.val,1]
  }
  ind.loc <- which(rels$srcuid %in% dups[,2])
  if (length(ind.loc)) {
    ind.val <- match(rels$srcuid[ind.loc],dups[,2])
    rels$srcuid[ind.loc] <- dups[ind.val,1]
  }
  
  
  if (anyDuplicated(rels[,c("srcuid","trguid")])) {
    ## Identify duplicated rows
    dup.rels = rels[duplicated(rels[,c(2,3)]), ]
    ## Take one example from each
    dup.rels.uniq = dup.rels[!duplicated(dup.rels[,c(2,3)]), ]
    
    for(i in 1:nrow(dup.rels.uniq)){
      ind = which(rels$srcuid == dup.rels.uniq$srcuid[i] & rels$trguid == dup.rels.uniq$trguid[i])
      if(all(c("increase","decrease") %in% unique(rels$type[ind]))){
        rels$type[ind] = 'conflict'
      }else if("increase" %in% unique(rels$type[ind])){
        rels$type[ind] = 'increase'
      }else if("decrease" %in% unique(rels$type[ind])){
        rels$type[ind] = 'decrease'
      }else{
        rels$type[ind] = 'conflict'
      }
    }
    rels$type[which(!(rels$type %in% c('increase', 'decrease', 'conflict')))] = 'conflict'
  }
  
  ## remove duplicated rels
  if (any(duplicated(rels[,c("srcuid","trguid", "type")]))) 
    rels <- rels[!duplicated(rels[,c("srcuid","trguid", "type")]),]
  
  if (verbose == TRUE) {
    cat("\n Processed network dimensions:")
    cat("\n ents:", dim(ents)[1])
    cat("\n rels:", dim(rels)[1])
  }
  
  L = list(ents = ents, rels = rels, id.map = id.map)
  
}

processKBOrig <- function(ents.file, rels.file, verbose = F)
{
  
  ## Generating the one-level network form the knowledge base
  ## Note: There are ' and # characters in KB that mess up the tables. The following will clean them up.
  ents = read.table(ents.file, header = T, stringsAsFactors = F, strip.white=T, sep = '\t', quote = NULL,
                    comment.char = '')
  colnames(ents) = c('uid', 'name', 'id', 'type')
  rels = read.table(rels.file, header = T, stringsAsFactors = F, strip.white=T, sep = '\t', quote = NULL,
                    comment.char = '')
  colnames(rels) = c('uid', 'srcuid', 'trguid', 'type', 'pmids', 'nls')
  
  ents = data.frame(apply(ents, 2, function(x) gsub('\"', '', x)), stringsAsFactors = F)
  ents = data.frame(apply(ents, 2, function(x) gsub('\'', 'p', x)), stringsAsFactors = F)
  ents = data.frame(apply(ents, 2, function(x) gsub('#', '_', x)), stringsAsFactors = F)
  ents = data.frame(apply(ents, 2, function(x) gsub(' ', '', x)), stringsAsFactors = F)
  
  rels = data.frame(apply(rels, 2, function(x) gsub('\"', '', x)), stringsAsFactors = F)
  rels = data.frame(apply(rels, 2, function(x) gsub('\'', 'p', x)), stringsAsFactors = F)
  rels = data.frame(apply(rels, 2, function(x) gsub('#', '_', x)), stringsAsFactors = F)
  rels = data.frame(apply(rels, 2, function(x) gsub(' ', '', x)), stringsAsFactors = F)
  
  rownames(rels) = 1:nrow(rels)
  
  ## Convert uids to integers
  uid.orig = ents$uid
  uid.new = seq(1, length(ents$uid))
  id.map = data.frame(uid.orig = uid.orig, uid.new = uid.new, stringsAsFactors = F)
  ents$uid = uid.new
  rels$srcuid = uid.new[match(rels$srcuid, uid.orig)]
  rels$trguid = uid.new[match(rels$trguid, uid.orig)]
  
  L = list(ents = ents, rels = rels, id.map = id.map)
  
  L
}
