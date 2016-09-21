getDEGs <- function(ents, x.train, y.train, de.sig = 0.05, verbose=TRUE)
{
  
  ents.mRNA = ents[which(ents$type == 'mRNA'), ]
  
  ## Compute pvalues of differential expression for all genes.
  if (verbose == TRUE) {
    cat("\n\n Computing Differentially Expressed genes and generating evidence file")
  }
  Group.train = factor(y.train , levels = levels(as.factor(y.train)))
  design.train = model.matrix(~Group.train)
  fit.train = limma::lmFit(t(x.train),design.train)
  fit.train = limma::eBayes(fit.train)
  tab.train = limma::topTable(fit.train, coef = 2, adjust = "fdr", n = ncol(x.train))
  if(!('ID' %in% colnames(tab.train))){
    tab.train = data.frame(ID=rownames(tab.train), tab.train, stringsAsFactors=F)
    rownames(tab.train) = 1:nrow(tab.train)
  }
  inds = match(colnames(x.train), tab.train[,1])
  
  ## Generate the evidence file for training examples
  evidence = tab.train[inds,c(1, 2, 5)]
  evidence = evidence[which(evidence$P.Value < de.sig & abs(evidence$logFC) > 0.55),c(1,2)]
  ##evidence = evidence[which(evidence$adj.P.Value < de.sig & abs(evidence$logFC) > 0.55),c(1,2)]

  if(nrow(evidence) < 20 & de.sig < 0.05){
    if(verbose){
      print('Number of DEGs too small.')
      print('Adjusting the de.seq cutoff to 0.05')
    }
    evidence = tab.train[,c(1, 2, 5)]
    evidence = evidence[which(evidence$P.Value < 0.05 & abs(evidence$logFC) > 0.55),c(1,2)]
  }
  
  
  evidence[,2] = ifelse(evidence[,2] > 0, 1, -1)
  evidence = data.frame(id = suppressWarnings(as.numeric(evidence[,1])), val = as.numeric(evidence[,2]),
                        stringsAsFactors = F)
  evidence = evidence[!is.na(evidence[,1]),]
  evidence = evidence[!duplicated(evidence[,1]), ]
  
  n.e1 = nrow(evidence)
  ## Make sure evidence id is in the ents. Same gene with multiple ids may exist
  ## in the original ents.
  evidence = evidence[which(evidence[,1] %in% ents.mRNA$id),]
  n.e2 = nrow(evidence)
  if (verbose == TRUE)
    cat("\n", (n.e1-n.e2), "evidence removed!")
  
  ## Change id back to uid
  evidence.tmp = merge(evidence, ents.mRNA, by.x = 1, by.y = 3)
  if (verbose == TRUE && nrow(evidence.tmp) > nrow(evidence)) {
    cat("\n Warning")
    cat("\n Entry ids in evidence file mapped to multiple uids in KB")
    cat("\n Total number of duplicates:" , (nrow(evidence.tmp) - nrow(evidence)))
    cat("\n Selecting one uid at random")
  }
  evidence = data.frame(uid = evidence.tmp$uid, val = evidence.tmp$val, stringsAsFactors = F)
  evidence = evidence[!duplicated(evidence),]
  
  if (verbose == TRUE)
    cat("\n Dimension of input file:\n", dim(evidence))
  
  
  L = list(evidence = evidence, tab.train = tab.train)
  
  L
}

getDEGsOrig <- function(ents, rels, x.train, y.train, de.sig = 0.05, verbose=TRUE)
{
  
  ents.mRNA = ents[which(ents$type == 'mRNA'), ]
  
  ## Compute pvalues of differential expression for all genes.
  if (verbose == TRUE) {
    cat("\n\n Computing Differentially Expressed genes and generating evidence file")
  }
  Group.train = factor(y.train , levels = levels(as.factor(y.train)))
  design.train = model.matrix(~Group.train)
  fit.train = limma::lmFit(t(x.train),design.train)
  fit.train = limma::eBayes(fit.train)
  tab.train = limma::topTable(fit.train, coef = 2, adjust = "fdr", n = ncol(x.train))
  if(!('ID' %in% colnames(tab.train))){
    tab.train = data.frame(ID=rownames(tab.train), tab.train, stringsAsFactors=F)
    rownames(tab.train) = 1:nrow(tab.train)
  }
  inds = match(colnames(x.train), tab.train[,1])
  
  ## Generate the evidence file for training examples
  evidence = tab.train[inds,c(1, 2, 5)]
  evidence = evidence[which(evidence$P.Value < de.sig & abs(evidence$logFC) > 0.55),c(1,2)]
  if(nrow(evidence) < 20 & de.sig < 0.05){
    if(verbose){
      print('Number of DEGs too small.')
      print('Adjusting the de.seq cutoff to 0.05')
    }
    evidence = tab.train[,c(1, 2, 5)]
    evidence = evidence[which(evidence$P.Value < 0.05 & abs(evidence$logFC) > 0.55),c(1,2)]
  }
  
  
  evidence[,2] = ifelse(evidence[,2] > 0, 1, -1)
  evidence = data.frame(id = suppressWarnings(as.numeric(evidence[,1])), val = as.numeric(evidence[,2]),
                        stringsAsFactors = F)
  evidence = evidence[!is.na(evidence[,1]),]
  evidence = evidence[!duplicated(evidence[,1]), ]
  
  n.e1 = nrow(evidence)
  ## Make sure evidence id is in the ents. Same gene with multiple ids may exist
  ## in the original ents.
  evidence = evidence[which(evidence[,1] %in% ents.mRNA$id),]
  n.e2 = nrow(evidence)
  if (verbose == TRUE)
    cat("\n", (n.e1-n.e2), "evidence removed!")
  
  ## Change id back to uid
  evidence.tmp = merge(evidence, ents.mRNA, by.x = 1, by.y = 3)
  if (verbose == TRUE && nrow(evidence.tmp) > nrow(evidence)) {
    cat("\n Warning")
    cat("\n Entry ids in evidence file mapped to multiple uids in KB")
    cat("\n Total number of duplicates:" , (nrow(evidence.tmp) - nrow(evidence)))
    cat("\n Selecting one uid at random")
  }
  evidence = data.frame(uid = evidence.tmp$uid, val = evidence.tmp$val, stringsAsFactors = F)
  evidence = evidence[!duplicated(evidence),]
  
  if (verbose == TRUE)
    cat("\n Dimension of input file:\n", dim(evidence))
  
  
  L = list(evidence = evidence, tab.train = tab.train)
  
  L
}

getDEGs2 <- function(ents, x.train, y.train, de.sig = 0.05, verbose=TRUE)
{

  ## Compute pvalues of differential expression for all genes.
  if (verbose) {
    cat("\n\n Computing Differentially Expressed genes and generating evidence file")
  }
  
  ## Merge columns associated with same gene (multiple probes)
  ##x.train <- mergeMultProbs(x.train)
  
  ## Fit linear model and get fdr-adjusted p-values
  Group.train = factor(y.train) # , levels = levels(as.factor(y.train)))
  design.train = model.matrix(~Group.train)
  fit.train = limma::lmFit(t(x.train),design.train)
  fit.train = limma::eBayes(fit.train)
  tab.train = limma::topTable(fit.train, coef = 2, adjust = "fdr", n = ncol(x.train))
  ## Add ID column as needed
  if(!('ID' %in% colnames(tab.train))){
    tab.train$ID = rownames(tab.train)
    rownames(tab.train) = 1:nrow(tab.train)
  }


  ## Generate the evidence file for training examples
  inds <- match(colnames(x.train), tab.train$ID)
  evidence <- tab.train[inds,c("ID","logFC","P.Value")]
  inds <- which(evidence$P.Value < de.sig & abs(evidence$logFC) > 0.55)
  if (length(inds) < 20 & de.sig < 0.05) {
    if (verbose) 
      cat("\n Number of DEGs too small \n Adjusting the de.seq cutoff to 0.05")
    inds <- which(tab.train$P.Value < 0.05 & abs(tab.train$logFC) > 0.55)
  }
  evidence <- evidence[inds,]  
  evidence$ID <- as.numeric(evidence$ID)
  evidence$logFC = ifelse(evidence$logFC > 0, 1, -1)
  evidence <- evidence[,c("ID","logFC","P.Value")]
  names(evidence) <- c("id","val","pval")

  ## Add unique id to evidence
  inds.mRNA <- which(ents$type == "mRNA")
  inds.matched <- match(evidence$id,ents$id[inds.mRNA])
  n <- nrow(evidence)
  matched <- !is.na(inds.matched)
  uid <- integer(n)
  uid[matched] <- ents$uid[inds.mRNA][inds.matched]
  uid[!matched] <- NA
  evidence$uid <- uid
  evidence$matched <- matched

  if (verbose) {
	if (any(!matched)) 
		cat("\n Warning \n",sum(!matched),"entry ids in evidence file not matched in KB")
    cat("\n Dimension of input file:\n", dim(evidence))
  }

  L = list(evidence = evidence, tab.train = tab.train)

  L
}
