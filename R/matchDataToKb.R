matchDataToKb <- function(ents, rels, D)
{
  
  gene.ids <- colnames(D)
  ents.mRNA <- ents[which(ents$id %in% gene.ids & ents$type == 'mRNA'),]
  
  ## First filter the data
  ind <- which(colnames(D) %in% ents.mRNA$id)
  unmatched.ind <- setdiff(1:length(colnames(D)), ind)
  
  D.matched <- D[, ind]
  D.unmatched <- D[, unmatched.ind]
  gene.ids = colnames(D.matched)
  ##leftout.genes <- colnames(D.unmatched)
  
  ## Prune the KB so it includes genes that are present in the data
  
  ents.mRNA <- ents.mRNA[which(ents.mRNA$id %in% gene.ids),]   ##update mRNA.ents
  ents <- ents[((ents$type == 'mRNA' & ents$uid %in% ents.mRNA$uid) | ents$type != 'mRNA'), ]
  rels <- rels[which((rels$trguid %in% ents.mRNA$uid) & (rels$srcuid %in% ents$uid)), ] ##update rels
  ents <- ents[which(ents$uid %in% c(rels$srcuid, rels$trguid)), ]
  
  L = list(ents = ents, rels = rels, D.matched = D.matched, D.unmatched = D.unmatched)
  L
}

matchDataToKb2 <- function(ents, rels, D)
{
	
	## If there are lines in relation file where source and target uid 
	## are not matched in entry file, remove them
	ind <- (rels$trguid %in% ents$uid & rels$srcuid %in% ents$uid)
	if (!all(ind)) rels <- rels[,ind]

	## If there are genes present in KB but not in data, remove them
	gene.ids <- colnames(D)
	ind.ents.out <- which(ents$type == "mRNA" & !(ents$id %in% gene.ids))
	if (length(ind.ents.out)) { 
		## Prune the entry file
		uid.out <- ents$uid[ind.ents.out]
		ents <- ents[-ind.ents.out,]
	 	## Prune the relation file
		ind.rels.out <- which(rels$trguid %in% uid.out)
		if (length(ind.rels.out)) rels <- rels[-ind.rels.out,] 
	} 
		
	## Determine matched and unmatched portion of data in relation file 
	ind.matched <- which(gene.ids %in% ents$id[ents$type == "mRNA"])
	n <- length(gene.ids)
	if (length(ind.matched) == n) {
		D.matched <- D
		D.unmatched <- NULL
	} else { 
		D.matched <- D[,ind.matched]
		D.unmatched <- D[,-ind.matched]
	}
	
	L = list(ents = ents, rels = rels, D.matched = D.matched, D.unmatched = D.unmatched)
	L
}
