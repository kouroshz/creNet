normalizeCEL <- function(celPath, phenoFile, annoFile, verbose = TRUE)
{

	eset = suppressWarnings(normalizeChipMas5(celPath, phenoFile))

	if (verbose == TRUE)
		cat("\n Annotations:", eset@annotation)

	annotations = read.table(annoFile,header = T, sep = '\t', stringsAsFactors = F, quote="", comment.char="#")

	annotations = annotations[,c('ID', 'ENTREZ_GENE_ID')]
	colnames(annotations) = c('id', 'entrez')
	grps = unique(pData(eset)$Cohort)

	train.eset = eset
	x.dat = data.frame(id = rownames(exprs(train.eset)), exprs(train.eset), stringsAsFactors = F)
	x.dat = merge(x.dat, annotations, by.x = 1, by.y = 1)

	ind.mult.entz = which(unlist(lapply(x.dat$entrez, function(x) length(strsplit(x, split = ' ')[[1]]))) != 1)
	x.dat = x.dat[-ind.mult.entz, ]

	rm.cols = {}
	dup.genes <- unique(x.dat$entrez[which(duplicated(x.dat$entrez))])
	for(ez in dup.genes){
	  ind = which(x.dat$entrez == ez)
	  vars <- apply(x.dat[ind, 2:(ncol(x.dat) - 1)], 1, var)
	  rm.cols <- c(rm.cols, ind[-which.max(vars)])
	}

	x.dat <- x.dat[-rm.cols, ]
	ids <- x.dat$entrez
	x.dat <- t(x.dat[,2:(ncol(x.dat) - 1)])
	colnames(x.dat) <- ids

	return(x.dat)
}
