mergeMultProbs <- function(x, method = c("var","cv","sum"))
{
  
  method <- match.arg(method)
  
  ## Remove constant columns
  x <- remConstCols(x, verbose = FALSE) 
  
  n <- ncol(x)
  grps <- split(1:n,colnames(x))
  ord <- match(unique(colnames(x)),names(grps))
  grps <- grps[ord]
  if (length(grps) == n) return(x)
  mergeFun <- switch(method, var= var, cv=function(x,na.rm) sd(x,na.rm)/mean(x,na.rm=na.rm), sum=sum)
  stat <- sapply(grps, function(g) if (length(g)==1) 0 else apply(x[,g],2,mergeFun,na.rm=TRUE))
  keep <- lapply(stat,function(g) which.max(g))
  ind <- unlist(lapply(1:length(keep), function(i) grps[[i]][keep[[i]]]))
  x[,ind]	
}

mergeMultProbs2 <- function(x, method = c("var","cv","sum"))
{
	method <- match.arg(method)
	n <- ncol(x)
	grps <- split(1:n,colnames(x))
	ord <- match(unique(colnames(x)),names(grps))
	grps <- grps[ord]
	if (length(grps) == n) return(x)
	mergeFun <- switch(method, var= var, cv=function(x,na.rm) sd(x,na.rm)/mean(x,na.rm=na.rm), sum=sum)
	stat <- sapply(grps, function(g) if (length(g)==1) 0 else apply(x[,g],2,mergeFun,na.rm=TRUE))
	ind <- sapply(stat,which.max(stat))
	x[,ind]	
}
