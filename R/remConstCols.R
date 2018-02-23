remConstCols <- function(x, verbose = TRUE) 
{
  if (!is.matrix(x)) dim(x) <- c(length(x),1)
  sdev <- apply(x,2,sd)
  const <- which(sdev < 1e-8)
  if (verbose)
      cat(paste('\n', length(const), 'columns removed\n'))
  if (length(const) == 0) 
  	return(x) else return(x[,-const])
}
