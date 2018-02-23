getGeneVals <- function(uids, evidence, as.factor = TRUE)
{
  val = numeric(length(uids))
  ind.loc = which(uids %in% evidence$uid)
  if (length(ind.loc)) {
  	ind.val = match(uids[ind.loc],evidence$uid)
  	val[ind.loc] = evidence$val[ind.val]
  }
  if (as.factor) val <- factor(val, levels = c(-1,0,1))
  val
}
