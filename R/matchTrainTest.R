matchTrainTest <- function(train.dat, test.dat){

  comm <- intersect(colnames(train.dat)[2:ncol(train.dat)], colnames(test.dat)[2:ncol(test.dat)])

  ind1 <- which(colnames(train.dat)[2:ncol(train.dat)] %in% comm) + 1
  train.dat <- train.dat[,c(1, ind1)]

  ind2 <- which(colnames(test.dat)[2:ncol(test.dat)] %in% comm) + 1
  test.dat <- test.dat[,c(1, ind2)]

  ind <- match(colnames(train.dat)[2:ncol(train.dat)], colnames(test.dat)[2:ncol(test.dat)]) + 1

  test.dat <- test.dat[,c(1, ind)]

  ##all(colnames(train.dat) == colnames(test.dat))

  return(L<- list(train.dat = train.dat, test.dat = test.dat))
}