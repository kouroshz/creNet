cvGlmnet <- function(x.train, y.train, x.test, y.test, num.iter, nfold = 5, alpha = 1){
  
  n <- nrow(x.test)
  m <- ncol(x.test)
  
  test.probs  <- matrix(NA, nrow = n, ncol = num.iter)
  labs        <- matrix(NA, nrow = n, ncol = num.iter)
  best.lams   <- matrix(NA, nrow = num.iter, ncol = 2)

  nonzero.genes = {}
  for(iter in 1:num.iter){
    if(verbose){
      cat('\n\n')
      cat(paste('running iteration', iter))
      cat('\n')
    }
    ## perform one cv on the whole data to get the best lambda
    #fit <- cv.glmnet(as.matrix(train.x), as.matrix(train.y), family = "binomial", type.measure="deviance", nfolds = dim(train.x)[1], grouped = F)
    #lambda = fit$lambda
    #lambda.min = fit$lambda.min
    cv.fit <- NULL
    while(is.null(cv.fit)){
      try(cv.fit<-cv.glmnet(x.train, y.train, family = "binomial", type.measure = "deviance", nfolds = nfold))
    }
    
    #cv.fit <- cv.glmnet(x.train, y.train, family = "binomial", type.measure = "deviance", nfolds = nfold)
    
    lambdas         <- cv.fit$lambda
    lam.min         <- cv.fit$lambda.min
    best.lam.ind    <- which(lambdas == lam.min)
    best.lams[iter, 1] <- best.lam.ind
    best.lams[iter, 2] <- lam.min
 
    ## Fitting the Full model
    main.fit <- glmnet(x.train, y.train, family="binomial", alpha = 1, lambda = lambdas)
    best.betas      <- main.fit$beta[,best.lam.ind]
    best.intercepts <- main.fit$a0[best.lam.ind]
    
    nonzero.genes = c(nonzero.genes, colnames(x.train[,which(best.betas != 0)]))
    
    #x.test  <- standardize(x.test)$x
    test.probs[,iter] <- predict(main.fit, newx = x.test, s = lam.min, type = "response", mode = "lambda")
    labs[,iter] <- ifelse(test.probs[,iter] > 0.5,1,0)
  }  
  L <- list(fit = main.fit, best.lams = best.lams, test.probs = test.probs, labs = labs, nonzero.genes = nonzero.genes)
  
  return(L)
}
