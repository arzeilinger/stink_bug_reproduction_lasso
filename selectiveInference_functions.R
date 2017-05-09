#### Functions for using selectiveInference package for lasso and lars analyses


#### Function to convert factor columns to numeric columns
factor2numeric <- function(factorx){
  numx <- as.numeric(levels(factorx))[factorx]
  return(numx)
}


########################################################################################
#### Function to get results from fixedLassoInf()
##########################################################################################

getLassoInfResults <- function(y = ylasso, x = xlasso, alpha = alphaBest, lambda = lamddaBest){
  require(glmnet)
  require(selectiveInference)
  # y = response variable vector
  # x = matrix of covariates
  # alpha = selected value of alpha from cross-validation
  # lambda = selected value of lambda from cross-validation
  ### Run glmnet::glmnet using bestAlpha
  bestModel <- glmnet::glmnet(y = ylasso, x = xlasso, alpha = alphaBest, standardize = FALSE)
  ### Extract best coefficients ("betas") using bestLambda divided by sample size (n), and remove intercept
  ### Based on description in help page for fixedLassoInf()
  n <- nrow(xlasso)
  betas <- coef(bestModel, s=bestLambda/n, exact = TRUE, y = ylasso, x = xlasso)[-1]
  ### Estimate the variance
  sigmaEst <- estimateSigma(x = xlasso, y = ylasso, standardize = FALSE)
  ### Run fixedLassoInf() to get p-values for each selected variable
  lassoTest <- fixedLassoInf(x = xlasso, y = ylasso, beta = betas, lambda = bestLambda, alpha = 0.05, type = "partial", sigma = sigmaEst$sigmahat)
  # Put results in nice data.frame with covariate names
  lassoTestResults <- data.frame(varName = names(lassoTest$vars),
                                 varIndex = as.numeric(lassoTest$vars),
                                 coef = lassoTest$coef0,
                                 pvalue = lassoTest$pv,
                                 cil = lassoTest$ci[,1],
                                 ciu = lassoTest$ci[,2])
  outlist <- list(bestModel, lassoTest, lassoTestResults)
  names(outlist) <- c("bestModel", "lassoTest", "lassoTestResults")
  return(outlist)
}

