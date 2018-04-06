#### Functions for using selectiveInference package for lasso and lars analyses

# Function to standardize covariates
# Used only for continuous covariates, not for factors/binary covariates
standardize <- function(covar){
  mean.covar <- mean(covar, na.rm = TRUE)
  sd.covar <- sd(covar[!is.na(covar)])
  return((covar - mean.covar)/sd.covar)
}


#### Function to convert factor columns to numeric columns
factor2numeric <- function(factorx){
  numx <- as.numeric(levels(factorx))[factorx]
  return(numx)
}

#### Function to convert all factor variables in a data frame to numeric
allFactor2numeric <- function(data){
  factorNames <- data %>% Filter(f = is.factor) %>% names()
  for(i in 1:length(factorNames)){
    name.i <- factorNames[i]
    data[,name.i] <- factor2numeric(data[,name.i])
  }
  return(data)
}

########################################################################################
#### Function to get results from fixedLassoInf()
##########################################################################################

getLassoInfResults <- function(y = ylasso, x = xlasso, alpha = alphaBest, lambda = lamddaBest){
  require(glmnet)
  require(selectiveInference)
  # y = response variable vector
  # x = matrix of covariates
  assign("y", y, envir = .GlobalEnv)
  assign("x", x, envir = .GlobalEnv)
  # alpha = selected value of alpha from cross-validation
  # lambda = selected value of lambda from cross-validation
  ### Run glmnet::glmnet using bestAlpha
  bestModel <- glmnet::glmnet(y = y, x = x, alpha = alphaBest, standardize = FALSE, thresh = 1E-20)
  ### Extract best coefficients ("betas") using bestLambda divided by sample size (n), and remove intercept
  ### Based on description in help page for fixedLassoInf()
  n <- nrow(x)
  betas <- coef(bestModel, s=bestLambda, exact = TRUE, y = y, x = x)[-1]
  ### Estimate the variance
  sigmaEst <- estimateSigma(x = x, y = y, standardize = FALSE)
  ### Run fixedLassoInf() to get p-values for each selected variable
  lassoTest <- fixedLassoInf(x = x, y = y, beta = betas, lambda = bestLambda, alpha = 0.05, type = "partial", sigma = sigmaEst$sigmahat)
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


########################################################################################
#### Function to fit lar model and get p-values/CI 
##########################################################################################

fitLAR <- function(y = ylasso, x = xlasso, infType = "aic"){
  require(dplyr)
  require(selectiveInference)
  # y = response variable vector
  # x = matrix of covariates
  # infType is an option passed to larInf() function
  assign("y", y, envir = .GlobalEnv)
  assign("x", x, envir = .GlobalEnv)
  assign("infType", infType, envir = .GlobalEnv)
  # Fit lar model
  larfit <- lar(x = x, y = y, normalize = FALSE)
  sigmaEst <- estimateSigma(x = x, y = y, standardize = FALSE)
  # Estimate p-values, using AIC stopping criterion
  larTest <- larInf(larfit, alpha = 0.05, type = infType, sigma = sigmaEst$sigmahat)
  # Clean up the results and combine with covariate names
  varNames <- attr(x, "dimnames")[[2]]
  larResults <- data.frame(varNames = varNames[larTest$vars],
                           varIndex = larTest$vars,
                           coef = coef(larfit, s = larTest$khat+1, mode = "step")[larTest$vars],
                           pvalue = larTest$pv,
                           cil = larTest$ci[,1],
                           ciu = larTest$ci[,2])
  # Return model fits and results
  outList <- list(larfit = larfit, # model fit object from lar()
                  larTest = larTest, # object from larInf()
                  larResults = larResults) # clean(er) table of lar results
  return(outList)
}


