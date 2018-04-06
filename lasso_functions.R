#### Functions for Elastic Net/LASSO analysis of stink bug reproduction

# Function to standardize covariates
# Used only for continuous covariates, not for factors/binary covariates
standardize <- function(covar){
  mean.covar <- mean(covar, na.rm = TRUE)
  sd.covar <- sd(covar[!is.na(covar)])
  return((covar - mean.covar)/sd.covar)
}


# Function to cross-validate alpha values in elastic net procedure
lassoAlphaCV <- function(x){
  lassoCVmin <- cv.glmnet(y = ylasso, x = xlasso, family = "gaussian", standardize = TRUE, alpha = x)$lambda.min
  return(c(x, lassoCVmin))
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

#### Function to run elastic net cross-validation (ENCV) multiple times for a given response variable
#### To get a distribution of the best alpha and corresponding estimates of lambda
#### Inputs are:
# data = data set; should be left as explVars
# y = response variable name
# times = number of CV runs
# alphaValues = vector of values of alpha to evaluate
#### Need to make a folder/directory with the name of the response variable before hand
#### Future version should make directory automatically

encv <- function(y = "lnlambda", times = 2000, data = explVars, alphaValues = seq(0.8,1,by=0.01)){
  #### Get data structured correctly to submit for cross-validation
  # Vector of response variable
  ylasso <- data[,y]
  assign("ylasso", ylasso, envir = .GlobalEnv)
  # Make as a matrix
  xlasso <- as.matrix(data[,-which(names(data) == y)])
  assign("xlasso", xlasso, envir = .GlobalEnv)
  # Output directory
  outdir <- paste("output/", y, "/", sep = "")
  ######################################################################################
  #### Elastic Net Cross-validation
  #### When alpha = 1, glmnet runs a lasso penalty; when 0 < alpha < 1, it uses the elastic net lasso penalty
  #### To cross-validate alpha, run cv.glmnet over a range of alpha, using lassoAlphaCV function
  alphaVec <- rep(alphaValues,times) # sequence of alpha values
  alphaCVResults <- as.data.frame(t(sapply(alphaVec, lassoAlphaCV, simplify = TRUE))) # run through alpha values
  names(alphaCVResults) <- c("alpha", "lambda.min")
  pdf(paste(outdir,"alpha-lambda_cross-validation_raw_plot.pdf", sep=""))
    plot(x = alphaCVResults$alpha, y = jitter(alphaCVResults$lambda.min))
  dev.off()
  ## summary of alpha values
  alphaSummary <- alphaCVResults %>% group_by(alpha) %>% summarise(mean = mean(lambda.min),
                                                                   median = median(lambda.min),
                                                                   sd = sd(lambda.min))
  pdf(paste(outdir,"alpha-lambda_cross-validation_median_plot.pdf", sep=""))
    plot(x = alphaSummary$alpha, y = alphaSummary$median, type = "b")
  dev.off()
  write.csv(as.data.frame(alphaSummary), file = paste(outdir,"alpha-lambda_cross-validation_summary.csv", sep=""), row.names = FALSE)
  #### Extract lambda.min values from all runs at best alpha
  alphaBest <- as.numeric(alphaSummary[alphaSummary$median == min(alphaSummary$median),"alpha"])
  lambdaBest <- as.numeric(min(alphaSummary$median))
  lambdas <- alphaCVResults[alphaCVResults$alpha == alphaBest,"lambda.min"]
  cvlist <- list(alphaCVResults = alphaCVResults,
                 alphaBest = alphaBest,
                 lambdaBest = lambdaBest,
                 # Summary of lambda values at best alpha
                 alphaSummary = alphaSummary[alphaSummary$median == min(alphaSummary$median),],
                 # Vector of lambda values at best alpha value
                 lambdas = lambdas)
  return(cvlist)
} 


########################################################################################
#### Function to run elastic net analysis for each combination alpha and lambda, estimated from cross-validation

ElasticNetFunction <- function(y = "lnlambda", data = explVars, alphaBest = alphaBest, lambdas = lambdas){
  #### Get data structured correctly to submit for cross-validation
  # Vector of response variable
  ylasso <- data[,y]
  assign("ylasso", ylasso, envir = .GlobalEnv)
  # Make covariate data as a matrix
  xlasso <- as.matrix(data[,-which(names(data) == y)])
  assign("xlasso", xlasso, envir = .GlobalEnv)
  # Output directory
  outdir <- paste("output/", y, "/", sep = "")
  #### Get model coefficients for every value of lambda at best alpha
  bsElasticNet <- glmnet(y = ylasso, x = xlasso, family = "gaussian", alpha = alphaBest, standardize = TRUE)
  # glmnet for best alpha and lambda
  #coef(bsElasticNet, s = lambdaBest)
  # Combine coefficients into a matrix and get mean coefficient for each term
  coefResults <- matrix(NA, nrow = ncol(xlasso)+1, ncol = length(lambdas))
  enPredict <- enResiduals <- matrix(NA, nrow = nrow(xlasso), ncol = length(lambdas))
  for(i in 1:length(lambdas)){
    lambda.i <- lambdas[i]
    coefResults[,i] <- as.numeric(coef(bsElasticNet, s = lambda.i))
    enPredict[,i] <- predict(bsElasticNet, s = lambda.i, newx = xlasso, type = "link")
    enResiduals[,i] <- ylasso - enPredict[,i]
  }
  # Mean coefficient estimates
  coefMeans <- data.frame(param = row.names(coef(bsElasticNet)),
                          estimate = signif(rowMeans(coefResults), digits = 3),
                          median = apply(coefResults, 1, median),
                          sd = signif(apply(coefResults, 1, sd), digits = 3),
                          ciu = signif(apply(coefResults, 1, function(x) quantile(x, 0.025)), digits = 3),
                          cil = signif(apply(coefResults, 1, function(x) quantile(x, 0.975)), digits = 3))
  write.csv(coefMeans, file = paste(outdir,"elastic_net_mean_model_coef.csv", sep=""), row.names = FALSE)
  # Mean residuals
  residMeans <- data.frame(mean = signif(rowMeans(enResiduals), digits = 3),
                           sd = signif(apply(enResiduals, 1, sd), digits = 3))
  write.csv(residMeans, file = paste(outdir,"elastic_net_mean_residuals.csv", sep=""), row.names = FALSE)
  enlist <- list(ylasso = ylasso,
                 xlasso = xlasso,
                 coefResults = coefResults,
                 enPredict = enPredict,
                 enResiduals = enResiduals,
                 coefMeans = coefMeans,
                 residMeans = residMeans)
  return(enlist)
}
