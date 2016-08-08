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
  lassoCVmin <- cv.glmnet(y = ylasso, x = xlasso, family = "gaussian", alpha = x)$lambda.min
  return(c(x, lassoCVmin))
}


#### Function to run elastic net multiple times for a given response variable
#### Inputs are:
# data = data set; should be left as explVars
# y = response variable name
# times = number of CV runs
# alphaValues = vector of values of alpha to evaluate
#### Need to make a folder/directory with the name of the response variable before hand
#### Future version should make directory automatically

ElasticNetFunction <- function(y = "lnlambda", times = 2000, data = explVars, alphaValues = seq(0.8,1,by=0.01)){
  #data <- explVars
  # Vector of response variable
  ylasso <- data[,y]
  assign("ylasso", ylasso, envir = .GlobalEnv)
  # Make as a matrix
  xlasso <- as.matrix(data[,-which(names(data) == y)])
  assign("xlasso", xlasso, envir = .GlobalEnv)
  # Output directory
  outdir <- paste("output/", y, "/", sep = "")
  ######################################################################################
  #### Elastic net LASSO
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
  #### Get model coefficients for every value of lambda at best alpha
  bsElasticNet <- glmnet(y = ylasso, x = xlasso, family = "gaussian", alpha = alphaBest)
  # glmnet for best alpha and lambda
  #coef(bsElasticNet, s = lambdaBest)
  # Combine coefficients into a matrix and get mean coefficient for each term
  bscoefResults <- matrix(0, nrow = ncol(xlasso)+1, ncol = length(lambdas))
  for(i in 1:length(lambdas)){
    lambda.i <- lambdas[i]
    bscoefResults[,i] <- as.numeric(coef(bsElasticNet, s = lambda.i))
  }
  # Mean coefficient estimates
  bscoefMeans <- data.frame(param = row.names(coef(bsElasticNet)),
                            estimate = signif(rowMeans(bscoefResults), digits = 3),
                            sd = signif(apply(bscoefResults, 1, sd), digits = 3))
  write.csv(bscoefMeans, file = paste(outdir,"elastic_net_mean_model_coef.csv", sep=""), row.names = FALSE)
  enlist <- list(alphaBest = alphaBest,
                 lambdaBest = lambdaBest,
                 alphaSummary = alphaSummary[alphaSummary$median == min(alphaSummary$median),],
                 # Summary of best lambdas at best alpha value
                 lambdaSummary = summary(lambdas))
  print(enlist)
  return(enlist)
} 