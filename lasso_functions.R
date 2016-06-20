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
