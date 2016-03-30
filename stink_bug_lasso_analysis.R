#### LASSO/Elastic Net Analysis for Dawn's landscape stink bug reproduction ms

# Preliminaries
rm(list = ls())
my.packages <- c("lattice", "tidyr", "ggplot2",  
                 "dplyr", "glmnet", "piecewiseSEM")
lapply(my.packages, require, character.only = TRUE)
# Note: using the sem.model.fits() function from the piecewiseSEM package to get R-squared for
# LMMs, which is based on the publication:
# Nakagawa, S., and H. Schielzeth. 2013. A general and simple method for obtaining R2 from generalized 
# linear mixed-effects models. Methods in Ecology and Evolution 4(2): 133-142. 
# DOI: 10.1111/j.2041-210x.2012.00261.x

setwd("C:/Users/Adam/Documents/Dawns stink bug reproduction ms/Multiple regression")

#### Brown stink bug reproduction data set
#### Explanatory variables can be grouped into three general categories:
# Landscape features:
# GV = proportion of woodland and pasture
# PA = perimeter to area ratio of cropland
# pmaize = proportion of maize in the landccape
# pcot = proportion of cotton in the landscape
# ppea = proportion of peanut in the landscape
# psoy = proportion of soybean in the landscape
# mdistall = mean distance to all crops = 1000 m from the sampled field
# mdistmaize = mean distance to maize fields = 1000 m from the sampled field
# mdistcot = means distance to cotton fields = 1000 m from the sampled field
# mdistpea = mean distance to peanut fields = 1000 m from the sampled field
# 
# Local woodland features
# pedge = proportion of non-crop hosts in the local area
# 
# Natural enemies:
# NE = mean no. of natural enemies (all species combined)
# mgeo = mean no. geocoris spp.
# mant = mean no. of fire ants

# Import data set
bsdata <- read.csv("Brown_lambda_multi-regression.csv", header = TRUE)

# make a new variable for total proportion of crop area
bsdata$pcrop <- with(bsdata, pmaize + pcot + ppea + psoy)

# Crop factor codes
bsdata$crop <- factor(ifelse(bsdata$crop == 1, "maize", 
                      ifelse(bsdata$crop == 2, "cotton",
                             ifelse(bsdata$crop == 3, "peanut", "soybean"))))
# Change factor order to maize, cotton, peanut, soybean.
# Maize will be base level when using treatment contrasts
bsdata$crop <- factor(bsdata$crop, levels(bsdata$crop)[c(2,1,3,4)]) 
levels(bsdata$crop)
# Create binary dummy variables for each crop factor level
for(i in levels(bsdata$crop)){
  crop.i <- i
  cropname.i <- paste("crop", crop.i, sep = "")
  bsdata[,cropname.i] <- ifelse(bsdata$crop == crop.i, 1, 0)
}

# Make year a factor and create binary dummy variables for each year
bsdata$year <- factor(bsdata$year)
levels(bsdata$year)
for(i in levels(bsdata$year)){
  year.i <- i
  yearname.i <- paste("year", year.i, sep = "")
  bsdata[,yearname.i] <- ifelse(bsdata$year == year.i, 1, 0)
}

str(bsdata)

##########################################################################################
#### Variable selection using LASSO
# For predictor variables, I'm including year and crop as binary dummy variables 
# Preparing data
datalasso <- bsdata[,c("lnlambda", "region", "gv", "pa", "ne", "mgeo", "mant", 
                    "pmaize", "pcot", "ppea", "psoy", "pcrop",
                    "mdistcorn", "mdistcot", "mdistpea", "mdistsoy", "mdistall",
                    "pedge", "cropmaize", "cropcotton", "croppeanut", "cropsoybean",
                    "year2009", "year2010", "year2011")]
# Need only numeric variables for cv.glmnet
# datalasso$region <- as.numeric(levels(datalasso$region))[datalasso$region]
# datalasso$crop <- as.numeric(levels(datalasso$crop))[datalasso$crop]
# Remove NAs
badRows.i <- unlist(sapply(1:ncol(datalasso), function(x) which(is.na(datalasso[,x])), simplify = TRUE))
explVars <- datalasso[-badRows.i,]
# Make matrix of explanatory variables
xlasso <- as.matrix(explVars[,-1])
ylasso <- explVars$lnlambda


######################################################################################
#### Elastic net LASSO
#### When alpha = 1, glmnet runs a lasso penalty; when 0 < alpha < 1, it uses the elastic net lasso penalty
#### To cross-validate alpha, run cv.glmnet over a range of alpha
# Function to loop through alpha values
lassoAlphaCV <- function(x){
  lassoCVmin <- cv.glmnet(y = ylasso, x = xlasso, family = "gaussian", alpha = x)$lambda.min
  return(c(x, lassoCVmin))
}
times <- 2000
alphaVec <- rep(seq(0.9,1,by=0.01),times) # sequence of alpha values
alphaCVResults <- as.data.frame(t(sapply(alphaVec, lassoAlphaCV, simplify = TRUE))) # run through alpha values
names(alphaCVResults) <- c("alpha", "lambda.min")
plot(x = alphaCVResults$alpha, y = jitter(alphaCVResults$lambda.min))
## summary of alpha values
alphaSummary <- alphaCVResults %>% group_by(alpha) %>% summarise(mean = mean(lambda.min), 
                                                                 median = median(lambda.min),
                                                                 sd = sd(lambda.min))
pdf("alpha-lambda cross-validation plot median.pdf")
  plot(x = alphaSummary$alpha, y = alphaSummary$median, type = "b")
dev.off()
alphaSummary[alphaSummary$median == min(alphaSummary$median),]
write.csv(as.data.frame(alphaSummary), file = "alpha-lambda cross-validation.csv", row.names = FALSE)
#### Extract lambda.min values from all runs at best alpha
alphaBest <- as.numeric(alphaSummary[alphaSummary$median == min(alphaSummary$median),"alpha"])
lambdaBest <- as.numeric(min(alphaSummary$median))
lambdas <- alphaCVResults[alphaCVResults$alpha == alphaBest,"lambda.min"]

#### Get model coefficients for every value of lambda at best alpha
bsElasticNet <- glmnet(y = ylasso, x = xlasso, family = "gaussian", alpha = alphaBest)
# glmnet for best alpha and lambda
coef(bsElasticNet, s = lambdaBest)
# Combine coefficients into a matrix and get mean coefficient for each term
bscoefResults <- matrix(0, nrow = ncol(xlasso)+1, ncol = length(lambdas))
for(i in 1:length(lambdas)){
  lambda.i <- lambdas[i]
  bscoefResults[,i] <- as.numeric(coef(bsElasticNet, s = lambda.i))
}
# Mean coefficient estimates
bscoefMeans <- data.frame(param = row.names(coef(bsElasticNet)),
                          estimate = rowMeans(bscoefResults),
                          sd = apply(bscoefResults, 1, sd))
write.csv(bscoefMeans, file = "elastic_net_mean_model_coef.csv", row.names = FALSE)

# Which terms are estimated at zero (dropped from the model) at least once
is.zero <- sapply(1:nrow(bscoefResults), function(x) any(bscoefResults[x,] == 0), simplify = TRUE) %>%
  which() %>% row.names(coef(bsElasticNet))[.]
selectedParams <- bscoefMeans[bscoefMeans$estimate != 0, "param"]
sometimes.zero <- selectedParams[which(selectedParams %in% is.zero)]


#### Ridge regression
cvRidge <- cv.glmnet(y = ylasso, x = xlasso, family = "gaussian", alpha = 0)
lambdaRidge <- cvRidge$lambda.min
bsRidge <- glmnet(y = ylasso, x = xlasso, family = "gaussian", alpha = 0)
coef(bsRidge, s = lambdaRidge)


