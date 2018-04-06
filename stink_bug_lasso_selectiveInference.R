#### LASSO/Elastic Net Analysis for Dawn's landscape stink bug reproduction ms

# Preliminaries
rm(list = ls())
my.packages <- c("openxlsx", "tidyr", "ggplot2",  
                 "dplyr", "glmnet", "lars", "covTest",
                 "selectiveInference", "glmnetUtils")
lapply(my.packages, require, character.only = TRUE)

source("selectiveInference_functions.R")
# source("functions_lars.R")

#### Brown stink bug reproduction data set
#### Explanatory variables can be grouped into three general categories:
# Landscape features:
# GV = proportion of woodland and pasture
# PA = perimeter to area ratio of cropland
# pmaize = proportion of maize in the landccape
# pcot = proportion of cotton in the landscape
# ppea = proportion of peanut in the landscape
# psoy = proportion of soybean in the landscape
# pall = proportion of all crop fields in the landscape
# ndist.... variables = no. of fields of a given crop within a given distance
# 
# Local woodland features
# pedge = proportion of non-crop hosts in the local area
# 
# Natural enemies:
# ne = no. of natural enemies (all species combined)
# geo = no. geocoris spp.
# ant = no. of fire ants

# Import data set
bsdata <- read.csv("data/Brown_lambda_lasso_data.csv", header = TRUE)

# Make region a 0/1 binary variable
bsdata$region <- bsdata$region - 1

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

# log transform mant and mgeo as they're highly skewed
bsdata$lngeo <- log(bsdata$geo + 1)
bsdata$lnant <- log(bsdata$ant + 1)
bsdata$lnne <- log(bsdata$ne + 1)

#bsdata$lnlambda <- factor2numeric(bsdata$lnlambda)
bsdata$ndist100m_all <- factor2numeric(bsdata$ndist100m_all)
bsdata$ndist500all <- factor2numeric(bsdata$ndist500all)

summary(bsdata)

##########################################################################################
#### Variable selection using LASSO
# For predictor variables, I'm including year and crop as binary dummy variables 
# Preparing data
datalasso <- bsdata %>% dplyr::select(., lnlambda, gv, pa, lngeo, lnant, lnne,
                    pmaize, pcot, ppea, psoy, pall,
                    starts_with("ndist"), # Selects for 15 variables -- each crop and all crops at 100, 500, and 1000 m distances
                    cropmaize, cropcotton, croppeanut, cropsoybean,
                    year2009, year2010, year2011) %>%
  scale(., center = TRUE, scale = TRUE) %>% 
  as.data.frame()
  # standardize all columns before running through glmnet, as recommended in selectiveInference manual
  # scale() converts datalasso into a matrix, converting back to data.frame for filter()

# Remove NAs
datalasso <- datalasso %>% dplyr::filter(., complete.cases(.)) 

# # Set of histograms for standardized covariates
# niceNames <- c("proportion green veining", "perimeter to area ratio", "mean predators", "mean geocoris",
#                "mean fire ants", "proportion maize", "proportion cotton", "proportion peanut", "proportion soybean",
#                "proportion total crops", "distance to maize", "distance to cotton", "distance to peanut", "distance to soybean",
#                "distance to all crops")
# for(i in 1:length(ccovars)){
#   covar.i <- ccovars[i]
#   fileName <- paste("figures/histogram_", covar.i, ".tif", sep="")
#   print(fileName)
#   tiff(file = fileName)
#     hist(explVars[,covar.i], xlab = niceNames[i], main = niceNames[i])
#   dev.off()
# }


######################################################################################
#### Elastic Net LASSO for stink bug lambda estimates
#### Cross-validation of alpha and lambda using the glmentUtils package

# cva.glment can take data in formula form, but submitting as matrix and vector to use in fixedLassoInf() too
ylasso <- as.numeric(datalasso$lnlambda)
xlasso <- datalasso %>% dplyr::select(., -lnlambda) %>% as.matrix()

alphaValues <- c(seq(0, 0.9, by=0.1), seq(0.92, 1, by=0.02))

sbcv <- cva.glmnet(y = ylasso, x = xlasso, alpha = alphaValues, standardize = FALSE)

minlossplot(sbcv, cv.type = "1se", type = "b") # How do I interpret this?

# Combine lambda.min and alpha for each model run
# sbcv$alpha = vector of alpha values
# sbcv$modlist = list of glmnet model objects
cvAlphaResults <- data.frame(alpha = sbcv$alpha,
                             lambda.min = sapply(1:length(sbcv$modlist), function(x) sbcv$modlist[[x]]$lambda.min, simplify = TRUE), 
                             lambda.1se = sapply(1:length(sbcv$modlist), function(x) sbcv$modlist[[x]]$lambda.1se, simplify = TRUE))

# plot alpha vs lambda.min
#tiff("output/selectiveInference/lnlambda/cross-validation_alpha-lambda_plot.tif")
  plot(x = cvAlphaResults$alpha, y = cvAlphaResults$lambda.1se, type = "b")
#dev.off()

# Alpha that produces best (lowest) lambda.1se
bestLambda <- min(cvAlphaResults$lambda.1se)
alphaBest <- cvAlphaResults[which(cvAlphaResults$lambda.1se == bestLambda), "alpha"]
alphaBest

# Run fixedLassoInf() and get clean results using getLassoInfResults()
sbLassoOutput <- getLassoInfResults(y = ylasso, x = xlasso, alpha = alphaBest, lambda = lambdaBest)
sbResults <- sbLassoOutput$lassoTestResults[order(sbLassoOutput$lassoTestResults$coef), c("varName", "coef", "pvalue")]
write.csv(sbResults, file = "output/selectiveInference/lnlambda/stink_bug_lasso_results.csv", row.names = FALSE)


#####################################################################################################################
#### selectiveInference using lars
#####################################################################################################################
# Fit lar model

sblar <- fitLAR(x = xlasso, y = ylasso)

write.csv(sblar$larResults, file = "output/selectiveInference/lnlambda/stink_bug_lasso_results.csv", row.names = FALSE)

# Combine LAR and Ordinary Least Squares coefficient estimates
OLSResults <- data.frame(varIndex = as.numeric(row.names(as.data.frame(sblar$larfit$bls))),
                         lscoef = sblar$larfit$bls)
larOLSResults <- sblar$larResults %>% dplyr::select(., -varNames) %>% full_join(., OLSResults, by = "varIndex")

# Get variable names
varNamesDF <- data.frame(varIndex = as.numeric(row.names(as.data.frame(varNames))),
                         varNames = varNames)
varNamesFullDF <- read.xlsx("data/lasso_full_variable_names.xlsx", sheet = "stink_bug_lambda")
varNamesDF <- full_join(varNamesDF, varNamesFullDF, by = "varIndex")
# Combine LAR results and variable names
larOLSResults <- full_join(varNamesDF, larOLSResults, by = "varIndex")
larOLSResults$coef[is.na(larOLSResults$coef)] <- 0
#larOLSResults <- larOLSResults %>% dplyr::arrange(., -abs(coef))
larOLSResults

# Make a pretty results table
resultsTable <- larOLSResults %>% mutate_at(., vars(coef, pvalue, cil, ciu, lscoef), funs(signif(., digits = 3))) %>% 
  mutate(., summary = ifelse(is.na(cil), coef, paste(coef, " [", cil, ", ", ciu, "]", sep = ""))) %>%
  dplyr::arrange(., -abs(coef)) %>%
  dplyr::select(., varNamesFull, summary, pvalue, lscoef)

write.csv(resultsTable, file = "output/selectiveInference/lnlambda/stink_bug_lasso_least-squares_results.csv", row.names = FALSE)



#####################################################################################################################
#### Examples from the web using fixedLassoInf()
#set.seed(43)
n = 50
p = 10
sigma = 1
x = matrix(rnorm(n*p),n,p)
x=scale(x,TRUE,TRUE)
beta = c(3,2,rep(0,p-2))
y = x%*%beta + sigma*rnorm(n)
# first run glmnet
gfit = glmnet::glmnet(x,y,standardize=FALSE)
# extract coef for a given lambda; note the 1/n factor!
# (and we don't save the intercept term)
lambda = .8
beta = coef(gfit, s=lambda/n, exact=TRUE)[-1]
# compute fixed lambda p-values and selection intervals
out = fixedLassoInf(x,y,beta,lambda,sigma=sigma)
out

data(diabetes)
G <- glmnet::glmnet(diabetes$x, diabetes$y)
cvG <- glmnet::cv.glmnet(diabetes$x, diabetes$y)
plot(cvG)
log(cvG$lambda.1se)
log(cvG$lambda.min)
G$lambda[20]
cvG$lambda.1se/nrow(diabetes$x)
beta.hat <- coef(G, s = G$lambda[20], exact = TRUE)[-1]
