#### LASSO/Elastic Net Analysis for Dawn's landscape stink bug reproduction ms

# Preliminaries
rm(list = ls())
my.packages <- c("lattice", "tidyr", "ggplot2",  
                 "dplyr", "glmnet", "lars", "covTest")
lapply(my.packages, require, character.only = TRUE)

source("lasso_functions.R")
source("functions_lars.R")

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
bsdata <- read.csv("data/Brown_lambda_lasso_data.csv", header = TRUE)

# # For combining data sets when Dawn's updates ndist data
# bsdata <- bsdata %>% dplyr::select(., -starts_with("ndist")) %>% unite(., col = uniqueField, year:field, remove = FALSE)
# n_distinct(bsdata$uniqueField)
# ndistdata <- read.csv("data/Brown_lambda_lasso_data_ndist.csv")
# ndistdata <- ndistdata %>% unite(., col = uniqueField, year:field, remove = TRUE)
# bsdata2 <- full_join(bsdata, ndistdata, by = "uniqueField")
# bsdata <- bsdata2
# summary(bsdata)
# write.csv(bsdata, file = "data/Brown_lambda_lasso_data.csv", row.names = FALSE)

# # make a new variable for total proportion of crop area
# bsdata$pcrop <- with(bsdata, pmaize + pcot + ppea + psoy)

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

summary(bsdata)

##########################################################################################
#### Variable selection using LASSO
# For predictor variables, I'm including year and crop as binary dummy variables 
# Preparing data
datalasso <- bsdata %>% dplyr::select(., lnlambda, gv, pa, lngeo, lnant, lnne,
                    pmaize, pcot, ppea, psoy, pall,
                    starts_with("ndist"), # Selects for 15 variables -- each crop and all crops at 100, 500, and 1000 m distances
                    cropmaize, cropcotton, croppeanut, cropsoybean,
                    year2009, year2010, year2011)
# Need only numeric variables for cv.glmnet
# datalasso$region <- as.numeric(levels(datalasso$region))[datalasso$region]
# datalasso$crop <- as.numeric(levels(datalasso$crop))[datalasso$crop]
# Remove NAs
goodRows <- datalasso %>% complete.cases() %>% which()
explVars <- datalasso[goodRows,]

# Don't need to standardize before hand, just use glmnet's standardize option
# # Make matrix of standardized covariates
# # Standardize continuous covariates
# ccovars <- c("lnlambda", "gv", "pa", "ne", "mgeo", "mant", 
#              "pmaize", "pcot", "ppea", "psoy", "pcrop",
#              "mdistcorn", "mdistcot", "mdistpea", "mdistsoy", "mdistall",
#              "lnmgeo", "lnmant")
# ccovars.i <- as.numeric(sapply(ccovars, function(x) which(names(explVars) == x), simplify = TRUE))
# # Replace covariates with standardized form
# for(i in ccovars.i){
#   var.i <- names(explVars)[i]
#   stdvar.i <- standardize(explVars[,var.i])
#   explVars[,var.i] <- stdvar.i
# }


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

explVars$ndist100m_all <- factor2numeric(explVars$ndist100m_all)
explVars$ndist500all <- factor2numeric(explVars$ndist500all)
str(explVars)

######################################################################################
#### Elastic Net LASSO for stink bug lambda estimates
print("Stink Bug Lambda")

ti <- Sys.time()
sbcv <- encv(y = "lnlambda", times = 2000,
             data = explVars,
             alphaValues = seq(0.9, 1, by = 0.02))
saveRDS(sbcv, file = "output/lnlambda/stink_bug_lnlambda_cross-validation_output.rds")
tf <- Sys.time()

# Time took for cross-validation 
tf-ti


sben <- ElasticNetFunction(y = "lnlambda", data = explVars, 
                           alphaBest = sbcv$alphaBest, 
                           lambdas = sbcv$lambdas)
saveRDS(sben, file = "output/lnlambda/stink_bug_lnlambda_elastic_net_output.rds")

sben$coefMeans

sbcv <- readRDS("output/lnlambda/stink_bug_lnlambda_cross-validation_output.rds")
sben <- readRDS("output/lnlambda/stink_bug_lnlambda_elastic_net_output.rds")

#### LASSO significance test with lars and covTest

xsb <- sben$xlasso
ysb <- sben$ylasso

sbTest <- larsLASSOFunction(y = ysb, x = xsb)

# Combine all results, from glmnet and lars
sbResults <- full_join(sben$coefMeans, sbTest, by = "param")
write.csv(sbResults, file = "output/stink_bug_lambda_lasso_results.csv", row.names = FALSE)

sbResults
