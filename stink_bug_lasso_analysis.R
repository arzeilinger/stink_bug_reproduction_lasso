#### LASSO/Elastic Net Analysis for Dawn's landscape stink bug reproduction ms

# Preliminaries
rm(list = ls())
my.packages <- c("lattice", "tidyr", "ggplot2",  
                 "dplyr", "glmnet")
lapply(my.packages, require, character.only = TRUE)

source("lasso_functions.R")

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
bsdata <- read.csv("data/Brown_lambda_multi-regression.csv", header = TRUE)

# make a new variable for total proportion of crop area
bsdata$pcrop <- with(bsdata, pmaize + pcot + ppea + psoy)

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

str(bsdata)

##########################################################################################
#### Variable selection using LASSO
# For predictor variables, I'm including year and crop as binary dummy variables 
# Preparing data
datalasso <- bsdata[,c("lnlambda", "region", "gv", "pa", "ne", "mgeo", "mant", 
                    "pmaize", "pcot", "ppea", "psoy", "pcrop",
                    "mdistcorn", "mdistcot", "mdistpea", "mdistsoy", "mdistall",
                    "cropmaize", "cropcotton", "croppeanut", "cropsoybean",
                    "year2009", "year2010", "year2011")]
# Need only numeric variables for cv.glmnet
# datalasso$region <- as.numeric(levels(datalasso$region))[datalasso$region]
# datalasso$crop <- as.numeric(levels(datalasso$crop))[datalasso$crop]
# Remove NAs
badRows.i <- unlist(sapply(1:ncol(datalasso), function(x) which(is.na(datalasso[,x])), simplify = TRUE))
explVars <- datalasso[-badRows.i,]

# log transform mant and mgeo as they're highly skewed
explVars$lnmgeo <- log(explVars$mgeo + 1)
explVars$lnmant <- log(explVars$mant + 1)

# Make matrix of standardized covariates
# Standardize continuous covariates
ccovars <- c("lnlambda", "gv", "pa", "ne", "mgeo", "mant", 
             "pmaize", "pcot", "ppea", "psoy", "pcrop",
             "mdistcorn", "mdistcot", "mdistpea", "mdistsoy", "mdistall",
             "lnmgeo", "lnmant")
ccovars.i <- as.numeric(sapply(ccovars, function(x) which(names(explVars) == x), simplify = TRUE))
# Replace covariates with standardized form
for(i in ccovars.i){
  var.i <- names(explVars)[i]
  stdvar.i <- standardize(explVars[,var.i])
  explVars[,var.i] <- stdvar.i
}


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
print("Stink Bug Lambda")
sb <- ElasticNetFunction(y = "lnlambda", times = 2000,
                         data = explVars[,-which(names(explVars) == c("lnmgeo", "lnmant"))],
                         alphaValues = seq(0.8,1,by=0.01))
rm(xlasso)
rm(ylasso)


######################################################################################
#### LASSO for Geocoris and fire ant densities
#### remove natural enemy densities as is likely correlated with these variables
print("Geocoris density")

geoVars <- select(explVars, -mgeo, -ne, -lnmant)

#### Cross-validation for geocoris
geocorisCV <- encv(y = "lnmgeo", times = 2000, 
                 data = geoVars, 
                 alphaValues = seq(0.90,1,by=0.01))
rm(xlasso)
rm(ylasso)

saveRDS(geocorisCV, file = "output/lnmgeo/geocoris_cross-validation_output.rds")


#### Elasitc Net results and residuals for geocoris
geocorisEN <- ElasticNetFunction(y = "lnmgeo", data = geoVars, 
                                 alphaBest = geocorisCV$alphaBest, 
                                 lambdas = geocorisCV$lambdas)
saveRDS(geocorisEN, file = "output/lnmgeo/geocoris_elastic_net_output.rds")

explVars$geoResiduals <- geocorisEN$residMeans$mean
explVars$geoResidualsSD <- geocorisEN$residMeans$sd
write.csv(explVars, file = "output/lnmgeo/lambda_data_with_residuals.csv", row.names = FALSE)

# Checking the residuals
geoResiduals <- geocorisEN$residMeans$mean
geoPredict <- rowMeans(geocorisEN$enPredict)
plot(geoPredict, geoResiduals)

plot(ylasso, geoPredict, xlim = c(-2,5), ylim = c(-2,5))
abline(a = 0, b =1)


######################################################################################
#### Elastic Net LASSO for stink bug lambda estimates
print("Fire ant density")

antVars <- select(explVars, -lnmgeo, -ne, -mant)

rm(xlasso)
rm(ylasso)

#### Cross-validation for fire ant density
antCV <- encv(y = "lnmant", times = 2000, 
              data = antVars, 
              alphaValues = seq(0.8,1,by=0.01))

saveRDS(antCV, file = "output/lnmant/ant_cross-validation_output.rds")


#### Elasitc Net results and residuals for fire ant density
antEN <- ElasticNetFunction(y = "lnmant", data = antVars, 
                                 alphaBest = antCV$alphaBest, 
                                 lambdas = antCV$lambdas)
saveRDS(antEN, file = "output/lnmant/ant_elastic_net_output.rds")


#########################################################################################





