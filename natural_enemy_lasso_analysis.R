#########################################################################################
#### Analyses of natural enemy densities and distance from field edges
#########################################################################################

#### Cross-validation conducted on Vector Computer Cluster, in parallel; 
#### processing of crass-validation results done here. Still need to clean up data set here for coef estimation


# Preliminaries
rm(list = ls())
my.packages <- c("tidyr", "dplyr", "glmnet")
lapply(my.packages, require, character.only = TRUE)

source("lasso_functions.R")

nedata <- read.csv("data/NE_all_year.csv", header = TRUE)

# Make region a 0/1 binary variable
nedata$region <- nedata$region - 1

# Crop factor codes
nedata$crop <- factor(ifelse(nedata$crop == 1, "maize", 
                             ifelse(nedata$crop == 2, "cotton",
                                    ifelse(nedata$crop == 3, "peanut", "soybean"))))
# Change factor order to maize, cotton, peanut, soybean.
# Maize will be base level when using treatment contrasts
nedata$crop <- factor(nedata$crop, levels(nedata$crop)[c(2,1,3,4)]) 
levels(nedata$crop)

# Select only variables of interest for LASSO
datalasso <- nedata %>% dplyr::select(., year, crop, ant, geo, dedge, gv, pa, pmaize, pcot, ppea, psoy, pall, starts_with("ndist"))

# Create binary dummy variables for each crop factor level
for(i in levels(datalasso$crop)){
  crop.i <- i
  cropname.i <- paste("crop", crop.i, sep = "")
  datalasso[,cropname.i] <- ifelse(datalasso$crop == crop.i, 1, 0)
}

# Make year a factor and create binary dummy variables for each year
datalasso$year <- factor(datalasso$year)
levels(datalasso$year)
for(i in levels(datalasso$year)){
  year.i <- i
  yearname.i <- paste("year", year.i, sep = "")
  datalasso[,yearname.i] <- ifelse(datalasso$year == year.i, 1, 0)
}

# Remove year and crop variables now that they've been split into dummy binary variables
datalasso <- datalasso %>% dplyr::select(., -year, -crop)
str(datalasso)

# Convert factor variables to numeric
datalasso <- allFactor2numeric(datalasso)

# Remove NAs
badRows.i <- unlist(sapply(1:ncol(datalasso), function(x) which(is.na(datalasso[,x])), simplify = TRUE))
explVars <- datalasso[-badRows.i,]

# log transform mant and mgeo as they're highly skewed
explVars$lngeo <- log(explVars$geo + 1)
explVars$lnant <- log(explVars$ant + 1)


#############################################################################################
#### Elastic Net for geocoris

geoVars <- dplyr::select(explVars, -geo, -lnant)
summary(geoVars)

#### Cross-validation results from parallel run (with doParallel package) on cluster
geocorisCV <- readRDS("output/lngeo/geocoris_cross-validation_output.rds")
str(geocorisCV)
# extract alphaSummary from each list
alphaSummary <- matrix(0, ncol = 4, nrow = length(geocorisCV)) %>% as.data.frame()
for(i in 1:length(geocorisCV)){
  alphaSummary[i,] <- geocorisCV[[i]]$alphaSummary
}
names(alphaSummary) <- names(geocorisCV[[1]]$alphaSummary)
# determine best alpha value
alphaBest <- alphaSummary %>% dplyr::filter(median == min(median)) %>% dplyr::select(alpha)
alphaBestIndex <- which(alphaSummary$alpha == as.numeric(alphaBest))
# Best alpha = 1

#### Elasitc Net results and residuals for geocoris
lambdasBest <- geocorisCV[[alphaBestIndex]]$alphaCVResults %>% dplyr::filter(., alpha == alphaBest) %>% dplyr::select(., lambda.min) %>% as.matrix()
geocorisEN <- ElasticNetFunction(y = "lngeo", data = geoVars, 
                                 alphaBest = alphaBest, 
                                 lambdas = lambdasBest)
saveRDS(geocorisEN, file = "output/lngeo/geocoris_elastic_net_output.rds")


#############################################################################################
#### Elastic Net for fire ant densities

antVars <- dplyr::select(explVars, -ant, -lngeo)

#### Cross-validation results from parallel run (with doParallel package) on cluster
antCV <- readRDS("output/lnant/ant_cross-validation_output.rds")
str(antCV)
# extract alphaSummary from each list
alphaSummary <- matrix(0, ncol = 4, nrow = length(antCV)) %>% as.data.frame()
for(i in 1:length(antCV)){
  alphaSummary[i,] <- antCV[[i]]$alphaSummary
}
names(alphaSummary) <- names(antCV[[1]]$alphaSummary)
# determine best alpha value
alphaBest <- alphaSummary %>% dplyr::filter(median == min(median)) %>% dplyr::select(alpha)
alphaBestIndex <- which(alphaSummary$alpha == as.numeric(alphaBest))
# Best alpha = 1

#### Elasitc Net results and residuals for ant
lambdasBest <- antCV[[alphaBestIndex]]$alphaCVResults %>% dplyr::filter(., alpha == alphaBest) %>% dplyr::select(., lambda.min) %>% as.matrix()
antEN <- ElasticNetFunction(y = "lnant", data = antVars, 
                                 alphaBest = alphaBest, 
                                 lambdas = lambdasBest)
saveRDS(antEN, file = "output/lnant/ant_elastic_net_output.rds")

#### Elasitc Net results and residuals for ant
antEN <- ElasticNetFunction(y = "lnant", data = antVars, 
                            alphaBest = antCV$alphaBest, 
                            lambdas = antCV$lambdas)
saveRDS(antEN, file = "output/lnant/ant_elastic_net_output.rds")

