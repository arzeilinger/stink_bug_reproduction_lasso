#########################################################################################
#### Analyses of natural enemy densities and distance from field edges
#########################################################################################

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
# Create binary dummy variables for each crop factor level
for(i in levels(nedata$crop)){
  crop.i <- i
  cropname.i <- paste("crop", crop.i, sep = "")
  nedata[,cropname.i] <- ifelse(nedata$crop == crop.i, 1, 0)
}

# Make year a factor and create binary dummy variables for each year
nedata$year <- factor(nedata$year)
levels(nedata$year)
for(i in levels(nedata$year)){
  year.i <- i
  yearname.i <- paste("year", year.i, sep = "")
  nedata[,yearname.i] <- ifelse(nedata$year == year.i, 1, 0)
}

str(nedata)

##########################################################################################
#### Variable selection using LASSO
# Preparing data
datalasso <- nedata %>% dplyr::select(., -date, -site, -field, -trans, -sample, -year, -crop)
# Need only numeric variables for cv.glmnet
# datalasso$region <- as.numeric(levels(datalasso$region))[datalasso$region]
# datalasso$crop <- as.numeric(levels(datalasso$crop))[datalasso$crop]
# Remove NAs
badRows.i <- unlist(sapply(1:ncol(datalasso), function(x) which(is.na(datalasso[,x])), simplify = TRUE))
explVars <- datalasso[-badRows.i,]

# log transform mant and mgeo as they're highly skewed
explVars$lngeo <- log(explVars$geo + 1)
explVars$lnant <- log(explVars$ant + 1)


#############################################################################################
#### Elastic Net for geocoris

geoVars <- dplyr::select(explVars, -geo, -ant, -lnant)

#### Cross-validation for geocoris
geocorisCV <- encv(y = "lngeo", times = 500, 
                   data = geoVars, 
                   alphaValues = seq(0,1,by=0.2))

saveRDS(geocorisCV, file = "output/lngeo/geocoris_cross-validation_output.rds")


#### Elasitc Net results and residuals for geocoris
geocorisEN <- ElasticNetFunction(y = "lngeo_edge", data = geoVars, 
                                 alphaBest = geocorisCV$alphaBest, 
                                 lambdas = geocorisCV$lambdas)
saveRDS(geocorisEN, file = "output/lngeo_edge/geocoris_edge_elastic_net_output.rds")

explVars$geoResiduals <- geocorisEN$residMeans$mean
explVars$geoResidualsSD <- geocorisEN$residMeans$sd
write.csv(explVars, file = "output/lngeo_edge/lambda_data_with_residuals.csv", row.names = FALSE)


#############################################################################################
#### Elastic Net for fire ant densities

antVars <- dplyr::select(explVars, -geo, -ant, -lngeo_edge)

#### Cross-validation for ant
antCV <- encv(y = "lnant_edge", times = 2000, 
              data = antVars, 
              alphaValues = seq(0.95,1,by=0.01))
rm(xlasso)
rm(ylasso)

saveRDS(antCV, file = "output/lnant_edge/ant_edge_cross-validation_output.rds")


#### Elasitc Net results and residuals for ant
antEN <- ElasticNetFunction(y = "lnant_edge", data = antVars, 
                            alphaBest = antCV$alphaBest, 
                            lambdas = antCV$lambdas)
saveRDS(antEN, file = "output/lnant_edge/ant_edge_elastic_net_output.rds")

explVars$antResiduals <- antEN$residMeans$mean
explVars$antResidualsSD <- antEN$residMeans$sd
write.csv(explVars, file = "output/lnant_edge/lambda_data_with_residuals.csv", row.names = FALSE)
