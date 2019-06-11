# laying out outcome variables and different econometric analyses

library(sandwich)
library(plm)
library(lmtest)
library(clubSandwich)
library(DataCombine)


################ addressing the dropping of observations in periods t+1 when pixel deforested in period t
# binary_coeffdist_fcn(n, nobs, years)
# function provides mean and variance of estimates
#col 1 and 2 drop
#col 3 and 4 keep
#col 1 and 3 simple DID
#col 2 and 4 two-way FE
binary_coeffdist_fcn(200, 10000, 2)
# find that two-way fe dropping identifies ATT+Diff
# DID seems to identify ATT, but slightly biased
# keeping obs biases estimates positively
# notice: keeping obs, twoway and DID are identical




################ aggregating with property level perturbations ########################

prop_deforsim(20000, 4, 5000,40)


################## creating the various outcome variables ############################

######################### defor_t/defor_0 ###########################
### creating deforrate1 outcome variable for both dataframes
#county level
year1 <- subset(countypert_df, year== 1)
colnames(year1)[colnames(year1)=="defor"] <- "defor0"
year1 <- subset(year1, select = c(countyid, defor0))
st_geometry(year1) <- NULL
countypert_df <-  merge(countypert_df, year1, by = "countyid")
countypert_df$ihsdeforrate1 <- sinh(countypert_df$defor / countypert_df$defor0)

countypert_df$deforrate2 <- countypert_df$defor / countypert_df$defor0


#property level
year1 <- subset(proppert_df, year== 1)
colnames(year1)[colnames(year1)=="defor"] <- "defor0"
year1 <- subset(year1, select = c(property, defor0))
st_geometry(year1) <- NULL
proppert_df <-  merge(proppert_df, year1, by = "property")

proppert_df$ihsdeforrate1 <- sinh(proppert_df$defor / proppert_df$defor0)

proppert_df$deforrate2 <- proppert_df$defor - proppert_df$defor0


######################### defor_t/defor_t-1 ###########################
### creating deforrate2 outcome variable for both dataframes

#county level
#order data for slide (lag)
countypert_df <- countypert_df[order(countypert_df$countyid, countypert_df$year),]

countypert_df <- slide(countypert_df, Var = "defor", GroupVar = "countyid", NewVar = "deforlag",
                   slideBy = -1)

countypert_df$deforrate2 <- countypert_df$defor / countypert_df$deforlag

#property level
proppert_df <- proppert_df[order(proppert_df$property, proppert_df$year),]

proppert_df <- slide(proppert_df, Var = "defor", GroupVar = "property", NewVar = "deforlag",
                       slideBy = -1)

proppert_df$deforrate2 <- proppert_df$defor / proppert_df$deforlag


######################### Busch et al (2014): defor_t-1-(defor_t/defor_t-1) ########################### 
### creating deforrate3 outcome variable for both dataframes

proppert_df$deforrate3 <- proppert_df$deforlag - (proppert_df$defor / proppert_df$deforlag)


countypert_df$deforrate3 <- countypert_df$deforlag- (countypert_df$defor / countypert_df$deforlag)


######################### Inverse Hyperbolic Lags ########################### 
### creating deforrate4 outcome variable for both dataframes




# agg_coeffdist_fcn(n, nobs, years, psize, cellsize)
agg_coeffdist_fcn(5, 10000, 1000, 25)
# function provides mean and variance of estimates for two way FE model
# col 1 is property level
# col 2 is county level




