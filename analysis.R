# laying out outcome variables and different econometric analyses

library(sandwich)
library(plm)
library(lmtest)
library(clubSandwich)
library(DataCombine)


################ addressing the dropping of observations in periods t+1 when pixel deforested in period t
# binary_coeffdist_fcn(n, nobs, years)
binary_coeffdist_fcn(200, 10000, 2)
# function provides mean and variance of estimates
#col 1 and 2 drop
#col 3 and 4 keep
#col 1 and 3 simple DID
#col 2 and 4 two-way FE


################ aggregating with property level perturbations ########################


####creating deforrate1 outcome variable for both dataframes
year1 <- subset(countyprop_df, year== 1)
colnames(year1)[colnames(year1)=="defor"] <- "defor0"
year1 <- subset(year1, select = c(countyid, defor0))
st_geometry(year1) <- NULL
countyprop_df <-  merge(countyprop_df, year1, by = "countyid")
countyprop_df$deforrate1 <- countyprop_df$defor / countyprop_df$defor0

year1 <- subset(proppert_df, year== 1)
colnames(year1)[colnames(year1)=="defor"] <- "defor0"
year1 <- subset(year1, select = c(property, defor0))
st_geometry(year1) <- NULL
proppert_df <-  merge(proppert_df, year1, by = "property")

proppert_df$deforrate1 <- proppert_df$defor / proppert_df$defor0
#### 

# agg_coeffdist_fcn(outcome, n, nobs, years, psize, cellsize)
agg_coeffdist_fcn(deforrate1, 100, 10000, 1000, 25)
# function provides mean and variance of estimates for two way FE model
# col 1 is property level
# col 2 is county level




