# laying out outcome variables and different econometric analyses

library(sandwich)
library(plm)
library(lmtest)
library(clubSandwich)
library(DataCombine)
library(plm)

#agg_coeffdist_fcn(n, nobs, years, psize, cellsize)
agg_coeffdist_fcn(200, 10000, 3, 1000, 25)




#binary outcomes without perturbations or counties/properties
# addressing the dropping of observations in periods t+1 when pixel deforested in period t
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
