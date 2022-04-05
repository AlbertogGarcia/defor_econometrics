library(tidyverse)
library(tictoc)
library(here)
library(DeclareDesign)
source(here::here('unbiased_dgp', 'full_landscape.R'))
source(here::here('unbiased_dgp', 'survival_did.R'))
source(here::here('unbiased_dgp', 'all_specifications.R'))

library(survival)
library(ggplot2)
library(dplyr)
library(ggfortify)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Parameterization (largely taken from analysis, but longer time periods) ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

years = 6
nobs = 100^2
n = 50

cellsize = 10
ppoints = 75
cpoints = 40

std_v = 0.5
std_a = 0
std_p = 0

# here are the landscape characteristics in this parameterization
# note that the bias for the TWFE model will be equal to the pre-treatment difference in deforestation rtes, which is 0.03
base_0 = .02
base_1 = .05
trend = -.005
ATT = -.01


###############################################################################################################
######## Baseline showing aggregation resolves pixel fixed effects issue
###############################################################################################################

std_avp = (std_a^2+std_v^2+std_p^2)^.5
b0 = qnorm(base_0, mean = 0, sd = std_avp)
b1 = qnorm(base_1, mean = 0, sd = std_avp) - b0
b2_0 = qnorm(trend + base_0, mean = 0, sd = std_avp) - b0
b2_1 = qnorm(trend + base_1, mean = 0, sd = std_avp) - b0 - b1
b3 = qnorm( pnorm(b0+b1+b2_1, mean = 0, sd = std_avp) + ATT , mean = 0, sd = std_avp) - (b0 + b1 + b2_1)

set.seed(0930)
# summary function that estimates all of the different specifications
aggregation <- all_specifications(n, nobs, years, b0, b1, b2_0, b2_1, b3, std_a, std_v, std_p, cellsize, ppoints, cpoints)

summary_long <- aggregation$summary_long

library(rio)
export(summary_long, "unbiased_dgp/results/summary_long.rds")

###############################################################################################################
######## Introduction pixel level unobservables, which impact non-random selection
###############################################################################################################
std_a = 0.1
# we'll need to recompute the parameters if we change the value of sigma_p

# we'll need to compute the parameters 
std_avp = (std_a^2+std_v^2+std_p^2)^.5
b0 = qnorm(base_0, mean = 0, sd = std_avp)
b1 = qnorm(base_1, mean = 0, sd = std_avp) - b0
b2_0 = qnorm(trend + base_0, mean = 0, sd = std_avp) - b0
b2_1 = qnorm(trend + base_1, mean = 0, sd = std_avp) - b0 - b1
b3 = qnorm( pnorm(b0+b1+b2_1, mean = 0, sd = std_avp) + ATT , mean = 0, sd = std_avp) - (b0 + b1 + b2_1)

#ATT = pnorm(b0+b1+b2_1+b3, 0, (std_a^2+std_v^2 )^.5) - pnorm(b0+b1+b2_1, 0, (std_a^2+std_v^2 )^.5)

# summary function that estimates all of the different specifications
aggregation_0 <- all_specifications(n, nobs, years, b0, b1, b2_0, b2_1, b3, std_a, std_v, std_p, cellsize, ppoints, cpoints)

summary_long_0 <- aggregation_0$summary_long 

export(summary_long_0, "unbiased_dgp/results/summary_selection.rds")
###############################################################################################################
######## Adding in property level disturbances
###############################################################################################################

#### 0.1
std_p = 0.1

# we'll need to compute the parameters 
std_avp = (std_a^2+std_v^2+std_p^2)^.5
b0 = qnorm(base_0, mean = 0, sd = std_avp)
b1 = qnorm(base_1, mean = 0, sd = std_avp) - b0
b2_0 = qnorm(trend + base_0, mean = 0, sd = std_avp) - b0
b2_1 = qnorm(trend + base_1, mean = 0, sd = std_avp) - b0 - b1
b3 = qnorm( pnorm(b0+b1+b2_1, mean = 0, sd = std_avp) + ATT , mean = 0, sd = std_avp) - (b0 + b1 + b2_1)

aggregation_1 <- all_specifications(n, nobs, years, b0, b1, b2_0, b2_1, b3, std_a, std_v, std_p, cellsize, ppoints, cpoints)

summary_long_1 <- aggregation_1$summary_long

##### 0.25
std_p = 0.25

# we'll need to compute the parameters 
std_avp = (std_a^2+std_v^2+std_p^2)^.5
b0 = qnorm(base_0, mean = 0, sd = std_avp)
b1 = qnorm(base_1, mean = 0, sd = std_avp) - b0
b2_0 = qnorm(trend + base_0, mean = 0, sd = std_avp) - b0
b2_1 = qnorm(trend + base_1, mean = 0, sd = std_avp) - b0 - b1
b3 = qnorm( pnorm(b0+b1+b2_1, mean = 0, sd = std_avp) + ATT , mean = 0, sd = std_avp) - (b0 + b1 + b2_1)

#ATT = pnorm(b0+b1+b2_1+b3, 0, (std_a^2+std_v^2 )^.5) - pnorm(b0+b1+b2_1, 0, (std_a^2+std_v^2 )^.5)

aggregation_2 <- all_specifications(n, nobs, years, b0, b1, b2_0, b2_1, b3, std_a, std_v, std_p, cellsize, ppoints, cpoints)

summary_long_2 <- aggregation_2$summary_long


#### 0.5
std_p = 0.5

std_avp = (std_a^2+std_v^2+std_p^2)^.5
b0 = qnorm(base_0, mean = 0, sd = std_avp)
b1 = qnorm(base_1, mean = 0, sd = std_avp) - b0
b2_0 = qnorm(trend + base_0, mean = 0, sd = std_avp) - b0
b2_1 = qnorm(trend + base_1, mean = 0, sd = std_avp) - b0 - b1
b3 = qnorm( pnorm(b0+b1+b2_1, mean = 0, sd = std_avp) + ATT , mean = 0, sd = std_avp) - (b0 + b1 + b2_1)

aggregation_3 <- all_specifications(n, nobs, years, b0, b1, b2_0, b2_1, b3, std_a, std_v, std_p, cellsize, ppoints, cpoints)

summary_long_3 <- aggregation_3$summary_long

summary_full <- rbind(summary_long_0, summary_long_1, summary_long_2, summary_long_3)

export(summary_full, "unbiased_dgp/results/summary_full.rds")


###############################################################################################################
######## alternative parameterization
###############################################################################################################

base_0 = .05
base_1 = .02
trend = -.005
ATT = -.01

# we'll need to compute the parameters 
std_avp = (std_a^2+std_v^2+std_p^2)^.5
b0 = qnorm(base_0, mean = 0, sd = std_avp)
b1 = qnorm(base_1, mean = 0, sd = std_avp) - b0
b2_0 = qnorm(trend + base_0, mean = 0, sd = std_avp) - b0
b2_1 = qnorm(trend + base_1, mean = 0, sd = std_avp) - b0 - b1
b3 = qnorm( pnorm(b0+b1+b2_1, mean = 0, sd = std_avp) + ATT , mean = 0, sd = std_avp) - (b0 + b1 + b2_1)

set.seed(0930)


aggregation_alt <- all_specifications(n, nobs, years, b0, b1, b2_0, b2_1, b3, std_a, std_v, std_p, cellsize, ppoints, cpoints)

summary_long_alt <- aggregation_alt$summary_long
export(summary_long_alt, "unbiased_dgp/results/summary_long_alt.rds")

base_0 = .02
base_1 = .05
trend = -.005
ATT = .01

# we'll need to compute the parameters 
std_avp = (std_a^2+std_v^2+std_p^2)^.5
b0 = qnorm(base_0, mean = 0, sd = std_avp)
b1 = qnorm(base_1, mean = 0, sd = std_avp) - b0
b2_0 = qnorm(trend + base_0, mean = 0, sd = std_avp) - b0
b2_1 = qnorm(trend + base_1, mean = 0, sd = std_avp) - b0 - b1
b3 = qnorm( pnorm(b0+b1+b2_1, mean = 0, sd = std_avp) + ATT , mean = 0, sd = std_avp) - (b0 + b1 + b2_1)

aggregation_alt2 <- all_specifications(n, nobs, years, b0, b1, b2_0, b2_1, b3, std_a, std_v, std_p, cellsize, ppoints, cpoints)

summary_long_alt2 <- aggregation_alt2$summary_long
export(summary_long_alt2, "unbiased_dgp/results/summary_long_alt2.rds")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# survival <- survival_did(n, nobs, years, b0, b1, b2_0, b2_1, b3, std_a, std_v, std_p, cellsize, ppoints, cpoints)
# 
# survival_long <- survival$summary_long

# library(rio)
# export(survival_long, "survival_long.rds")

