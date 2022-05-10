library(tidyverse)
library(ggplot2)
library(here)
source(here::here('unbiased_dgp', 'propFE.R'))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Parameterization 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

years = 6
nobs = 100^2
n = 100

cellsize = 20
ppoints = 100
cpoints = 25

avg_parea  = nobs/ppoints
avg_carea = nobs/cpoints

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
aggregation <- propFE(n, nobs, years, b0, b1, b2_0, b2_1, b3, std_a, std_v, std_p, cellsize_list = as.list(cellsize), ppoints, cpoints, proptreatassign = FALSE)

summary_long <- aggregation$summary_long

aggregation_prop <- propFE(n, nobs, years, b0, b1, b2_0, b2_1, b3, std_a, std_v, std_p, cellsize_list = as.list(cellsize), ppoints, cpoints, proptreatassign = TRUE)

summary_long_prop <- aggregation_prop$summary_long


###############################################################################################################
######## Introduction pixel level unobservables, which impact non-random selection
###############################################################################################################
std_a = 0.1
std_p = 0.5
# we'll need to recompute the parameters if we change the value of sigma_p

# we'll need to compute the parameters 
std_avp = (std_a^2+std_v^2+std_p^2)^.5
b0 = qnorm(base_0, mean = 0, sd = std_avp)
b1 = qnorm(base_1, mean = 0, sd = std_avp) - b0
b2_0 = qnorm(trend + base_0, mean = 0, sd = std_avp) - b0
b2_1 = qnorm(trend + base_1, mean = 0, sd = std_avp) - b0 - b1
b3 = qnorm( pnorm(b0+b1+b2_1, mean = 0, sd = std_avp) + ATT , mean = 0, sd = std_avp) - (b0 + b1 + b2_1)

# summary function that estimates all of the different specifications
aggregation_5 <- propFE(n, nobs, years, b0, b1, b2_0, b2_1, b3, std_a, std_v, std_p, cellsize_list = as.list(cellsize), ppoints, cpoints, proptreatassign = FALSE)

summary_long_5 <- aggregation_5$summary_long 

aggregation_prop5 <- propFE(n, nobs, years, b0, b1, b2_0, b2_1, b3, std_a, std_v, std_p, cellsize_list = as.list(cellsize), ppoints, cpoints, proptreatassign = TRUE)

summary_long_prop5 <- aggregation_prop5$summary_long

final_county <- rbind(summary_long, summary_long_5) %>%
  mutate_at(vars(bias, cover, std_p), as.numeric)%>%
  group_by(std_p, treatment.var) %>%
  summarize(Bias = mean(bias))%>%
  mutate(assignment = "county")

final_prop <- rbind(summary_long_prop, summary_long_prop5) %>%
  mutate_at(vars(bias, cover, std_p), as.numeric)%>%
  group_by(std_p, treatment.var) %>%
  summarize(Bias = mean(bias))%>%
  mutate(assignment = "property")

final <- rbind(final_county, final_prop)

library(rio)
export(final, "propFE_explore.rds")