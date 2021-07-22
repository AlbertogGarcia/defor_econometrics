# script to perform simulation runs and export csv files

source(here::here('unbiased_dgp', 'aggregate_complete.R'))
source(here::here('unbiased_dgp', 'TWFE_fcn.R'))

# we start with our base parameterization without property level perturbations
std_a = 0.1
std_v = 0.25
years = 2
nobs = 14400
n = 200

cellsize = 10
ppoints = 70
std_p = 0
cpoints = 40

# here are the landscape characteristics in this parameterization
# note that the bias for the TWFE model will be equal to the pre-treatment difference in deforestation rtes, which is 0.03
base_0 = .02
base_1 = .05
trend = -.005
ATT = -.01

# we'll need to compute the parameters 
std_avp = (std_a^2+std_v^2+std_p^2)^.5
b0 = qnorm(base_0, mean = 0, sd = std_avp)
b1 = qnorm(base_1, mean = 0, sd = std_avp) - b0
b2_0 = qnorm(trend + base_0, mean = 0, sd = std_avp) - b0
b2_1 = qnorm(trend + base_1, mean = 0, sd = std_avp) - b0 - b1
b3 = qnorm( pnorm(b0+b1+b2_1, mean = 0, sd = std_avp) + ATT , mean = 0, sd = std_avp) - (b0 + b1 + b2_1)

# runs simulation for all models
aggregation_0.03 <- aggregate_complete(n, nobs, years, b0, b1, b2_0, b2_1, b3, std_a, std_v, std_p, cellsize, ppoints, cpoints)
#set long dataframe with summary for each iteration
summary_long_0.03 <- aggregation_0.03$summary_long

library(rio)
export(summary_long_0.03, "specchart_long.rds")

TWFE_0.03 <- TWFE_fcn(n, nobs, years, b0, b1, b2_0, b2_1, b3, std_a, std_v)
TWFE_long_0.03 <- TWFE_0.03$summary_long



######################################################################################################################
# in order to explore how twfe bias changes with pre-treatment difference in deforestation rates, we adjust starting landscape parameters
######################################################################################################################
base_0 = .02
base_1 = .04
trend = -.005
ATT = -.01

# we'll need to recompute the parameters 
std_avp = (std_a^2+std_v^2+std_p^2)^.5
b0 = qnorm(base_0, mean = 0, sd = std_avp)
b1 = qnorm(base_1, mean = 0, sd = std_avp) - b0
b2_0 = qnorm(trend + base_0, mean = 0, sd = std_avp) - b0
b2_1 = qnorm(trend + base_1, mean = 0, sd = std_avp) - b0 - b1
b3 = qnorm( pnorm(b0+b1+b2_1, mean = 0, sd = std_avp) + ATT , mean = 0, sd = std_avp) - (b0 + b1 + b2_1)

aggregation_0.02 <- aggregate_complete(n, nobs, years, b0, b1, b2_0, b2_1, b3, std_a, std_v, std_p, cellsize, ppoints, cpoints)
summary_long_0.02 <- aggregation_0.02$summary_long



######################################################################################
# convert
#####################################################################################
TWFE_0.02 <- TWFE_fcn(n, nobs, years, b0, b1, b2_0, b2_1, b3, std_a, std_v)
TWFE_long_0.02 <- TWFE_0.02$summary_long

schart_df <- summary_long_0.02 %>%
  mutate(bias = as.numeric(bias),
         cover = as.numeric(cover))%>%
  group_by(pixel, grid, property, county, pixel.fe, grid.fe, property.fe, county.fe, treatment.fe, se_pixel, se_grid, se_property, se_county, notes) %>%
  summarise(q05 = quantile(bias, probs= .05),
            q95 = quantile(bias, probs= .95),
            rmse = rmse(bias, rep(0, length(bias))),
            bias = mean(bias),
            coverage = mean(cover)
  )%>%
  mutate_at(vars(pixel:se_county), .funs = as.numeric)%>%
  mutate_at(vars(pixel:se_county), .funs = as.logical)%>%
  select(bias, everything()) 

######################################################################################################################
#
######################################################################################################################

base_0 = .02
base_1 = .03
trend = -.005
ATT = -.01

# we'll need to recompute the parameters 
std_avp = (std_a^2+std_v^2+std_p^2)^.5
b0 = qnorm(base_0, mean = 0, sd = std_avp)
b1 = qnorm(base_1, mean = 0, sd = std_avp) - b0
b2_0 = qnorm(trend + base_0, mean = 0, sd = std_avp) - b0
b2_1 = qnorm(trend + base_1, mean = 0, sd = std_avp) - b0 - b1
b3 = qnorm( pnorm(b0+b1+b2_1, mean = 0, sd = std_avp) + ATT , mean = 0, sd = std_avp) - (b0 + b1 + b2_1)

# 
# aggregation_0.01 <- aggregate_complete(n, nobs, years, b0, b1, b2_0, b2_1, b3, std_a, std_v, std_p, cellsize, ppoints, cpoints)
# summary_long_0.01 <- aggregation_0.01$summary_long


TWFE_0.01 <- TWFE_fcn(n, nobs, years, b0, b1, b2_0, b2_1, b3, std_a, std_v)
TWFE_long_0.01 <- TWFE_0.01$summary_long

######################################################################################################################
#
######################################################################################################################

base_0 = .02
base_1 = .02
trend = -.005
ATT = -.01

# we'll need to recompute the parameters 
std_avp = (std_a^2+std_v^2+std_p^2)^.5
b0 = qnorm(base_0, mean = 0, sd = std_avp)
b1 = qnorm(base_1, mean = 0, sd = std_avp) - b0
b2_0 = qnorm(trend + base_0, mean = 0, sd = std_avp) - b0
b2_1 = qnorm(trend + base_1, mean = 0, sd = std_avp) - b0 - b1
b3 = qnorm( pnorm(b0+b1+b2_1, mean = 0, sd = std_avp) + ATT , mean = 0, sd = std_avp) - (b0 + b1 + b2_1)


# aggregation_0.00 <- aggregate_complete(n, nobs, years, b0, b1, b2_0, b2_1, b3, std_a, std_v, std_p, cellsize, ppoints, cpoints)
# summary_long_0.00 <- aggregation_0.00$summary_long


TWFE_0 <- TWFE_fcn(n, nobs, years, b0, b1, b2_0, b2_1, b3, std_a, std_v)
TWFE_long_0 <- TWFE_0$summary_long


######################################################################################################################
#
######################################################################################################################

base_0 = .05
base_1 = .02
trend = -.005
ATT = -.01

# we'll need to recompute the parameters 
std_avp = (std_a^2+std_v^2+std_p^2)^.5
b0 = qnorm(base_0, mean = 0, sd = std_avp)
b1 = qnorm(base_1, mean = 0, sd = std_avp) - b0
b2_0 = qnorm(trend + base_0, mean = 0, sd = std_avp) - b0
b2_1 = qnorm(trend + base_1, mean = 0, sd = std_avp) - b0 - b1
b3 = qnorm( pnorm(b0+b1+b2_1, mean = 0, sd = std_avp) + ATT , mean = 0, sd = std_avp) - (b0 + b1 + b2_1)

aggregation_n0.03 <- aggregate_complete(n, nobs, years, b0, b1, b2_0, b2_1, b3, std_a, std_v, std_p, cellsize, ppoints, cpoints)
summary_long_n0.03 <- aggregation_n0.03$summary_long


TWFE_n0.03 <- TWFE_fcn(n, nobs, years, b0, b1, b2_0, b2_1, b3, std_a, std_v)
TWFE_long_n0.03 <- TWFE_n0.03$summary_long


######################################################################################################################
#
######################################################################################################################


######################################################################################################################
#
######################################################################################################################

# bind dataframes together with different parameterizations and write summary_long as csv

# summary_long <- rbind(summary_long_0.00, summary_long_0.01, summary_long_0.02, summary_long_0.03)

TWFE_long <- rbind(TWFE_long_0, TWFE_long_0.01, TWFE_long_0.02, TWFE_long_0.03#, TWFE_long_n0.03
                   )%>%
  group_by(b0, b1, b2_0, b2_1, b3)%>%
  mutate(parameterization = cur_group_id())

write.csv(TWFE_long, "TWFE_long.csv")