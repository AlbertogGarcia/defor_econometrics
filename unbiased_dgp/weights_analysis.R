
source(here::here('unbiased_dgp', 'heterogeneous_propertyarea.R'))

# we start with our base parameterization without property level perturbations
std_a = 0.1
std_v = 0.5
std_p = 0.0
std_b3 = .05
years = 2
nobs = 120^2
n = 250

ppoints = 70
cellsize = 10
cpoints = 40

# here are the landscape characteristics in this parameterization
# note that the bias for the TWFE model will be equal to the pre-treatment difference in deforestation rtes, which is 0.03
base_0 = .02
base_1 = .05
trend = -.005
ATT = -.01

# we'll need to compute the parameters 
std_avp = (std_a^2+std_v^2+std_p^2)^.5
std_avpt = (std_a^2+std_v^2+std_p^2+std_b3^2)^.5
b0 = qnorm(base_0, mean = 0, sd = std_avp)
b1 = qnorm(base_1, mean = 0, sd = std_avp) - b0
b2_0 = qnorm(trend + base_0, mean = 0, sd = std_avp) - b0
b2_1 = qnorm(trend + base_1, mean = 0, sd = std_avp) - b0 - b1

set.seed(0930)
testing_weights <- heterogeneous_propertyarea(n, nobs, years, b0, b1, b2_0, b2_1, std_a, std_v, std_p, std_b3, given_ATT = ATT, cellsize, ppoints, cpoints, rm.selection = TRUE)
summary_pweights <- testing_weights$summary_long

library(rio)
export(summary_pweights, "summary_pweights.rds")
