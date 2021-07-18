source(here::here('multigroup_dgp', 'multipleGT.R'))
base_a = .07
base_b = .04
base_c = .02
trend1 = .00
trend2 = .00
trend3 = .00
ATT = -.02
dyn_ATT = 0

std_a = 0.0
std_v = 0.25
std_p = 0.0

cpoints = 50
cellsize=10 
ppoints=50

nobs = 10000
n=100

multiGT <- multipleGT(n, nobs, base_a, base_b, base_c, trend1, trend2, trend3, ATT, dyn_ATT = 0, std_a = 0.0, std_v = 0.25, std_p = 0.0, cellsize=10, ppoints=50, cpoints)
  
library(rio)
export(multiGT$county_es_long, "county_es_long.rds")
export(multiGT$pixel_es_long, "pixel_es_long.rds")

county_es <- multiGT$county_es_long %>%
  group_by(term, estimator)%>%
  mutate(estimate = as.numeric(estimate))%>%
  summarise(q05 = quantile(estimate, probs = 0.05),
            q95 = quantile(estimate, probs = 0.95),
            estimate = mean(estimate))

pixel_es <- multiGT$pixel_es_long %>%
  group_by(term, estimator)%>%
  summarise(q05 = quantile(estimate, probs = 0.05),
            q95 = quantile(estimate, probs = 0.95),
            estimate = mean(estimate))

num_cores = 2
iterations_per_core=1
library(snow)

clus <- makeCluster(num_cores)
# character vector of needed pieces
these_vars <- c("iterations_per_core", "pixloc", "std_a", "std_p", "b0a", "b1a", "b2a", "b3a", "b4a", "b0b", "b1b", "b2b", "b3b", "b4b", "b0c", "b1c", "b2c", "b3c", "b4c", "tau_a", "tau_a2", "tau_a3", "tau_b", "tau_b2")
# need to export all objects to each node
clusterExport(clus,these_vars)


x <- parLapply(clus,
               1:num_cores,
       event_simulation(iterations_per_core, pixloc, std_a, std_p, b0a, b1a, b2a, b3a, b4a, b0b, b1b, b2b, b3b, b4b, b0c, b1c, b2c, b3c, b4c, tau_a, tau_a2, tau_a3, tau_b, tau_b2)
)
       
       