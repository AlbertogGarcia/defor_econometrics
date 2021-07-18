library(tidyverse)
library(tictoc)
library(fixest)
library(here)
library(DeclareDesign)
library(did2s)
source(here::here('multigroup_dgp', 'multi_group_landscape.R'))
source(here::here('multigroup_dgp', 'my_event_study.R'))
source(here::here('multigroup_dgp', 'event_simulation.R'))

#begin function
multipleGT_parallel <- function(n, nobs, base_a, base_b, base_c, trend1, trend2, trend3, ATT, dyn_ATT = 0, std_a = 0.0, std_v = 0.25, std_p = 0.0, cellsize, ppoints, cpoints, num_cores=1){
  
  countyscape = multi_group_landscape(nobs, cellsize, ppoints, cpoints)
  pixloc_df = countyscape$pixloc_df
  pixloc <- pixloc_df
  
  #create data frame with 0 rows and 5 columns
  pixel_es_long <- data.frame(matrix(ncol = 5, nrow = 0))
  county_es_long <- pixel_es_long
  
  #provide column names
  cnames <- c('term', 'estimate', 'std.error', 'estimator', 'iteration')
  colnames(pixel_es_long) <- cnames
  colnames(county_es_long) <- cnames
  
  
  std_avp = (std_a^2+std_v^2+std_p^2)^.5
  
  
  ################################################################################################
  ####### early group
  ################################################################################################
  b0a = qnorm(base_a, mean = 0, sd = std_avp)
  
  b1a = qnorm(trend1 + base_a, mean = 0, sd = std_avp) -b0a
  b2a = qnorm(trend2 + pnorm(b0a+b1a, mean = 0, sd = std_avp), mean = 0, sd = std_avp) -b1a - b0a
  b3a = qnorm(trend3 + pnorm(b0a+b1a+b2a, mean = 0, sd = std_avp), mean = 0, sd = std_avp) -b1a - b0a - b2a
  b4a = qnorm(trend3 + pnorm(b0a+b1a+b2a+b3a, mean = 0, sd = std_avp), mean = 0, sd = std_avp) -b1a - b0a - b2a - b3a
  
  tau_a = qnorm( pnorm(b0a+b1a+b2a, mean = 0, sd = std_avp) + ATT , mean = 0, sd = std_avp) - (b0a+b1a)
  tau_a2 = qnorm( pnorm(b0a+b1a+b2a+b3a+tau_a, mean = 0, sd = std_avp) + dyn_ATT , mean = 0, sd = std_avp) - (b0a+b1a+b2a+b3a+tau_a)
  tau_a3 = qnorm( pnorm(b0a+b1a+b2a+b3a+b4a+tau_a2, mean = 0, sd = std_avp) + dyn_ATT , mean = 0, sd = std_avp) - (b0a+b1a+b2a+b3a++b4a+tau_a2)
  
  ################################################################################################
  ####### late group
  ################################################################################################
  b0b = qnorm(base_b, mean = 0, sd = std_avp)
  
  b1b = qnorm(trend1 + base_b, mean = 0, sd = std_avp) -b0b
  b2b = qnorm(trend2 + pnorm(b0b+b1b, mean = 0, sd = std_avp), mean = 0, sd = std_avp) -b1b - b0b
  b3b = qnorm(trend3 + pnorm(b0b+b1b+b2b, mean = 0, sd = std_avp), mean = 0, sd = std_avp) -b1b - b0b - b2b
  b4b = qnorm(trend3 + pnorm(b0b+b1b+b2b+b3b, mean = 0, sd = std_avp), mean = 0, sd = std_avp) -b1b - b0b - b2b - b3b
  
  tau_b = qnorm( pnorm(b0b+b1b+b2b+b3b, mean = 0, sd = std_avp) + ATT , mean = 0, sd = std_avp) - (b0b+b1b+b2b+b3b)
  tau_b2 = qnorm( pnorm(b0b+b1b+b2b+b3b+b4b+tau_b, mean = 0, sd = std_avp) + dyn_ATT , mean = 0, sd = std_avp) - (b0b+b1b+b2b+b3b+b4b+tau_b)
  ################################################################################################
  ###########      never treated group
  ################################################################################################
  b0c = qnorm(base_c, mean = 0, sd = std_avp)
  
  b1c = qnorm(trend1 + base_c, mean = 0, sd = std_avp) -b0c
  b2c = qnorm(trend2 + pnorm(b0c+b1c, mean = 0, sd = std_avp), mean = 0, sd = std_avp) -b1c - b0c
  b3c = qnorm(trend3 + pnorm(b0c+b1c+b2c, mean = 0, sd = std_avp), mean = 0, sd = std_avp) -b1c - b0c - b2c
  b4c = qnorm(trend3 + pnorm(b0c+b1c+b2c+b3c, mean = 0, sd = std_avp), mean = 0, sd = std_avp) -b1c - b0c - b2c - b3c
  
  #--- number of polygons in a group ---#
  iterations_per_core <- floor(n/num_cores)
  
  library(doParallel)
  library(snow)
  
  x <- mclapply(1:num_cores,
                     function(i) event_simulation(iterations_per_core, pixloc, std_a, std_v, std_p, b0a, b1a, b2a, b3a, b4a, b0b, b1b, b2b, b3b, b4b, b0c, b1c, b2c, b3c, b4c, tau_a, tau_a2, tau_a3, tau_b, tau_b2)
                , mc.cores = num_cores)
  
  es_long <- x %>%
    rbindlist(idcol = TRUE)
  
  return{es_long}
  
}
