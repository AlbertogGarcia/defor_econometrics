source(here::here('unbiased_dgp', 'aggregate_complete.R')) #aggreagate_complete is the function we'll use most, as it runs all combinations and generates summary dataframe to be used for specification chart
source(here::here('unbiased_dgp', 'TWFE_expost.R'))
source(here::here('unbiased_dgp', 'TWFE_fcn.R'))

library(tidyverse)

base_0 = .02
base_1 = .05
trend = -.005
ATT = -.01

std_a = 0
std_v = 0.5
years = 2
nobs = 100^2
n = 500

cellsize = 10
ppoints = 70
cpoints = 40

###############################################################################################################
######## set seed
###############################################################################################################

set.seed(0930)

###############################################################################################################
######## Baseline showing aggregation resolves pixel fixed effects issue
###############################################################################################################
std_p = 0
std_avp = (std_a^2+std_v^2+std_p^2)^.5
b0 = qnorm(base_0, mean = 0, sd = std_avp)
b1 = qnorm(base_1, mean = 0, sd = std_avp) - b0
b2_0 = qnorm(trend + base_0, mean = 0, sd = std_avp) - b0
b2_1 = qnorm(trend + base_1, mean = 0, sd = std_avp) - b0 - b1
b3 = qnorm( pnorm(b0+b1+b2_1, mean = 0, sd = std_avp) + ATT , mean = 0, sd = std_avp) - (b0 + b1 + b2_1)

# summary function that estimates all of the different specifications
aggregation <- aggregate_complete(n, nobs, years, b0, b1, b2_0, b2_1, b3, std_a, std_v, std_p, cellsize, ppoints, cpoints, rm.selection = FALSE)

summary_long <- aggregation$summary_long

library(rio)
export(summary_long, "summary_long.rds")


###############################################################################################################
######## show TWFE is equivalent to dropping all pixels deforested in first period
###############################################################################################################
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

#ATT = pnorm(b0+b1+b2_1+b3, 0, (std_a^2+std_v^2 )^.5) - pnorm(b0+b1+b2_1, 0, (std_a^2+std_v^2 )^.5)

estimator_comp <- TWFE_expost(n, nobs, years, b0, b1, b2_0, b2_1, b3, std_a, std_v)

summary_coeff <- estimator_comp$summary_long %>%
  mutate_at(vars(bias), as.numeric)

summary_wide  <- summary_coeff %>%
    group_by(model)%>%
    summarise(RMSE = rmse(bias, 0),
              q25 = quantile(bias, probs = .25),
              q75 = quantile(bias, probs = .75),
              Bias = mean(bias))

export(summary_wide, "twfe_comp.rds")

###############################################################################################################
######## Introduction pixel level unobservables, which impact non-random selection
###############################################################################################################
std_a = 0.25
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
aggregation_sel <- aggregate_complete(n, nobs, years, b0, b1, b2_0, b2_1, b3, std_a, std_v, std_p, cellsize, ppoints, cpoints, rm.selection = FALSE)

summary_selection <- aggregation_sel$summary_long 

export(summary_selection, "summary_selection.rds")
###############################################################################################################
######## Adding in property level disturbances
###############################################################################################################

std_a = 0.1

#### 0.0
std_p = 0

# we'll need to compute the parameters 
std_avp = (std_a^2+std_v^2+std_p^2)^.5
b0 = qnorm(base_0, mean = 0, sd = std_avp)
b1 = qnorm(base_1, mean = 0, sd = std_avp) - b0
b2_0 = qnorm(trend + base_0, mean = 0, sd = std_avp) - b0
b2_1 = qnorm(trend + base_1, mean = 0, sd = std_avp) - b0 - b1
b3 = qnorm( pnorm(b0+b1+b2_1, mean = 0, sd = std_avp) + ATT , mean = 0, sd = std_avp) - (b0 + b1 + b2_1)

#ATT = pnorm(b0+b1+b2_1+b3, 0, (std_a^2+std_v^2 )^.5) - pnorm(b0+b1+b2_1, 0, (std_a^2+std_v^2 )^.5)

# summary function that estimates all of the different specifications
aggregation_0 <- aggregate_complete(n, nobs, years, b0, b1, b2_0, b2_1, b3, std_a, std_v, std_p, cellsize, ppoints, cpoints, rm.selection = FALSE)

summary_long_0 <- aggregation_0$summary_long

#### 0.1
std_p = 0.1

# we'll need to compute the parameters 
std_avp = (std_a^2+std_v^2+std_p^2)^.5
b0 = qnorm(base_0, mean = 0, sd = std_avp)
b1 = qnorm(base_1, mean = 0, sd = std_avp) - b0
b2_0 = qnorm(trend + base_0, mean = 0, sd = std_avp) - b0
b2_1 = qnorm(trend + base_1, mean = 0, sd = std_avp) - b0 - b1
b3 = qnorm( pnorm(b0+b1+b2_1, mean = 0, sd = std_avp) + ATT , mean = 0, sd = std_avp) - (b0 + b1 + b2_1)

#ATT = pnorm(b0+b1+b2_1+b3, 0, (std_a^2+std_v^2 )^.5) - pnorm(b0+b1+b2_1, 0, (std_a^2+std_v^2 )^.5)

aggregation_1 <- aggregate_complete(n, nobs, years, b0, b1, b2_0, b2_1, b3, std_a, std_v, std_p, cellsize, ppoints, cpoints, rm.selection = FALSE)

summary_long_1 <- aggregation_1$summary_long

##### 0.2

std_p = 0.25

# we'll need to compute the parameters 
std_avp = (std_a^2+std_v^2+std_p^2)^.5
b0 = qnorm(base_0, mean = 0, sd = std_avp)
b1 = qnorm(base_1, mean = 0, sd = std_avp) - b0
b2_0 = qnorm(trend + base_0, mean = 0, sd = std_avp) - b0
b2_1 = qnorm(trend + base_1, mean = 0, sd = std_avp) - b0 - b1
b3 = qnorm( pnorm(b0+b1+b2_1, mean = 0, sd = std_avp) + ATT , mean = 0, sd = std_avp) - (b0 + b1 + b2_1)

#ATT = pnorm(b0+b1+b2_1+b3, 0, (std_a^2+std_v^2 )^.5) - pnorm(b0+b1+b2_1, 0, (std_a^2+std_v^2 )^.5)

aggregation_25 <- aggregate_complete(n, nobs, years, b0, b1, b2_0, b2_1, b3, std_a, std_v, std_p, cellsize, ppoints, cpoints, rm.selection = FALSE)

summary_long_25 <- aggregation_25$summary_long

# 
# #### 0.3
# 
# std_p = 0.3
# std_avp = (std_a^2+std_v^2+std_p^2)^.5
# b0 = qnorm(base_0, mean = 0, sd = std_avp)
# b1 = qnorm(base_1, mean = 0, sd = std_avp) - b0
# b2_0 = qnorm(trend + base_0, mean = 0, sd = std_avp) - b0
# b2_1 = qnorm(trend + base_1, mean = 0, sd = std_avp) - b0 - b1
# b3 = qnorm( pnorm(b0+b1+b2_1, mean = 0, sd = std_avp) + ATT , mean = 0, sd = std_avp) - (b0 + b1 + b2_1)
# 
# aggregation_5 <- aggregate_complete(n, nobs, years, b0, b1, b2_0, b2_1, b3, std_a, std_v, std_p, cellsize, ppoints, cpoints, rm.selection = TRUE)
# 
# summary_long_5 <- aggregation_5$summary_long

summary_full <- rbind(summary_long_0, summary_long_1, summary_long_25)

export(summary_full, "summary_full2.rds")



full_summary <- summary_full %>%
  mutate_at(vars(bias, cover), as.numeric
  )%>%
  mutate_at(vars(pixel, grid, property, county, pixel.fe, grid.fe, property.fe, county.fe, treatment.fe, weights, se_pixel, se_grid, se_property, se_county), ~ as.logical(as.integer(.))
  )%>%
  filter(is.na(notes) & pixel.fe == FALSE)%>%
  group_by(std_p, pixel, grid, property, county, pixel.fe, grid.fe, property.fe, county.fe, treatment.fe, weights, se_pixel, se_grid, se_property, se_county)%>%
  summarise(cover = mean(cover))%>%
  spread(std_p, cover)

###############################################################################################################
######## constructing coverage tables
###############################################################################################################

summary <- readRDS("summary_full2.rds") %>%
  group_by(std_p, pixel, grid, property, county, pixel.fe, grid.fe, property.fe , county.fe, treatment.fe, weights, se_pixel, se_property, se_county, se_grid, notes) %>%
  summarise(std_p = mean(std_p),
            cover = mean(as.numeric(cover))
            )%>%
  ungroup()%>%
  filter(is.na(notes))

did_cover <- summary %>%
  filter(pixel == 1 & treatment.fe == 1)%>%
  mutate(`unit of analysis` = "pixel",
         `fixed effects` = "treatment",
         se = ifelse(se_pixel == 1, "clustered at pixel", ifelse(se_county == 1, "clustered at county", ifelse(se_property == 1, "clustered at property", ifelse(se_grid == 1, "clustered at grid", NA))))
  )%>%
  spread(key = std_p,
         value = cover)%>%
  select(`unit of analysis`, `fixed effects`, se, `0`, `0.1`, `0.25`)

fix_cover <- summary %>%
  filter(pixel==1 & treatment.fe != 1 & pixel.fe != 1)%>%
  mutate(`unit of analysis` = "pixel",
         `fixed effects` = ifelse(county.fe == 1, "county", ifelse(property.fe == 1, "property", ifelse(grid.fe == 1, "grid", NA))),
         se = ifelse(se_pixel == 1, "clustered at pixel", ifelse(se_county == 1, "clustered at county", ifelse(se_property == 1, "clustered at property", ifelse(se_grid == 1, "clustered at grid", NA))))
  )%>%
  spread(key = std_p,
         value = cover)%>%
  select(`unit of analysis`, `fixed effects`, se, `0`, `0.1`, `0.25`)

agg_cover <- summary %>%
  filter(property==1 | county ==1 | grid == 1)%>%
  mutate(`unit of analysis` = ifelse(county == 1, "county", ifelse(property == 1, "property", ifelse(grid == 1, "grid", NA))),
         `fixed effects` = ifelse(county.fe == 1, "county", ifelse(property.fe == 1, "property", ifelse(grid.fe == 1, "grid", NA))),
         weighted = ifelse(weights == 1, "yes", "no"),
         se = ifelse(se_pixel == 1, "clustered at pixel", ifelse(se_county == 1, "clustered at county", ifelse(se_property == 1, "clustered at property", ifelse(se_grid == 1, "clustered at grid", NA))))
  )%>%
  spread(key = std_p,
         value = cover)%>%
  select(`unit of analysis`, `fixed effects`, weighted, se, `0`, `0.1`, `0.25`)

write.csv(agg_cover, "agg_cover.csv")
write.csv(fix_cover, "fix_cover.csv")
write.csv(did_cover, "did_cover.csv")


