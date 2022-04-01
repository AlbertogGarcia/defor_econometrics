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
std_v = 0.5
years = 10
nobs = 100^2
n = 500

cellsize = 10
ppoints = 70
std_a = 0
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

set.seed(0930)

survival_full <- all_specifications(n, nobs, years, b0, b1, b2_0, b2_1, b3, std_a, std_v, std_p, cellsize, ppoints, cpoints)
survival_long <- survival_full$summary_long

# survival <- survival_did(n, nobs, years, b0, b1, b2_0, b2_1, b3, std_a, std_v, std_p, cellsize, ppoints, cpoints)
# 
# survival_long <- survival$summary_long

library(rio)
export(survival_long, "survival_long.rds")

survival_summary <- survival_long %>%
  group_by(model)%>%
  summarize(HRR = mean(as.numeric(HRR)),
            Bias = mean(as.numeric(bias)),
            RMSE = rmse(as.numeric(bias), 0),
            cover = mean(as.numeric(cover), na.rm = TRUE)
  )











#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Create panel dataframe using code from aggregate_complete ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

countyscape = full_landscape(nobs, cellsize, ppoints, cpoints)
pixloc_df = countyscape$pixloc_df

pixloc <- pixloc_df

Nobs <- length(pixloc$treat)  
panels <- fabricate(
  pixels = add_level(N = Nobs, a_i = rnorm(N, 0, std_a), treat = pixloc$treat),
  year = add_level(N = (years*2), nest = FALSE),
  obs = cross_levels(
    by = join(pixels, year),
    post = ifelse(year > years, 1, 0),
    v_it = rnorm(N, 0, std_v),
    ystar = b0 + b1*treat + b2_0*post*(1-treat) + b2_1*post*treat + b3*treat*post + a_i + v_it,
    ystar_cf = b0 + b1*treat + b2_0*post*(1-treat) + b2_1*post*treat + a_i + v_it
  )
)

#generate random 
errortable_property <- data.frame(property = as.character(unique(pixloc$property)), p_err = rnorm(length(unique(pixloc$property)), 0, std_p))
#errortable_county <- data.frame(county = as.character(unique(pixloc$county)), c_err = rnorm(length(unique(pixloc$county)), 0, std_c))

panels$pixels <- gsub("(?<![0-9])0+", "", panels$pixels, perl = TRUE)

panels <- panels %>%
  inner_join(pixloc, by = c("pixels", "treat")) %>%
  inner_join(errortable_property, by = "property") %>%
  #inner_join(errortable_county, by = "county") %>%
  mutate(ystar = ystar + p_err,
         ystar_cf = ystar_cf + p_err,
         y = (ystar > 0)*1 ,
         y_cf = (ystar_cf > 0)*1 )


panels <- panels %>%
  mutate(pixels = as.numeric(pixels),
         year = as.numeric(year),
         defor_indic = ifelse(y==1, year, 99),
         defor_indic_cf = ifelse(y_cf==1, year, 99))%>%
  group_by(pixels)%>%
  mutate(defor_year = min(defor_indic),
         defor_year_cf = min(defor_indic_cf),
         defor = (year>=defor_year)*1,
         y_it = ifelse(year>defor_year, NA, defor),
         defor_cf = (year>=defor_year_cf)*1,
         y_it_counter = ifelse(year>defor_year_cf, NA, defor_cf)
  )%>%
  ungroup()


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Create counterfactual panel dataframe ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# this dataframe only contains treated units, both with and without treatment

panels_long_cf <- subset(panels, treat == 1)%>%
  gather(key = "group", value = "outcome", c(y_it, y_it_counter))%>%
  mutate(observed = ifelse(group == "y_it", 1, 0)) # observed = 1 means that treatment did actually happen



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Set up survival dataframe ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Ref: https://cran.r-project.org/web/packages/survival/vignettes/timedep.pdf
surv_df <- panels %>% 
  mutate(t_start = year - 1,
         t_end = year,
         # t_end = ifelse((t_end==20) & (y_it==0), Inf, t_end),
         outcome = y_it,
         treat_now = treat * post) %>% 
  select(pixels, t_start, t_end, outcome, treat, treat_now, post, year) %>% 
  drop_na()


### counterfactual survival dataframe
surv_cf_df <- panels_long_cf %>% 
  mutate(t_start = year - 1,
         t_end = year,
         # t_end = ifelse((t_end==20) & (y_it==0), Inf, t_end),
         outcome = outcome,
         treat_now = observed * post) %>% 
  select(pixels, t_start, t_end, outcome, treat, group, treat_now, post, year, observed) %>% 
  drop_na()




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Run Cox proportional hazards model ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# counterfactual cox proportional hazards model ----

# we use the counterfactual survival dataframe, which only includes the treated observations and either their treated or untreated outcomes
cox_counterfactual <- coxph(Surv(t_start, t_end, outcome) ~ treat_now 
                            , data = surv_cf_df)


#ggforest(cox)
print(summary(cox_counterfactual))
coefs <- cox_counterfactual$coefficients
hr_cf <- coefs[[1]] %>% exp()
hr_cf_list[i] <- hr_cf
# this hr corresponds really well to what we'd expect if the hazard ratio we want is 
haz_rat <- pnorm(b0+b1+b2_1+b3, 0, (std_a^2+std_v^2 +std_p^2)^.5) / pnorm(b0+b1+b2_1, 0, (std_a^2+std_v^2 + std_p^2)^.5)
haz_rat


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# running main cox model ----

# km_AG_fit <- survfit(Surv(tstart, tend, outcome) ~ treat + treat_now, data=surv_df)
# autoplot(km_AG_fit)
cox <- coxph(Surv(t_start, t_end, outcome) ~ treat + treat_now 
             , data = surv_df)


#ggforest(cox)
print(summary(cox))

coefs <- cox$coefficients
hr <- coefs[[2]] %>% exp()

## Testing conversion of coefficients from multiple cox regressions into
## hazard ratio of interest

surv_df_correction <- surv_df %>% 
  mutate(t_start = ifelse(post==1, t_start -10, t_start),
         t_end = ifelse(post==1, t_end -10, t_end))

# cox_pre <- coxph(Surv(t_start, t_end, outcome) ~ treat
#                         , data = surv_df_correction %>% filter(post == 0))
# summary(cox_pre)

cox_post <- coxph(Surv(t_start, t_end, outcome) ~ treat
                  , data = surv_df_correction %>% filter(post==1))
summary(cox_post)
hr_11_01 <- cox_post$coefficients %>% exp()

cox_treat <- coxph(Surv(t_start, t_end, outcome) ~ post
                   , data = surv_df_correction %>% filter(treat==1))
summary(cox_treat)
hr_11_10 <- cox_treat$coefficients %>% exp()

cox_notreat <- coxph(Surv(t_start, t_end, outcome) ~ post
                     , data = surv_df_correction %>% filter(treat==0))
summary(cox_notreat)
hr_01_00 <- cox_notreat$coefficients %>% exp()

hr_11_cf <- 1/(1/hr_11_10 + 1/hr_11_01 - (1/(hr_11_01*hr_01_00)))
hr_11_cf

hr_2 <- hr_11_10/hr_01_00

stratcox_treat <- coxph(Surv(t_start, t_end, outcome) ~  treat_now + post + strata(as.factor(treat))
                        , data = surv_df_correction )
summary(stratcox_treat)

summary(cox_counterfactual)
haz_rat <- pnorm(b0+b1+b2_1+b3, 0, (std_a^2+std_v^2 +std_p^2)^.5) / pnorm(b0+b1+b2_1, 0, (std_a^2+std_v^2 + std_p^2)^.5)
haz_rat
haz_smpl <- (base_1 + trend + ATT) / (base_1 + trend)
haz_smpl

## this model is more similar to the traditional DID setup but does not actually recover the desired HRR
cox_interaction <- coxph(Surv(t_start, t_end, outcome) ~ post*treat
                         , data = surv_df_correction )
summary(cox_interaction)

defor_summary <- panels %>%
  group_by(treat, post) %>%
  summarise(mean_y_it = mean(y_it, na.rm = TRUE),
            mean_y_it_counter = mean(y_it_counter, na.rm = TRUE))

defor_summary[4,3] / defor_summary[4,4]


# converting into ATT ---- 
d_obs = defor_summary[4,3]$mean_y_it
ATT_hat = d_obs - d_obs/hr_11_cf 


# summary(cox)
# summary(cox_test)
# summary(cox_counterfactual)
# 
# 

zph <- cox.zph(cox_int)
print(zph)
zph_test <- cox.zph(cox_test)
plot(zph_test)
zph_cf <- cox.zph(cox_counterfactual)
print(zph_cf)
plot(zph_cf)

# 
# 
# 
# surv_df <- surv_df %>% 
#   arrange(pixels, year)
# 
# cox <- coxph(Surv(t_start, t_end, outcome) ~ treat + treat_now + year 
#              , data = surv_df)
# summary(cox)
# 
# 
# 
# 
# surv_df_test <- surv_df %>% 
#   mutate(treat_level = ifelse(post==1 & treat==1, 4, ifelse(post==1 & treat==0, 2, ifelse(post==0 & treat==1, 3, 1))),
#          treat_level = as.factor(treat_level))
# cox <- coxph(Surv(t_start, t_end, outcome) ~ treat_level
#              , data = surv_df_test)
# summary(cox)
# 
# 
# 
# # base_0 = .02
# # base_1 = .05
# # trend = 0
# # ATT = -.01
# #
# 
# 
# panels <- panels %>% 
#   arrange(pixels, year)
# 
# 
# defor_summary <- panels %>% 
#   group_by(treat, post) %>% 
#   # group_by(treat, post, year) %>% 
#   summarise(mean_y_it = mean(y_it, na.rm = TRUE),
#             mean_y_it_counter = mean(y_it_counter, na.rm = TRUE))
# 
# defor_summary
# summary(cox)
# summary(cox_counterfactual)