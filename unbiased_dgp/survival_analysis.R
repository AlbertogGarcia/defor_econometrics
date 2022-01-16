library(ggplot2)
library(clubSandwich)
library(reshape2)
library(matrixStats)
library(ggplot2)
# library(plm)
library(Metrics)
library(DataCombine)
library(dplyr)
library(tidyverse)
library(tictoc)
library(fixest)
library(here)
library(DeclareDesign)
source(here::here('unbiased_dgp', 'full_landscape.R'))

library(survival)
library(ggplot2)
library(dplyr)
library(ggfortify)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Parameterization (largely taken from analysis, but longer time periods) ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
std_a = 0
std_v = 0.5
years = 10
nobs = 14400
n = 25

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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Create panel dataframe using code from aggregate_complete ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hr_list <- numeric(n)
ci_list <- numeric(n)
hr_cf_list <- numeric(n)

countyscape = full_landscape(nobs, cellsize, ppoints, cpoints)
pixloc_df = countyscape$pixloc_df

pixloc <- pixloc_df

# COMMENTING OUT LOOP TEMPORARILY ----
for(i in 1:n){
  tic("loop")
  
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
  
  
  
  # these_years <- unique(panels$year)
  # 
  # haz_ratios <- data.frame()
  # 
  # for(k in these_years){
  #   
  #   
  #   haz_ratios <- rbind(haz_ratios, 
  #                       data.frame("time" = k, 
  #                                  "observed_prob" = mean(subset(panels, year == k & treat == 1)$y_it, na.rm = TRUE),
  #                                  "counterfactual_prob" =
  #                                    mean(subset(panels, year == k & treat == 1)$y_it_counter, na.rm = TRUE)
  #                       )
  #   )
  #   
  # }
  # 
  # haz_ratios <- haz_ratios %>%
  #   mutate(haz_rat = observed_prob/counterfactual_prob,
  #          post = (time > max(these_years)/2)*1 )%>%
  #   filter(post == 1)
  # 
  # truth_list[i] <- mean(haz_ratios$haz_rat)
  
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
  
  
  ### In the DID world, this hazard ratio would be the ATT parameter, and we'd expect the DID to recover this
  ### Is there a reason not to expect this here?
  
  
  
  # # Create the new data
  # survplot_df <- with(surv_cf_df,
  #                data.frame(observed = c(0,1)
  #                )
  # )
  # res.cox <- coxph(Surv(t_start, t_end, outcome) ~ observed
  #                             , data = surv_cf_df)
  # # Survival curves
  # fit <- survfit(res.cox, data = surv_cf_df, newdata = survplot_df)
  # ggsurvplot(fit, conf.int = TRUE, 
  #            ggtheme = theme_minimal())
  # 
  # 
  # 
  # these_years <- unique(panels_long_cf$year)
  # 
  # surv <- data.frame()
  # 
  # for(k in these_years){
  # 
  #   surv <- rbind(surv,
  #                 data.frame("observed" = c(1, 0), "time" = k, "model" = "truth",
  #                            "surv" = c(1 - mean(subset(panels_long_cf, year == k & observed == 1)$defor_cf, na.rm = TRUE),
  #                                       1 - mean(subset(panels_long_cf, year == k & observed == 0)$defor, na.rm = TRUE))
  #                 )
  #   )
  # 
  # 
  # }
  # 
  # survival_df <- rbind(data.frame("observed" = c(rep(0, 20), rep(1, 20)), "time" = c(seq(from = 1, to = 20, by = 1), seq(from = 1, to = 20, by = 1)),
  #                           "surv" = c(fit$surv[,2], fit$surv[,1]), "model" = "cox ph fit"
  # ), surv)%>%
  #   mutate(observed = as.factor(observed))
  # 
  # ggplot(data = survival_df, aes(x = time, y = surv, color = observed, linetype = model))+
  #   geom_line(size = 1.2)+theme_minimal()
  # 
  # 
  
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
  ci_lwr <- confint(cox)[[2,1]] %>% exp()
  ci_upr <- confint(cox)[[2,2]] %>% exp()
  ci <- c(ci_lwr, ci_upr)
  
  hr_list[i] <- hr
  ci_list[i] <- ci
  print(i)
  toc()
}

mean(hr_list)
mean(hr_cf_list)

haz_rat <- pnorm(b0+b1+b2_1+b3, 0, (std_a^2+std_v^2 +std_p^2)^.5) / pnorm(b0+b1+b2_1, 0, (std_a^2+std_v^2 + std_p^2)^.5)
haz_rat



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

summary(cox_counterfactual)
haz_rat <- pnorm(b0+b1+b2_1+b3, 0, (std_a^2+std_v^2 +std_p^2)^.5) / pnorm(b0+b1+b2_1, 0, (std_a^2+std_v^2 + std_p^2)^.5)
haz_rat
haz_smpl <- (base_1 + trend + ATT) / (base_1 + trend)
haz_smpl


defor_summary <- panels %>%
  group_by(treat, post) %>%
  summarise(mean_y_it = mean(y_it, na.rm = TRUE),
            mean_y_it_counter = mean(y_it_counter, na.rm = TRUE))

defor_summary[4,3] / defor_summary[4,4]


# summary(cox)
# summary(cox_test)
# summary(cox_counterfactual)
# 
# 
# 
# zph <- cox.zph(cox)
# plot(zph)
# zph_test <- cox.zph(cox_test)
# plot(zph_test)
# zph_cf <- cox.zph(cox_counterfactual)
# plot(zph_cf)
# 
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
