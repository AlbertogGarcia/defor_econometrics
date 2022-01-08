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
library(ranger)
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
n = 100

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

rm.selection = FALSE

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Create panel dataframe using code from aggregate_complete ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hr_list <- numeric(n)

countyscape = full_landscape(nobs, cellsize, ppoints, cpoints)
pixloc_df = countyscape$pixloc_df


ATT <- pnorm(b0+b1+b2_1+b3, 0, (std_a^2+std_v^2 +std_p^2)^.5) - pnorm(b0+b1+b2_1, 0, (std_a^2+std_v^2 + std_p^2)^.5)

pixloc <- pixloc_df

did_covermat <- matrix(nrow = n, ncol = 5)
agg_covermat <- matrix(nrow = n, ncol = 3)
weight_covermat <- matrix(nrow = n, ncol = 2)
bad_covermat <- matrix(nrow = n, ncol = 2)
fix_covermat <- matrix(nrow = n, ncol = 3)


coeffmatrix <- matrix(nrow = n, ncol = 3)
weight_coeffmatrix <- matrix(nrow = n, ncol = 3)
bad_coeffmatrix <- matrix(nrow = n, ncol = 2)
fix_coeffmatrix <- matrix(nrow = n, ncol = 4)

selection_bias <- data.frame('iteration' = rep(NA, n), 
                             'sample_sel_bias' = rep(NA, n))

n_mod = 14
summ_row <- n_mod * n

summary_long <- data.frame('b0'= rep(b0, summ_row), 'b1'= rep(b1, summ_row), 'b2_0'= rep(b2_0, summ_row), 'b2_1'= rep(b2_1, summ_row), 'b3'= rep(b3, summ_row), 
                           'std_a'= rep(std_a, summ_row), 'std_v'= rep(std_v, summ_row), 'std_p'= rep(std_p, summ_row),
                           'iteration' = rep(NA, summ_row), 
                           'pixel'=rep(NA, summ_row),'grid'=rep(NA, summ_row),'property'=rep(NA, summ_row),'county'=rep(NA, summ_row),
                           'pixel fe'=rep(NA, summ_row),'grid fe'=rep(NA, summ_row),'property fe'=rep(NA, summ_row),'county fe'=rep(NA, summ_row),'treatment fe'=rep(NA, summ_row),
                           'weights'=rep(NA, summ_row),
                           'se_pixel'=rep(NA, summ_row), 'se_grid'=rep(NA, summ_row), 'se_property'=rep(NA, summ_row), 'se_county'=rep(NA, summ_row),
                           'bias'=rep(NA, summ_row), 'cover'=rep(NA, summ_row),'notes'=rep(NA, summ_row),
                           stringsAsFactors=FALSE)

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
  
  
  condition_1 <- subset(panels, treat==1)%>%
    group_by(pixels)%>%
    mutate(conditional = (y== 0)*1*(1-post))%>%
    filter(max(conditional==1) & post ==1)%>%
    mutate(
      if_term_1 = (y==1)*1*post
    )
  
  term_1 = mean(condition_1$if_term_1)
  
  condition_2 <- subset(panels, treat==0 )%>%
    group_by(pixels)%>%
    mutate(conditional = (y== 0)*1*(1-post))%>%
    filter(max(conditional==1) & post ==1)%>%
    mutate(
      if_term_2 = (y== 1)*1*post
    )
  
  term_2 = mean(condition_2$if_term_2)
  
  term_3 = mean( subset(panels, treat==1&post==1)$y
  )
  
  term_4 = mean( subset(panels, treat==0&post==1)$y
  )
  
  sel_bias = term_1 - term_2 -(term_3 - term_4)
  
  
  selection_bias[i, 1] <- i
  selection_bias[i, 2] <- sel_bias
  
  ATT = mean( subset(panels, treat==1&post==1)$y)-mean( subset(panels, treat==1&post==1)$y_cf)
  ATT = pnorm(b0+b1+b2_1+b3, 0, (std_a^2+std_v^2 +std_p^2)^.5) - pnorm(b0+b1+b2_1, 0, (std_a^2+std_v^2 + std_p^2)^.5)
  
  if(rm.selection){
    ATT = ATT+sel_bias  
  }
  
  #need to determine which year deforestation occurred
  year_df <- panels %>%
    dplyr::select(pixels, year, y) %>%
    dcast(pixels ~ year , value.var = "y")
  
  rownames(year_df) <- year_df$pixels
  
  year_df <- year_df %>%
    dplyr::select(- pixels)
  
  #creating variable for the year a pixel is deforested
  not_defor <- rowSums(year_df)<1 *1
  defor_year <- max.col(year_df, ties.method = "first") 
  defor_df <- transform(year_df, defor_year = ifelse(not_defor==1, years*2+1, defor_year))
  defor_df <- tibble::rownames_to_column(defor_df)
  names(defor_df)[1] <- paste("pixels")
  
  panels <- defor_df %>%
    dplyr::select(pixels, defor_year) %>%
    inner_join(panels, by = "pixels")
  
  
  cols.num <- c("pixels", "grid", "property", "county", "year")
  panels[cols.num] <- sapply(panels[cols.num],as.numeric)
  
  panels <- panels %>%
    mutate(indic = year - defor_year) %>%
    mutate(defor = ifelse(indic > 0, 1, y))%>%
    mutate(y_it = ifelse(indic > 0, NA, y))
    
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # Set up survival dataframe ----
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  surv_df <- panels %>% 
    mutate(t_start = year - 1,
           t_end = year,
           # t_end = ifelse((t_end==20) & (y_it==0), Inf, t_end),
           outcome = y_it,
           treat_now = treat * post) %>% 
    select(pixels, t_start, t_end, outcome, treat, treat_now, post) %>% 
    drop_na()
  
  
  ## Old structure
  # y_defor <- panels %>% 
  #   filter(y_it == 1) %>% 
  #   mutate(year_defor = year) %>% 
  #   select(pixels, year_defor)
  # 
  # surv_df <- panels %>% 
  #   select(pixels, treat) %>% 
  #   distinct() %>% 
  #   merge(y_defor, by = "pixels", all = TRUE)
  # 
  # surv_df <- surv_df %>% 
  #   mutate(censor = is.na(year_defor),
  #          year_defor = replace_na(year_defor, years * 2))
  
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # Run Cox proportional hazards model ----
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # km_AG_fit <- survfit(Surv(tstart, tend, outcome) ~ treat + treat_now, data=surv_df)
  # autoplot(km_AG_fit)
  cox <- coxph(Surv(t_start, t_end, outcome) ~ treat + treat_now, data = surv_df)
  print(summary(cox))
  coefs <- cox$coefficients
  hr <- coefs[[2]] %>% exp()
  
  hr_list[i] <- hr
  toc()
}
