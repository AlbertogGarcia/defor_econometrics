#this function uses the property_scape gen fcn to generate a landscape with properties. Then we introduce property level perturbations and compare estimates when the data is aggregateed to the grid vs. property level

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
library(msm)
library(survival)

source(here::here('unbiased_dgp', 'full_landscape.R'))

#begin function
all_specifications <- function(n, nobs, years, b0, b1, b2_0, b2_1, b3, std_a = 0.1, std_v = 0.25, std_p = 0.0, cellsize, ppoints, cpoints){
  
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
  
  n_mod = 16
  summ_row <- n_mod * n
  
  summary_long <- data.frame('b0'= rep(b0, summ_row), 'b1'= rep(b1, summ_row), 'b2_0'= rep(b2_0, summ_row), 'b2_1'= rep(b2_1, summ_row), 'b3'= rep(b3, summ_row), 
                             'std_a'= rep(std_a, summ_row), 'std_v'= rep(std_v, summ_row), 'std_p'= rep(std_p, summ_row), "years" = years,
                             'iteration' = rep(NA, summ_row), 
                             'pixel'=rep(NA, summ_row),'grid'=rep(NA, summ_row),'property'=rep(NA, summ_row),'county'=rep(NA, summ_row),
                             'pixel fe'=rep(NA, summ_row),'grid fe'=rep(NA, summ_row),'property fe'=rep(NA, summ_row),'county fe'=rep(NA, summ_row),'treatment fe'=rep(NA, summ_row),
                             'weights'=rep(NA, summ_row), 'cox'=rep(NA, summ_row),
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
        post = ifelse(as.numeric(year) > years, 1, 0),
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
      mutate(year = as.numeric(year),
             ystar = ystar + p_err,
             ystar_cf = ystar_cf + p_err,
             y = (ystar > 0)*1 ,
             y_cf = (ystar_cf > 0)*1 )
    
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
    
    
    # aggregate up to county in each year 
    suppressWarnings(
      gridlevel_df <-  aggregate(panels, by = list(panels$year, panels$grid), FUN = mean, drop = TRUE)[c("grid", "treat", "post", "year","defor", "garea")]
    )
    
    gridlevel_df <- gridlevel_df[order(gridlevel_df$grid, gridlevel_df$year),]
    gridlevel_df <- slide(gridlevel_df, Var = "defor", GroupVar = "grid", NewVar = "deforlag",
                          slideBy = -1, reminder = FALSE)
    
    ##### creating forested share variable #####
    gridlevel_df$forshare <- 1 -  gridlevel_df$defor
    gridlevel_df$forsharelag <- 1 -  gridlevel_df$deforlag
    
    #generate outcome var
    gridlevel_df$deforrate <- ((gridlevel_df$forsharelag- gridlevel_df$forshare) / gridlevel_df$forsharelag)
    #remove any infinite values
    #gridlevel_df <- subset(gridlevel_df, select = -c(geometry))
    gridlevel_df <- 
      gridlevel_df %>% 
      filter_all(all_vars(!is.infinite(.)))
    
    # aggregate up to property in each year 
    suppressWarnings(
      proplevel_df <-  aggregate(panels, by = list(panels$property, panels$year), FUN = mean, drop = TRUE)[c("property", "treat", "post", "year","defor", "parea")]
    )
    
    
    proplevel_df <- proplevel_df[order(proplevel_df$property, proplevel_df$year),]
    proplevel_df <- slide(proplevel_df, Var = "defor", GroupVar = "property", NewVar = "deforlag",
                          slideBy = -1, reminder = FALSE)
    
    ##### creating forested share variable #####
    proplevel_df$forshare <- 1 -  proplevel_df$defor
    proplevel_df$forsharelag <- 1 -  proplevel_df$deforlag
    
    #generate outcome var
    proplevel_df$deforrate <- ((proplevel_df$forsharelag- proplevel_df$forshare) / proplevel_df$forsharelag)
    #remove any infinite values
    
    proplevel_df <- 
      proplevel_df %>% 
      filter_all(all_vars(!is.infinite(.)))
    
    # aggregate up to county in each year 
    suppressWarnings(
      countylevel_df <-  aggregate(panels, by = list(panels$county, panels$treat, panels$year), FUN = mean, drop = TRUE)[c("county", "property", "treat", "post", "year","defor", "carea")]
    )
    
    
    countylevel_df <- countylevel_df[order(countylevel_df$county, countylevel_df$year),]
    countylevel_df <- slide(countylevel_df, Var = "defor", GroupVar = "county", NewVar = "deforlag",
                            slideBy = -1, reminder = FALSE)
    
    ##### creating forested share variable #####
    countylevel_df$forshare <- 1 -  countylevel_df$defor
    countylevel_df$forsharelag <- 1 -  countylevel_df$deforlag
    
    #generate outcome var
    countylevel_df$deforrate <- ((countylevel_df$forsharelag- countylevel_df$forshare) / countylevel_df$forsharelag)
    #remove any infinite values
    
    countylevel_df <- 
      countylevel_df %>% 
      filter_all(all_vars(!is.infinite(.)))
    
    # aggregated units of analysis
    
    agg_DID1 <- feols(deforrate ~  post*treat|year+grid, data = gridlevel_df)
    
    
    agg_DID2 <- feols(deforrate ~  post*treat|year+property, data = proplevel_df)
    
    
    agg_DID3 <- feols(deforrate ~  post*treat|year+county, data = countylevel_df)
    
    # weighting by county and property area
    
    weight_DIDp <- feols(deforrate ~  post*treat|year+property, data = proplevel_df, weights = proplevel_df$parea)
    
    
    weight_DIDc <- feols(deforrate ~  post*treat|year+county, , weights = countylevel_df$carea, data = countylevel_df)
    
    # problematic specifications
    
    bad_DID1 <- feols(defor ~  post*treat, data = panels)
    
    bad_DID2 <- feols(y_it ~  post*treat|year+pixels, data = panels)
    
    
    
    ### TWFE regressions with aggregated fixed effects
    
    fix_DID1 <- feols(y_it ~  post*treat|year+grid, data = panels)
    
    fix_DID2 <- feols(y_it ~  post*treat|year+property, data = panels)
    
    fix_DID3 <- feols(y_it ~  post*treat|year+county, data = panels)
    
    # simple DID
    
    DID <- feols(y_it ~  post*treat, data = panels)
    
    
    #calculating bias from each method
    
    # aggregated units of analysis
    coeffmatrix[i,1] <- agg_DID1$coefficients - ATT
    coeffmatrix[i,2] <- agg_DID2$coefficients - ATT
    coeffmatrix[i,3] <- agg_DID3$coefficients - ATT
    
    # with weights
    weight_coeffmatrix[i,1] <- weight_DIDp$coefficients - ATT
    weight_coeffmatrix[i,2] <- weight_DIDc$coefficients - ATT
    
    # problematic specifications
    bad_coeffmatrix[i,1] <- bad_DID1$coefficients[4] - ATT
    bad_coeffmatrix[i,2] <- bad_DID2$coefficients - ATT
    
    # aggregated fixed effects
    fix_coeffmatrix[i,1] <- tail(fix_DID1$coefficients, n=1) - ATT
    fix_coeffmatrix[i,2] <- tail(fix_DID2$coefficients, n=1) - ATT
    fix_coeffmatrix[i,3] <- tail(fix_DID3$coefficients, n=1) - ATT
    
    # regular did
    fix_coeffmatrix[i,4] <- DID$coefficients[4] - ATT
    
    #regular DID standard errors
    # robust
    DID_se <- summary(DID, se = "hetero")$se[4]
    # clustered at pixel
    DID_clse_pixel  <- tail(summary(DID, cluster = ~pixels)$se, n=1)
    # clustered pixel level standard errors at other units of aggregation
    DID_clse_grid  <- tail(summary(DID, cluster = ~grid)$se, n=1)
    DID_clse_property <- tail(summary(DID, cluster = ~property)$se, n=1)
    DID_clse_county <- tail(summary(DID, cluster = ~county)$se, n=1)
    
    did_covermat[i,1] <- dplyr::between(ATT, DID$coefficients[4] - 1.96 * DID_se, DID$coefficients[4] + 1.96 * DID_se)*1
    did_covermat[i,2] <- dplyr::between(ATT, DID$coefficients[4] - 1.96 * DID_clse_pixel, DID$coefficients[4] + 1.96 * DID_clse_pixel)*1
    did_covermat[i,3] <- dplyr::between(ATT, DID$coefficients[4] - 1.96 * DID_clse_grid, DID$coefficients[4] + 1.96 * DID_clse_grid)*1
    did_covermat[i,4] <- dplyr::between(ATT, DID$coefficients[4] - 1.96 * DID_clse_property, DID$coefficients[4] + 1.96 * DID_clse_property)*1
    did_covermat[i,5] <- dplyr::between(ATT, DID$coefficients[4] - 1.96 * DID_clse_county, DID$coefficients[4] + 1.96 * DID_clse_county)*1
    
    
    #clustering at group level for aggregated analyses
    agg_clse1    <- tail(summary(agg_DID1, cluster = ~grid)$se, n=1)
    agg_clse2    <- tail(summary(agg_DID2, cluster = ~property)$se, n=1)
    agg_clse3    <- tail(summary(agg_DID3, cluster = ~county)$se, n=1)
    
    # same for weighted
    weight_clsep    <- tail(summary(weight_DIDp, cluster = ~property)$se, n=1)
    weight_clsec   <- tail(summary(weight_DIDc, cluster = ~county)$se, n=1)
    
    # coverage with clustered standard errors at group level
    agg_covermat[i,1] <- dplyr::between(ATT, agg_DID1$coefficients - 1.96 * agg_clse1, agg_DID1$coefficients + 1.96 * agg_clse1)*1
    agg_covermat[i,2] <- dplyr::between(ATT, agg_DID2$coefficients - 1.96 * agg_clse2, agg_DID2$coefficients + 1.96 * agg_clse2)*1
    agg_covermat[i,3] <- dplyr::between(ATT, agg_DID3$coefficients - 1.96 * agg_clse3, agg_DID3$coefficients + 1.96 * agg_clse3)*1
    
    #weighted
    weight_covermat[i,1] <- dplyr::between(ATT, weight_DIDp$coefficients - 1.96 * weight_clsep, weight_DIDp$coefficients + 1.96 * weight_clsep)*1
    weight_covermat[i,2] <- dplyr::between(ATT, weight_DIDc$coefficients - 1.96 * weight_clsec, weight_DIDc$coefficients + 1.96 * weight_clsec)*1
    
    
    ### clustered standard errors for aggregated fixed effects
    fix_clse1    <- tail(summary(fix_DID1, cluster = ~grid)$se, n=1)
    fix_clse2    <- tail(summary(fix_DID2, cluster = ~property)$se, n=1)
    fix_clse3    <- tail(summary(fix_DID3, cluster = ~county)$se, n=1)
    
    #whether att is within CI
    
    
    fix_covermat[i,1] <- dplyr::between(ATT, tail(fix_DID1$coefficients, n=1) - 1.96 * fix_clse1, tail(fix_DID1$coefficients, n=1) + 1.96 * fix_clse1)*1
    fix_covermat[i,2] <- dplyr::between(ATT, tail(fix_DID2$coefficients, n=1) - 1.96 * fix_clse2, tail(fix_DID2$coefficients, n=1) + 1.96 * fix_clse2)*1
    fix_covermat[i,3] <- dplyr::between(ATT, tail(fix_DID3$coefficients, n=1) - 1.96 * fix_clse3, tail(fix_DID3$coefficients, n=1) + 1.96 * fix_clse3)*1
    
    # se for bad specifications clustered at pixel level
    bad_clse1    <- tail(summary(bad_DID1, cluster = ~pixels)$se, n=1)
    bad_clse2    <- tail(summary(bad_DID2, cluster = ~pixels)$se, n=1)
    
    #coverage for bad specifications
    bad_covermat[i,1] <- dplyr::between(ATT, tail(bad_DID1$coefficients, n=1) - 1.96 * bad_clse1, tail(bad_DID1$coefficients, n=1) + 1.96 * bad_clse1)*1
    bad_covermat[i,2] <- dplyr::between(ATT, bad_DID2$coefficients - 1.96 * bad_clse2, bad_DID2$coefficients + 1.96 * bad_clse2)*1
    
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ### survival analysis
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    haz_rat <- pnorm(b0+b1+b2_1+b3, 0, (std_a^2+std_v^2 +std_p^2)^.5) / pnorm(b0+b1+b2_1, 0, (std_a^2+std_v^2 + std_p^2)^.5)
    
    surv_df <- panels %>% 
      mutate(t_start = year - 1,
             t_end = year,
             # t_end = ifelse((t_end==20) & (y_it==0), Inf, t_end),
             outcome = y_it,
             treat_now = treat * post) %>% 
      select(pixels, t_start, t_end, outcome, treat, treat_now, post, year) %>% 
      drop_na()%>%
      mutate(t_start = ifelse(post==1, t_start -10, t_start),
             t_end = ifelse(post==1, t_end -10, t_end))
    
    
    ## this model is more similar to the traditional DID setup but does not actually recover the desired HRR
    cox_did <- coxph(Surv(t_start, t_end, outcome) ~ post*treat
                             , data = surv_df )
    #summary(cox_interaction)
    coefs <- cox_did$coefficients
    hr_did <- coefs[[3]] %>% exp()
    
    vcoeff <- cox_did$var[3,3]
    se_cox_did <- deltamethod(~ exp(x1), coefs[[3]], vcoeff)
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ### now we'll try to estimate separately the components of the desired hazard ratio
    
    cox_post <- coxph(Surv(t_start, t_end, outcome) ~ treat
                      , data = surv_df %>% filter(post==1))
    #summary(cox_post)
    hr_11_01 <- cox_post$coefficients %>% exp()
    
    cox_treat <- coxph(Surv(t_start, t_end, outcome) ~ post
                       , data = surv_df %>% filter(treat==1))
    #summary(cox_treat)
    hr_11_10 <- cox_treat$coefficients %>% exp()
    
    cox_notreat <- coxph(Surv(t_start, t_end, outcome) ~ post
                         , data = surv_df %>% filter(treat==0))
    #summary(cox_notreat)
    hr_01_00 <- cox_notreat$coefficients %>% exp()
    
    hr_11_cf <- 1/(1/hr_11_10 + 1/hr_11_01 - (1/(hr_11_01*hr_01_00)))
    
    vcov_hr <- matrix(0, nrow = 3, ncol = 3)
    vcov_hr[1,1] <- cox_treat$var
    vcov_hr[2,2] <- cox_post$var
    vcov_hr[3,3] <- cox_notreat$var
    
    se_hr_11_cf <- deltamethod(
      ~ 1/(1/exp(x1) + 1/exp(x2) - (1/(exp(x2)*exp(x3)))),
      c(cox_treat$coefficients, cox_post$coefficients, cox_notreat$coefficients),
      vcov_hr
    )
    
    #### calculating observed deforestation rate to transition to ATT estimate
    defor_summary <- panels %>%
      group_by(treat, post) %>%
      summarise(mean_y_it = mean(y_it, na.rm = TRUE))%>%
      ungroup()
    
    d_obs = defor_summary[4,3]$mean_y_it
    
    ATT_11_cf <- d_obs - d_obs/hr_11_cf
    ATT_cox <- d_obs - d_obs/hr_did
    
    hr_11_cf_cover <- dplyr::between(haz_rat, hr_11_cf - 1.96 * se_hr_11_cf, hr_11_cf + 1.96 * se_hr_11_cf)*1
    cox_cover <- dplyr::between(haz_rat, hr_did - 1.96 * se_cox_did, hr_did + 1.96 * se_cox_did)*1
    
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ### filling in summary long dataframe output
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    firstcol = which(colnames(summary_long)=="iteration")
    lastcol = which(colnames(summary_long)=="notes")
    
    summary_long[i,c(firstcol:lastcol)] <- c(
      i,
      1,0,0,0,
      0,0,0,0,1,
      0,0,
      1,0,0,0,
      bad_coeffmatrix[i,1],
      bad_covermat[i,1],
      "keeping pixels after deforestation event"
    )
    
    firstcol = which(colnames(summary_long)=="iteration")
    lastcol = which(colnames(summary_long)=="cover")
    
    summary_long[i+n,c(firstcol:lastcol)] <- c(
      i,
      1,0,0,0,
      1,0,0,0,0,
      0,0,
      1,0,0,0,
      bad_coeffmatrix[i,2],
      bad_covermat[i,2]
    )
    
    summary_long[i+n*2,c(firstcol:lastcol)] <- c(
      i,
      1,0,0,0,
      0,0,0,0,1,
      0,0,
      1,0,0,0,
      fix_coeffmatrix[i,4],
      did_covermat[i,2]
    )
    
    summary_long[i+n*3,c(firstcol:lastcol)] <- c(
      i,
      1,0,0,0,
      0,0,0,0,1,
      0,0,
      0,1,0,0,
      fix_coeffmatrix[i,4],
      did_covermat[i,3]
    )
    
    summary_long[i+n*4,c(firstcol:lastcol)] <- c(
      i,
      1,0,0,0,
      0,0,0,0,1,
      0,0,
      0,0,1,0,
      fix_coeffmatrix[i,4],
      did_covermat[i,4]
    )
    
    summary_long[i+n*5,c(firstcol:lastcol)] <- c(
      i,
      1,0,0,0,
      0,0,0,0,1,
      0,0,
      0,0,0,1,
      fix_coeffmatrix[i,4],
      did_covermat[i,5]
    )
    
    summary_long[i+n*6,c(firstcol:lastcol)] <- c(
      i,
      0,1,0,0,
      0,1,0,0,0,
      0,0,
      0,1,0,0,
      coeffmatrix[i,1],
      agg_covermat[i,1]
    )
    
    summary_long[i+n*7,c(firstcol:lastcol)] <- c(
      i,
      0,0,1,0,
      0,0,1,0,0,
      0,0,
      0,0,1,0,
      coeffmatrix[i,2],
      agg_covermat[i,2]
    )
    
    summary_long[i+n*8,c(firstcol:lastcol)] <- c(
      i,
      0,0,0,1,
      0,0,0,1,0,
      0,0,
      0,0,0,1,
      coeffmatrix[i,3],
      agg_covermat[i,3]
    )
    
    summary_long[i+n*9,c(firstcol:lastcol)] <- c(
      i,
      1,0,0,0,
      0,1,0,0,0,
      0,0,
      0,1,0,0,
      fix_coeffmatrix[i,1],
      fix_covermat[i,1]
    )
    
    summary_long[i+n*10,c(firstcol:lastcol)] <- c(
      i,
      1,0,0,0,
      0,0,1,0,0,
      0,0,
      0,0,1,0,
      fix_coeffmatrix[i,2],
      fix_covermat[i,2]
    )
    
    summary_long[i+n*11,c(firstcol:lastcol)] <- c(
      i,
      1,0,0,0,
      0,0,0,1,0,
      0,0,
      0,0,0,1,
      fix_coeffmatrix[i,3],
      fix_covermat[i,3]
    )
    
    summary_long[i+n*12,c(firstcol:lastcol)] <- c(
      i,
      0,0,1,0,
      0,0,1,0,0,
      1,0,
      0,0,1,0,
      weight_coeffmatrix[i,1],
      weight_covermat[i,1]
    )
    
    summary_long[i+n*13,c(firstcol:lastcol)] <- c(
      i,
      0,0,0,1,
      0,0,0,1,0,
      1,0,
      0,0,0,1,
      weight_coeffmatrix[i,2],
      weight_covermat[i,2]
    )
    
    firstcol = which(colnames(summary_long)=="iteration")
    lastcol = which(colnames(summary_long)=="notes")
    
    summary_long[i+n*14,c(firstcol:lastcol)] <- c(
      i,
      1,0,0,0,
      0,0,0,0,1,
      0,1,
      1,0,0,0,
      ATT_cox- ATT,
      cox_cover,
      "Cox DID regression"
    )
    
    summary_long[i+n*15,c(firstcol:lastcol)] <- c(
      i,
      1,0,0,0,
      0,0,0,0,0,
      0,1,
      1,0,0,0,
      ATT_11_cf- ATT,
      hr_11_cf_cover,
      "Cox HRATT estimator"
    )
    
    
    print(i)
    toc()
  }
  
  outputs = list(
    "summary_long" = summary_long)
  
  return(outputs)
  
  #end function  
}  