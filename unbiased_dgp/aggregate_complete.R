#this function uses the property_scape gen fcn to generate a landscape with properties. Then we introduce property level perturbations and compare estimates when the data is aggregateed to the grid vs. property level


library(ggplot2)
library(clubSandwich)
library(reshape2)
library(matrixStats)
library(ggplot2)
library(plm)
library(Metrics)
library(DataCombine)
library(dplyr)
library(tidyverse)
library(tictoc)
library(fixest)
library(DeclareDesign)
source('full_landscape.R')

#begin function
aggregate_complete <- function(n, nobs, years, b0, b1, b2_0, b2_1, b3, std_a = 0.1, std_v = 0.25, std_p = 0.0, cellsize, ppoints, cpoints){
  
  countyscape = full_landscape(nobs, cellsize, ppoints, cpoints)
  pixloc_df = countyscape$pixloc_df
  
  ATT <- pnorm(b0+b1+b2_1+b3, 0, (std_a^2+std_v^2 +std_p^2)^.5) - pnorm(b0+b1+b2_1, 0, (std_a^2+std_v^2 + std_p^2)^.5)
  
  pixloc <- pixloc_df
  
  did_covermat <- matrix(nrow = n, ncol = 5)
  agg_covermat <- matrix(nrow = n, ncol = 3)
  bad_covermat <- matrix(nrow = n, ncol = 2)
  fix_covermat <- matrix(nrow = n, ncol = 3)
  
  coeffmatrix <- matrix(nrow = n, ncol = 4)
  bad_coeffmatrix <- matrix(nrow = n, ncol = 2)
  fix_coeffmatrix <- matrix(nrow = n, ncol = 3)
  
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
        ystar = b0 + b1*treat + b2_0*post*(1-treat) + b2_1*post*treat + b3*treat*post + a_i + v_it
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
      mutate(ystar = ystar + p_err) %>%
      mutate(y = (ystar > 0)*1 )
    
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
    
    DID <- feols(y_it ~  post*treat, data = panels)
    
    # problematic specifications
    
    bad_DID1 <- feols(defor ~  post*treat, data = panels)
    
    bad_DID2 <- feols(y_it ~  post*treat|year+pixels, data = panels)
    
    
    
    ### TWFE regressions with aggregated fixed effects
    
    fix_DID1 <- feols(y_it ~  post*treat|year+grid, data = panels)
    
    fix_DID2 <- feols(y_it ~  post*treat|year+property, data = panels)
    
    fix_DID3 <- feols(y_it ~  post*treat|year+county, data = panels)
    
    
    #calculating bias from each method
    
    # aggregated units of analysis
    coeffmatrix[i,1] <- agg_DID1$coefficients - ATT
    coeffmatrix[i,2] <- agg_DID2$coefficients - ATT
    coeffmatrix[i,3] <- agg_DID3$coefficients - ATT
    
    # regular did
    coeffmatrix[i,4] <- DID$coefficients[4] - ATT
    
    # problematic specifications
    bad_coeffmatrix[i,1] <- bad_DID1$coefficients[4] - ATT
    bad_coeffmatrix[i,2] <- bad_DID2$coefficients - ATT
    
    # aggregated fixed effects
    fix_coeffmatrix[i,1] <- tail(fix_DID1$coefficients, n=1) - ATT
    fix_coeffmatrix[i,2] <- tail(fix_DID2$coefficients, n=1) - ATT
    fix_coeffmatrix[i,3] <- tail(fix_DID3$coefficients, n=1) - ATT
    
    #regular DID standard errors
    # robust
    DID_se <- summary(DID, se = "hetero")$se[4]
    # clustered at pixel
    DID_clse_pixel  <- tail(summary(DID, cluster = ~pixels)$se, n=1)
    # clustered pixel level standard errors at other units of aggregation
    DID_clse_grid  <- tail(summary(DID, cluster = ~grid)$se, n=1)
    DID_clse_property <- tail(summary(DID, cluster = ~property)$se, n=1)
    DID_clse_county <- tail(summary(DID, cluster = ~county)$se, n=1)
    
    did_covermat[i,1] <- between(ATT, DID$coefficients[4] - 1.96 * DID_se, DID$coefficients[4] + 1.96 * DID_se)*1
    did_covermat[i,2] <- between(ATT, DID$coefficients[4] - 1.96 * DID_clse_pixel, DID$coefficients[4] + 1.96 * DID_clse_pixel)*1
    did_covermat[i,3] <- between(ATT, DID$coefficients[4] - 1.96 * DID_clse_grid, DID$coefficients[4] + 1.96 * DID_clse_grid)*1
    did_covermat[i,4] <- between(ATT, DID$coefficients[4] - 1.96 * DID_clse_property, DID$coefficients[4] + 1.96 * DID_clse_property)*1
    did_covermat[i,5] <- between(ATT, DID$coefficients[4] - 1.96 * DID_clse_county, DID$coefficients[4] + 1.96 * DID_clse_county)*1
    
    
    #clustering at group level for aggregated analyses
    agg_clse1    <- tail(summary(agg_DID1, cluster = ~grid)$se, n=1)
    agg_clse2    <- tail(summary(agg_DID2, cluster = ~property)$se, n=1)
    agg_clse3    <- tail(summary(agg_DID3, cluster = ~county)$se, n=1)
    
    # coverage with clustered standard errors at group level
    agg_covermat[i,1] <- between(ATT, agg_DID1$coefficients - 1.96 * agg_clse1, agg_DID1$coefficients + 1.96 * agg_clse1)*1
    agg_covermat[i,2] <- between(ATT, agg_DID2$coefficients - 1.96 * agg_clse2, agg_DID2$coefficients + 1.96 * agg_clse2)*1
    agg_covermat[i,3] <- between(ATT, agg_DID3$coefficients - 1.96 * agg_clse3, agg_DID3$coefficients + 1.96 * agg_clse3)*1
    
    
    ### clustered standard errors for aggregated fixed effects
    fix_clse1    <- tail(summary(fix_DID1, cluster = ~grid)$se, n=1)
    fix_clse2    <- tail(summary(fix_DID2, cluster = ~property)$se, n=1)
    fix_clse3    <- tail(summary(fix_DID3, cluster = ~county)$se, n=1)
    
    #whether att is within CI
    ### coverage matrix for aggregated fixed effects
    
    fix_covermat[i,1] <- between(ATT, tail(fix_DID1$coefficients, n=1) - 1.96 * fix_clse1, tail(fix_DID1$coefficients, n=1) + 1.96 * fix_clse1)*1
    fix_covermat[i,2] <- between(ATT, tail(fix_DID2$coefficients, n=1) - 1.96 * fix_clse2, tail(fix_DID2$coefficients, n=1) + 1.96 * fix_clse2)*1
    fix_covermat[i,3] <- between(ATT, tail(fix_DID3$coefficients, n=1) - 1.96 * fix_clse3, tail(fix_DID3$coefficients, n=1) + 1.96 * fix_clse3)*1
    
    # se for bad specifications clustered at pixel level
    bad_clse1    <- tail(summary(bad_DID1, cluster = ~pixels)$se, n=1)
    bad_clse2    <- tail(summary(bad_DID2, cluster = ~pixels)$se, n=1)
    
    #coverage for bad specifications
    bad_covermat[i,1] <- between(ATT, tail(bad_DID1$coefficients, n=1) - 1.96 * bad_clse1, tail(bad_DID1$coefficients, n=1) + 1.96 * bad_clse1)*1
    bad_covermat[i,2] <- between(ATT, bad_DID2$coefficients - 1.96 * bad_clse2, bad_DID2$coefficients + 1.96 * bad_clse2)*1
    
    
    print(i)
    toc()
  }
  
  coeff_bias <- as.data.frame(coeffmatrix)
  bad_bias <- as.data.frame(bad_coeffmatrix)
  actual <- rep(0, times = n)
  names(coeff_bias)[1] <- paste("grid")
  names(coeff_bias)[2] <- paste("property")
  names(coeff_bias)[3] <- paste("county")
  names(coeff_bias)[4] <- paste("pixel")
  suppressWarnings(cbias <- melt(coeff_bias, value.name = "bias"))
  suppressWarnings(bad_cbias <- melt(bad_bias, value.name = "bias"))
  
  fix_bias <- as.data.frame(fix_coeffmatrix)
  names(fix_bias)[1] <- paste("grid fe")
  names(fix_bias)[2] <- paste("property fe")
  names(fix_bias)[3] <- paste("county fe")
  suppressWarnings(fix_cbias <- melt(fix_bias, value.name = "bias"))
  
  ##################################################################
  ###### Bias distribution Plots #######
  
  plot <- ggplot(data = cbias, aes(x = bias, fill=variable)) +
    geom_density(alpha = .2) +
    guides(fill=guide_legend(title=NULL))+
    scale_fill_discrete(breaks=c("grid", "property", "county", "pixel"), labels=c("aggregated to grids", "aggregated to properties", "aggregated to counties", "pixel"))+
    geom_vline(xintercept = 0, linetype = "dashed")+
    #geom_vline(aes(xintercept= (DID_estimand - ATT), color="DID estimand - ATT"), linetype="dashed")+
    #theme(plot.margin = unit(c(1,1,3,1), "cm"))+
    theme(plot.caption = element_text(hjust = 0.5))+
    labs(x= "Bias", caption = paste("Mean grid:", round(mean(coeff_bias$grid), digits = 4),
                                    ", RMSE:", round(rmse(actual, coeff_bias$grid), digits = 4), "\n", 
                                    "Mean property:", round(mean(coeff_bias$property), digits = 4),
                                    ", RMSE:", round(rmse(actual, coeff_bias$property), digits = 4),"\n",
                                    "Mean county:", round(mean(coeff_bias$county), digits = 4),
                                    ", RMSE:", round(rmse(actual, coeff_bias$county), digits = 4), "\n",
                                    "Mean pixel:", round(mean(coeff_bias$pixel), digits = 4),
                                    ", RMSE:", round(rmse(actual, coeff_bias$pixel), digits = 4)) 
    )
  
  fe_plot <- ggplot(data = fix_cbias, aes(x = bias, fill=variable)) +
    geom_density(alpha = .2) +
    guides(fill=guide_legend(title=NULL))+
    #scale_fill_discrete(breaks=c("grid fe", "property fe", "county fe"), labels=c("grid fe", "property fe", "county fe"))+
    geom_vline(xintercept = 0, linetype = "dashed")+
    #geom_vline(aes(xintercept= (DID_estimand - ATT), color="DID estimand - ATT"), linetype="dashed")+
    #theme(plot.margin = unit(c(1,1,3,1), "cm"))+
    theme(plot.caption = element_text(hjust = 0.5))+
    labs(x= "Bias", caption = paste("Mean grid fe:", round(mean(fix_bias[,1]), digits = 4),
                                    ", RMSE:", round(rmse(actual, fix_bias[,1]), digits = 4), "\n", 
                                    "Mean property fe:", round(mean(fix_bias[,2]), digits = 4),
                                    ", RMSE:", round(rmse(actual, fix_bias[,2]), digits = 4),"\n",
                                    "Mean county fe:", round(mean(fix_bias[,3]), digits = 4),
                                    ", RMSE:", round(rmse(actual, fix_bias[,3]), digits = 4))
    )
  
  
  did_coverages_df <- data.frame( 
    aggregation = c('pixel',
                    'pixel',
                    'pixel',
                    'pixel',
                    'pixel'),
    std_error = c('robust',
                  'clustered at pixel',
                  'clustered at grid',
                  'clustered at property'
                  , 'clustered at county'),
    coverage = colMeans(did_covermat)
  )
  
  agg_coverages_df <- data.frame( 
    aggregation = c('grid', 
                    'property',
                    'county'),
    std_error = c('clustered at grid', 
                  'clustered at property',
                  'clustered at county'),
    coverage = colMeans(agg_covermat)
  )
  
  fix_coverages_df <- data.frame( 
    fix_effects = c('grid', 
                    'property',
                    'county',
                    'grid',
                    'property',
                    'county'),
    std_error = c('clustered at grid', 
                  'clustered at property',
                  'clustered at county'),
    coverage = colMeans(fix_covermat)
  )
  
  summary_df <- data.frame('pixel'=rep(NA, 12),'grid'=rep(NA, 12),'property'=rep(NA, 12),'county'=rep(NA, 12),
                           'pixel fe'=rep(NA, 12),'grid fe'=rep(NA, 12),'property fe'=rep(NA, 12),'county fe'=rep(NA, 12),'treatment fe'=rep(NA, 12),
                           'se_pixel'=rep(NA, 12), 'se_grid'=rep(NA, 12), 'se_property'=rep(NA, 12), 'se_county'=rep(NA, 12),
                           'mean_bias'=rep(NA, 12), 'RMSE'=rep(NA, 12),
                           'q05'=rep(NA, 12), 'q95'=rep(NA, 12),
                           'cover'=rep(NA, 12),
                           stringsAsFactors=FALSE)
  
  # DID1 <- feols(deforrate ~  post*treat|year+grid, data = gridlevel_df)
  summary_df[7,] <- c(
                      0,1,0,0,
                      0,1,0,0,0,
                      0,1,0,0,
                      mean(coeff_bias$grid), rmse(actual, coeff_bias$grid),
                      quantile(coeff_bias$grid, 0.05), quantile(coeff_bias$grid, 0.95),
                      #quantile(coeff_bias$grid, 0.01), quantile(coeff_bias$grid, 0.99),
                      mean(agg_covermat[,1]))
  
  # DID2 <- feols(deforrate ~  post*treat|year+property, data = proplevel_df)
  summary_df[8,] <- c(
                      0,0,1,0,
                      0,0,1,0,0,
                      0,0,1,0,
                      mean(coeff_bias$property), rmse(actual, coeff_bias$property),
                      quantile(coeff_bias$property, 0.05), quantile(coeff_bias$property, 0.95),
                      #quantile(coeff_bias$property, 0.01), quantile(coeff_bias$property, 0.99),
                      mean(agg_covermat[,2]))
  
  # DID3 <- feols(deforrate ~  post*treat|year+county, data = countylevel_df)
  summary_df[9,] <- c(
                      0,0,0,1,
                      0,0,0,1,0,
                      0,0,0,1,
                      mean(coeff_bias$county), rmse(actual, coeff_bias$county),
                      quantile(coeff_bias$county, 0.05), quantile(coeff_bias$county, 0.95),
                      #quantile(coeff_bias$county, 0.01), quantile(coeff_bias$county, 0.99),
                      mean(agg_covermat[,3]))
  
  
  # DID4 <- feols(y_it ~  post*treat, data = panels)
  summary_df[3,] <- c(
                      1,0,0,0,
                      0,0,0,0,1,
                      1,0,0,0,
                      mean(coeff_bias$pixel), rmse(actual, coeff_bias$pixel),
                      quantile(coeff_bias$pixel, 0.05), quantile(coeff_bias$pixel, 0.95),
                      #quantile(coeff_bias$pixel, 0.01), quantile(coeff_bias$pixel, 0.99),
                      mean(did_covermat[,2]))
  
  summary_df[4,] <- c(
    1,0,0,0,
    0,0,0,0,1,
    0,1,0,0,
    mean(coeff_bias$pixel), rmse(actual, coeff_bias$pixel),
    quantile(coeff_bias$pixel, 0.05), quantile(coeff_bias$pixel, 0.95),
    #quantile(coeff_bias$pixel, 0.01), quantile(coeff_bias$pixel, 0.99),
    mean(did_covermat[,3]))
  
  summary_df[5,] <- c(
    1,0,0,0,
    0,0,0,0,1,
    0,0,1,0,
    mean(coeff_bias$pixel), rmse(actual, coeff_bias$pixel),
    quantile(coeff_bias$pixel, 0.05), quantile(coeff_bias$pixel, 0.95),
    #quantile(coeff_bias$pixel, 0.01), quantile(coeff_bias$pixel, 0.99),
    mean(did_covermat[,4]))
  
  summary_df[6,] <- c(
    1,0,0,0,
    0,0,0,0,1,
    0,0,0,1,
    mean(coeff_bias$pixel), rmse(actual, coeff_bias$pixel),
    quantile(coeff_bias$pixel, 0.05), quantile(coeff_bias$pixel, 0.95),
    #quantile(coeff_bias$pixel, 0.01), quantile(coeff_bias$pixel, 0.99),
    mean(did_covermat[,5]))
  
  
  summary_df[1,] <- c(
                      1,0,0,0,
                      0,0,0,0,1,
                      1,0,0,0,
                      mean(bad_bias[,1]), rmse(actual, bad_bias[,1]),
                      quantile(bad_bias[,1], 0.05), quantile(bad_bias[,1], 0.95),
                      #quantile(bad_bias[,1], 0.01), quantile(bad_bias[,1], 0.99),
                      mean(bad_covermat[,1]))
  
  # DID6 <- feols(y_it ~  post*treat|year+pixels, data = panels)
  summary_df[2,] <- c(
                      1,0,0,0,
                      1,0,0,0,0,
                      1,0,0,0,
                      mean(bad_bias[,2]), rmse(actual, bad_bias[,2]),
                      quantile(bad_bias[,2], 0.05), quantile(bad_bias[,2], 0.95),
                      #quantile(bad_bias[,2], 0.01), quantile(bad_bias[,2], 0.99),
                      mean(bad_covermat[,2]))
  
  ### TWFE regressions with aggregated fixed effects
  
  # fix_DID1 <- feols(y_it ~  post*treat|year+grid, data = panels)
  summary_df[10,] <- c(
                      1,0,0,0,
                      0,1,0,0,0,
                      0,1,0,0,
                      mean(fix_bias[,1]), rmse(actual, fix_bias[,1]),
                      quantile(fix_bias[,1], 0.05), quantile(fix_bias[,1], 0.95),
                      #quantile(fix_bias[,1], 0.01), quantile(fix_bias[,1], 0.99),
                      mean(fix_covermat[,1]))
  
  # fix_DID2 <- feols(y_it ~  post*treat|year+property, data = panels)
  summary_df[11,] <- c(
                      1,0,0,0,
                      0,0,1,0,0,
                      0,0,1,0,
                      mean(fix_bias[,2]), rmse(actual, fix_bias[,2]),
                      quantile(fix_bias[,2], 0.05), quantile(fix_bias[,2], 0.95),
                      #quantile(fix_bias[,2], 0.01), quantile(fix_bias[,2], 0.99),
                      mean(fix_covermat[,2]))
  
  # fix_DID3 <- feols(y_it ~  post*treat|year+county, data = panels)
  summary_df[12,] <- c(
                      1,0,0,0,
                      0,0,0,1,0,
                      0,0,0,1,
                      mean(fix_bias[,3]), rmse(actual, fix_bias[,3]),
                      quantile(fix_bias[,3], 0.05), quantile(fix_bias[,3], 0.95),
                      #quantile(fix_bias[,3], 0.01), quantile(fix_bias[,3], 0.99),
                      mean(fix_covermat[,3]))
  
  summary_df <- summary_df %>%
    mutate_at(1:13, as.logical)%>%
    dplyr::select(mean_bias, everything())
  
  
  outputs = list("plot" = plot, "fe_plot" = fe_plot, "biases" = coeff_bias, "fe_biases" = fix_bias, "did_coverages_df" = did_coverages_df, "agg_coverages_df" = agg_coverages_df, "fe_coverages_df" = fix_coverages_df, "summary_df" = summary_df)
  return(outputs)
  
  #end function  
}  
