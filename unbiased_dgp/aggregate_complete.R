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
  
  covermat <- matrix(nrow = n, ncol = 10)
  fix_covermat <- matrix(nrow = n, ncol = 6)
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

    
    DID1 <- feols(deforrate ~  post*treat|year+grid, data = gridlevel_df)
    
    
    DID2 <- feols(deforrate ~  post*treat|year+property, data = proplevel_df)

    
    DID3 <- feols(deforrate ~  post*treat|year+county, data = countylevel_df)
    
    
    DID4 <- feols(y_it ~  post*treat, data = panels)
    
    DID5 <- feols(defor ~  post*treat, data = panels)
    
    DID6 <- feols(y_it ~  post*treat|year+pixels, data = panels)
    
    
    
    ### TWFE regressions with aggregated fixed effects
    
    fix_DID1 <- feols(y_it ~  post*treat|year+grid, data = panels)
    
    fix_DID2 <- feols(y_it ~  post*treat|year+property, data = panels)
    
    fix_DID3 <- feols(y_it ~  post*treat|year+county, data = panels)
    
    
    #calculating bias from each aggregation method
    coeffmatrix[i,1] <- DID1$coefficients - ATT
    coeffmatrix[i,2] <- DID2$coefficients - ATT
    coeffmatrix[i,3] <- DID3$coefficients - ATT
    coeffmatrix[i,4] <- DID4$coefficients[4] - ATT
    
    bad_coeffmatrix[i,1] <- DID5$coefficients[4] - ATT
    bad_coeffmatrix[i,2] <- DID6$coefficients - ATT
    
    fix_coeffmatrix[i,1] <- tail(fix_DID1$coefficients, n=1) - ATT
    fix_coeffmatrix[i,2] <- tail(fix_DID2$coefficients, n=1) - ATT
    fix_coeffmatrix[i,3] <- tail(fix_DID3$coefficients, n=1) - ATT
    
    # classical standard errors
    se1 <- summary(DID1, se = "hetero")$se[1]
    se2 <- summary(DID2, se = "hetero")$se[1]
    se3 <- summary(DID3, se = "hetero")$se[1]
    se4 <- summary(DID4, se = "hetero")$se[4]
    
    #clustering at group level for aggregated analyses
    cluster_se1    <- tail(summary(DID1, cluster = ~grid)$se, n=1)
    cluster_se2    <- tail(summary(DID2, cluster = ~property)$se, n=1)
    cluster_se3    <- tail(summary(DID3, cluster = ~county)$se, n=1)
    
    # clustered pixel level standard errors
    cluster_prop_pix  <- tail(summary(DID4, cluster = ~property)$se, n=1)
    cluster_grid_pix <- tail(summary(DID4, cluster = ~grid)$se, n=1)
    cluster_county_pix <- tail(summary(DID4, cluster = ~county)$se, n=1)
    
    
    ### typical and clustered standard errors for aggregated fixed effects
    
    fix_se1 <- tail(summary(fix_DID1, se = "hetero")$se, n=1)
    fix_se2 <- tail(summary(fix_DID2, se = "hetero")$se, n=1)
    fix_se3 <- tail(summary(fix_DID3, se = "hetero")$se, n=1)
    
    fix_clse1    <- tail(summary(fix_DID1, cluster = ~grid)$se, n=1)
    fix_clse2    <- tail(summary(fix_DID2, cluster = ~property)$se, n=1)
    fix_clse3    <- tail(summary(fix_DID3, cluster = ~county)$se, n=1)
    
    #whether att is within CI
    
    # coverage with clustered standard errors at group level
    covermat[i,1] <- between(ATT, DID1$coefficients - 1.96 * cluster_se1, DID1$coefficients + 1.96 * cluster_se1)*1
    covermat[i,2] <- between(ATT, DID2$coefficients - 1.96 * cluster_se2, DID2$coefficients + 1.96 * cluster_se2)*1
    covermat[i,3] <- between(ATT, DID3$coefficients - 1.96 * cluster_se3, DID3$coefficients + 1.96 * cluster_se3)*1
    
    # classical standard errors
    covermat[i,4] <- between(ATT, DID1$coefficients - 1.96 * se1, DID1$coefficients + 1.96 * se1)*1
    covermat[i,5] <- between(ATT, DID2$coefficients - 1.96 * se2, DID2$coefficients + 1.96 * se2)*1
    covermat[i,6] <- between(ATT, DID3$coefficients - 1.96 * se3, DID3$coefficients + 1.96 * se3)*1
    covermat[i,7] <- between(ATT, DID4$coefficients[4] - 1.96 * se4, DID4$coefficients[4] + 1.96 * se4)*1
    
    # coverage with pixel level analyses clustered at different levels
    covermat[i,8] <- between(ATT, DID4$coefficients[4] - 1.96 * cluster_prop_pix, DID4$coefficients[4] + 1.96 * cluster_prop_pix)*1
    covermat[i,9] <- between(ATT, DID4$coefficients[4] - 1.96 * cluster_grid_pix, DID4$coefficients[4] + 1.96 * cluster_grid_pix)*1
    covermat[i,10] <- between(ATT, DID4$coefficients[4] - 1.96 * cluster_county_pix, DID4$coefficients[4] + 1.96 * cluster_county_pix)*1
    
    ### coverage matrix for aggregated fixed effects
    
    fix_covermat[i,1] <- between(ATT, tail(fix_DID1$coefficients, n=1) - 1.96 * fix_se1, tail(fix_DID1$coefficients, n=1) + 1.96 * fix_se1)*1
    fix_covermat[i,2] <- between(ATT, tail(fix_DID2$coefficients, n=1) - 1.96 * fix_se2, tail(fix_DID2$coefficients, n=1) + 1.96 * fix_se2)*1
    fix_covermat[i,3] <- between(ATT, tail(fix_DID3$coefficients, n=1) - 1.96 * fix_se3, tail(fix_DID3$coefficients, n=1) + 1.96 * fix_se3)*1
    
    fix_covermat[i,4] <- between(ATT, tail(fix_DID1$coefficients, n=1) - 1.96 * fix_clse1, tail(fix_DID1$coefficients, n=1) + 1.96 * fix_clse1)*1
    fix_covermat[i,5] <- between(ATT, tail(fix_DID2$coefficients, n=1) - 1.96 * fix_clse2, tail(fix_DID2$coefficients, n=1) + 1.96 * fix_clse2)*1
    fix_covermat[i,6] <- between(ATT, tail(fix_DID3$coefficients, n=1) - 1.96 * fix_clse3, tail(fix_DID3$coefficients, n=1) + 1.96 * fix_clse3)*1
    
    
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
  
  
  coverages_df <- data.frame( 
    aggregation = c('grid', 
                   'property',
                   'county',
                   'grid',
                   'property',
                   'county',
                   'pixel',
                   'pixel',
                   'pixel',
                   'pixel'),
    std_error = c('clustered at grid', 
                 'clustered at property',
                 'clustered at county',
                 'classical',
                 'classical',
                 'classical',
                 'classical'
                 , 'clustered at property',
                 'clustered at grid'
                 , 'clustered at county'),
    coverage = colMeans(covermat)
  )
  
  fix_coverages_df <- data.frame( 
    fix_effects = c('grid', 
                    'property',
                    'county',
                    'grid',
                    'property',
                    'county'),
    std_error = c('classical',
                  'classical',
                  'classical',
                  'clustered at grid', 
                  'clustered at property',
                  'clustered at county'),
    coverage = colMeans(fix_covermat)
  )
    
  summary_df <- data.frame('DID'=rep(NA, 9), 'TWFE'=rep(NA, 9), 'pixels dropped' =rep(NA, 9),
                     'pixel'=rep(NA, 9),'grid'=rep(NA, 9),'property'=rep(NA, 9),'county'=rep(NA, 9),
                     'pixel fe'=rep(NA, 9),'grid fe'=rep(NA, 9),'property fe'=rep(NA, 9),'county fe'=rep(NA, 9), 
                     'mean_bias'=rep(NA, 9), 'RMSE'=rep(NA, 9),
                     'q05'=rep(NA, 9), 'q95'=rep(NA, 9),# as many cols as you need
                     'q01'=rep(NA, 9), 'q99'=rep(NA, 9),
                   stringsAsFactors=FALSE)
  
  # DID1 <- feols(deforrate ~  post*treat|year+grid, data = gridlevel_df)
  summary_df[4,] <- c(0,1,0,
            0,1,0,0,
            0,1,0,0,
            mean(coeff_bias$grid), rmse(actual, coeff_bias$grid),
            quantile(coeff_bias$grid, 0.05), quantile(coeff_bias$grid, 0.95),
            quantile(coeff_bias$grid, 0.01), quantile(coeff_bias$grid, 0.99))
  
  # DID2 <- feols(deforrate ~  post*treat|year+property, data = proplevel_df)
  summary_df[5,] <- c(0,1,0,
            0,0,1,0,
            0,0,1,0,
            mean(coeff_bias$property), rmse(actual, coeff_bias$property),
            quantile(coeff_bias$property, 0.05), quantile(coeff_bias$property, 0.95),
            quantile(coeff_bias$property, 0.01), quantile(coeff_bias$property, 0.99))
  
  # DID3 <- feols(deforrate ~  post*treat|year+county, data = countylevel_df)
  summary_df[6,] <- c(0,1,0,
            0,0,0,1,
            0,0,0,1,
            mean(coeff_bias$county), rmse(actual, coeff_bias$county),
            quantile(coeff_bias$county, 0.05), quantile(coeff_bias$county, 0.95),
            quantile(coeff_bias$county, 0.01), quantile(coeff_bias$county, 0.99))
  
  
  # DID4 <- feols(y_it ~  post*treat, data = panels)
  summary_df[3,] <- c(1,0,1,
            1,0,0,0,
            0,0,0,0,
            mean(coeff_bias$pixel), rmse(actual, coeff_bias$pixel),
            quantile(coeff_bias$pixel, 0.05), quantile(coeff_bias$pixel, 0.95),
            quantile(coeff_bias$pixel, 0.01), quantile(coeff_bias$pixel, 0.99))
  
  
  summary_df[1,] <- c(1,0,0,
            1,0,0,0,
            0,0,0,0,
            mean(bad_bias[,1]), rmse(actual, bad_bias[,1]),
            quantile(bad_bias[,1], 0.05), quantile(bad_bias[,1], 0.95),
            quantile(bad_bias[,1], 0.01), quantile(bad_bias[,1], 0.99))
  
  # DID6 <- feols(y_it ~  post*treat|year+pixels, data = panels)
  summary_df[2,] <- c(0,1,1,
                      1,0,0,0,
                      1,0,0,0,
                      mean(bad_bias[,2]), rmse(actual, bad_bias[,2]),
                      quantile(bad_bias[,2], 0.05), quantile(bad_bias[,2], 0.95),
                      quantile(bad_bias[,2], 0.01), quantile(bad_bias[,2], 0.99))
  
  ### TWFE regressions with aggregated fixed effects
  
  # fix_DID1 <- feols(y_it ~  post*treat|year+grid, data = panels)
  summary_df[7,] <- c(0,1,1,
            1,0,0,0,
            0,1,0,0,
            mean(fix_bias[,1]), rmse(actual, fix_bias[,1]),
            quantile(fix_bias[,1], 0.05), quantile(fix_bias[,1], 0.95),
            quantile(fix_bias[,1], 0.01), quantile(fix_bias[,1], 0.99))
  
  # fix_DID2 <- feols(y_it ~  post*treat|year+property, data = panels)
  summary_df[8,] <- c(0,1,1,
            1,0,0,0,
            0,0,1,0,
            mean(fix_bias[,2]), rmse(actual, fix_bias[,2]),
            quantile(fix_bias[,2], 0.05), quantile(fix_bias[,2], 0.95),
            quantile(fix_bias[,2], 0.01), quantile(fix_bias[,2], 0.99))
  
  # fix_DID3 <- feols(y_it ~  post*treat|year+county, data = panels)
  summary_df[9,] <- c(0,1,1,
            1,0,0,0,
            0,0,0,1,
            mean(fix_bias[,3]), rmse(actual, fix_bias[,3]),
            quantile(fix_bias[,3], 0.05), quantile(fix_bias[,3], 0.95),
            quantile(fix_bias[,3], 0.01), quantile(fix_bias[,3], 0.99))
  
  summary_df <- summary_df %>%
    mutate_at(1:10, as.logical)%>%
    dplyr::select(mean_bias, everything())
  
  
  outputs = list("plot" = plot, "fe_plot" = fe_plot, "biases" = coeff_bias, "fe_biases" = fix_bias, "coverages_df" = coverages_df, "fe_coverages_df" = fix_coverages_df, "summary_df" = summary_df)
  return(outputs)
  
  #end function  
}  
