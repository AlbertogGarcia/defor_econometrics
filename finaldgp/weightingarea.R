#this function uses the property_scape gen fcn to generate a landscape with properties. Then we introduce property level perturbations and compare estimates when the data is aggregateed to the grid vs. property level


library(ggplot2)
library(clubSandwich)
library(matrixStats)
library(ggplot2)
library(plm)
library(Metrics)
library(DataCombine)
library(dplyr)
library(tictoc)
source('county_scapegen.R')

#begin function
weightingarea <- function(n, nobs, years, b0, b1, b2, b3, std_a = 0.1, std_v = 0.25, std_p = .1, cellsize, ppoints, cpoints){
  
  countyscape = county_scapegen(nobs, cellsize, ppoints, cpoints)
  pixloc_df = countyscape$pixloc_df
  ATT <- pnorm(b0+b1+b2+b3, 0, (std_a^2+std_v^2 +std_p^2)^.5) - pnorm(b0+b1+b2, 0, (std_a^2+std_v^2 + std_p^2)^.5)
  DID_estimand <- (pnorm(b0+b1+b2+b3, 0, (std_a^2+std_v^2+std_p^2)^.5) - pnorm(b0+b1, 0, (std_a^2+std_v^2+std_p^2)^.5)
                   - (pnorm(b0+b2, 0, (std_a^2+std_v^2+std_p^2)^.5) - pnorm(b0, 0, (std_a^2+std_v^2+std_p^2)^.5)) )
  
  pixloc <- pixloc_df[order(pixloc_df$pixels),]
  
  covermat <- matrix(nrow = n, ncol = 10)

  coeffmatrix <- matrix(nrow = n, ncol = 6)
  
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
        ystar = b0 + b1*treat + b2*post + b3*treat*post + a_i + v_it
      )
    )
    
    #generate random 
    error_table <- data.frame(property = as.character(unique(pixloc$property)), p_err = rnorm(length(unique(pixloc$property)), 0, std_p))
    
    panels$pixels <- gsub("(?<![0-9])0+", "", panels$pixels, perl = TRUE)
    
    panels <- panels %>%
      inner_join(pixloc, by = c("pixels", "treat")) %>%
      inner_join(error_table, by = "property") %>%
      mutate(ystar = ystar + p_err) %>%
      mutate(y = (ystar > 0)*1 )
    
    #need to determine which year deforestation occurred
    year_df <- panels %>%
      select(pixels, year, y) %>%
      dcast(pixels ~ year , value.var = "y")
    
    rownames(year_df) <- year_df$pixels
    
    year_df <- year_df %>%
      select(- pixels)
    
    #creating variable for the year a pixel is deforested
    not_defor <- rowSums(year_df)<1 *1
    defor_year <- max.col(year_df, ties.method = "first") 
    defor_df <- transform(year_df, defor_year = ifelse(not_defor==1, years*2+1, defor_year))
    defor_df <- tibble::rownames_to_column(defor_df)
    names(defor_df)[1] <- paste("pixels")
    
    panels <- defor_df %>%
      select(pixels, defor_year) %>%
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
      proplevel_df <-  aggregate(panels, by = list(panels$property, panels$treat, panels$year), FUN = mean, drop = TRUE)[c("property", "treat", "post", "year","defor", "parea")]
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
      countylevel_df <-  aggregate(panels, by = list(panels$county, panels$treat, panels$year), FUN = mean, drop = TRUE)[c("county", "treat", "post", "year","defor", "carea")]
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
    
    
    
    
    
    
    DID1 <- plm(deforrate ~  post*treat, 
                data   = gridlevel_df, 
                #weights=(gridlevel_df$garea),
                method = "within", #fixed effects model
                effect = "twoway", #property and year fixed effects
                index  = c("grid", "year")
    )
    
    # run two-way fixed effects with outcome 1 
    DID2 <- plm(deforrate ~  post*treat, 
                data   = proplevel_df, 
                #weights=(proplevel_df$parea),
                method = "within", #fixed effects model
                effect = "twoway", #property and year fixed effects
                index  = c("property", "year")
    )
    
    DID3 <- plm(deforrate ~  post*treat, 
                data   = countylevel_df, 
                #weights=(countylevel_df$carea),
                method = "within", #fixed effects model
                effect = "twoway", #county and year fixed effects
                index  = c("county", "year")
    )
    
    DID4 <- plm(deforrate ~  post*treat, 
                data   = gridlevel_df, 
                weights=(gridlevel_df$garea),
                method = "within", #fixed effects model
                effect = "twoway", #property and year fixed effects
                index  = c("grid", "year")
    )
    
    DID7 <- lm(y_it ~  post*treat, 
                data   = panels) 
               
    
    # run two-way fixed effects with outcome 1 
    DID5 <- plm(deforrate ~  post*treat, 
                data   = proplevel_df, 
                weights=(proplevel_df$parea),
                method = "within", #fixed effects model
                effect = "twoway", #property and year fixed effects
                index  = c("property", "year")
    )
    
    DID6 <- plm(deforrate ~  post*treat, 
                data   = countylevel_df, 
                weights=(countylevel_df$carea),
                method = "within", #fixed effects model
                effect = "twoway", #county and year fixed effects
                index  = c("county", "year")
    )
    
    #calculating bias from each aggregation method
    coeffmatrix[i,1] <- DID1$coefficients - ATT
    coeffmatrix[i,2] <- DID2$coefficients - ATT
    coeffmatrix[i,3] <- DID3$coefficients - ATT
    
    coeffmatrix[i,4] <- DID4$coefficients - ATT
    coeffmatrix[i,5] <- DID5$coefficients - ATT
    coeffmatrix[i,6] <- DID6$coefficients - ATT
    
    #calculating standard errors and whether att is within CI
    cluster_se1    <- sqrt(diag(vcovHC(DID1, type = "HC0", cluster = "group")))
    cluster_se2    <- sqrt(diag(vcovHC(DID2, type = "HC0", cluster = "group")))
    cluster_se3    <- sqrt(diag(vcovHC(DID3, type = "HC0", cluster = "group")))
    covermat[i,1] <- between(ATT, DID1$coefficients - 1.96 * cluster_se1, DID1$coefficients + 1.96 * cluster_se1)*1
    covermat[i,2] <- between(ATT, DID2$coefficients - 1.96 * cluster_se2, DID2$coefficients + 1.96 * cluster_se2)*1
    covermat[i,3] <- between(ATT, DID3$coefficients - 1.96 * cluster_se3, DID3$coefficients + 1.96 * cluster_se3)*1
    
    se1 <- sqrt(DID1$vcov)
    se2 <- sqrt(DID2$vcov)
    se3 <- sqrt(DID3$vcov)
    se7 <- sqrt(diag(vcovHC(DID7)))[4]
    covermat[i,4] <- between(ATT, DID1$coefficients - 1.96 * se1, DID1$coefficients + 1.96 * se1)*1
    covermat[i,5] <- between(ATT, DID2$coefficients - 1.96 * se2, DID2$coefficients + 1.96 * se2)*1
    covermat[i,6] <- between(ATT, DID3$coefficients - 1.96 * se3, DID3$coefficients + 1.96 * se3)*1
    covermat[i,10] <- between(ATT, DID7$coefficients[4] - 1.96 * se7, DID7$coefficients[4] + 1.96 * se7)*1
    
    se4 <- sqrt(DID4$vcov)
    se5 <- sqrt(DID5$vcov)
    se6 <- sqrt(DID6$vcov)
    covermat[i,7] <- between(ATT, DID4$coefficients - 1.96 * se4, DID4$coefficients + 1.96 * se4)*1
    covermat[i,8] <- between(ATT, DID5$coefficients - 1.96 * se5, DID5$coefficients + 1.96 * se5)*1
    covermat[i,9] <- between(ATT, DID6$coefficients - 1.96 * se6, DID6$coefficients + 1.96 * se6)*1
    
   
    print(i)
    toc()
  }
  
  coeff_bias <- as.data.frame(coeffmatrix)
  actual <- rep(ATT, times = n)
  names(coeff_bias)[1] <- paste("ugrid")
  names(coeff_bias)[2] <- paste("uproperty")
  names(coeff_bias)[3] <- paste("ucounty")
  names(coeff_bias)[4] <- paste("wgrid")
  names(coeff_bias)[5] <- paste("wproperty")
  names(coeff_bias)[6] <- paste("wcounty")
  suppressWarnings(cbias <- melt(coeff_bias, value.name = "bias"))
  
  plot <- ggplot(data = cbias, aes(x = bias, fill=variable)) +
    geom_density(alpha = .2) +
    guides(fill=guide_legend(title=NULL))+
    scale_fill_discrete(breaks=c("ugrid", "uproperty", "ucounty", "wgrid", "wproperty", "wcounty"), labels=c("grid", "property", "county", "weighted grid", "weighted property", "weighted county"))+
    geom_vline(xintercept = 0, linetype = "dashed")+
    geom_vline(aes(xintercept= (DID_estimand - ATT), color="DID estimand - ATT"), linetype="dashed")+
    #theme(plot.margin = unit(c(1,1,3,1), "cm"))+
    theme(plot.caption = element_text(hjust = 0.5))+
    labs(x= "Bias", caption = paste("Mean grid:", round(mean(coeff_bias$ugrid), digits = 4),
                                    ", RMSE:", round(rmse(actual, coeff_bias$ugrid), digits = 4), "\n", 
                                    "Mean weighted grid:", round(mean(coeff_bias$wgrid), digits = 4),
                                    ", RMSE:", round(rmse(actual, coeff_bias$wgrid), digits = 4), "\n", 
                                    "Mean property:", round(mean(coeff_bias$uproperty), digits = 4),
                                    ", RMSE:", round(rmse(actual, coeff_bias$uproperty), digits = 4),"\n",
                                    "Mean weighted property:", round(mean(coeff_bias$wproperty), digits = 4),
                                    ", RMSE:", round(rmse(actual, coeff_bias$wproperty), digits = 4),"\n",
                                    "Mean county:", round(mean(coeff_bias$ucounty), digits = 4),
                                    ", RMSE:", round(rmse(actual, coeff_bias$ucounty), digits = 4),"\n", 
                                    "Mean weighted county:", round(mean(coeff_bias$wcounty), digits = 4),
                                    ", RMSE:", round(rmse(actual, coeff_bias$wcounty), digits = 4)
                                    ) 
    )
  
  grid_clustercover <- mean(covermat[,1]) 
  prop_clustercover <- mean(covermat[,2])
  county_clustercover <- mean(covermat[,3])  
  
  grid_cover <- mean(covermat[,4]) 
  prop_cover <- mean(covermat[,5])
  county_cover <- mean(covermat[,6]) 
  
  wgrid_cover <- mean(covermat[,7]) 
  wprop_cover <- mean(covermat[,8])
  wcounty_cover <- mean(covermat[,9])
  
  pix_cover <- mean(covermat[,10]) 
 
  outputs = list("plot" = plot, "biases" = coeff_bias, 
                 "grid_cover" = grid_cover, "prop_cover" = prop_cover, "county_cover" = county_cover, 
                 "grid_clustercover" = grid_clustercover, "prop_clustercover" = prop_clustercover, "county_clustercover" = county_clustercover, 
                 "wgrid_cover" = wgrid_cover, "wprop_cover" = wprop_cover, "wcounty_cover" = wcounty_cover
                 )
  
  return(outputs)
  
  #end function  
}  
