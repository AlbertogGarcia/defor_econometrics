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
xy_proppert <- function(n, nobs, years, ATT, base_0, base_1, trend, std_a = 0.1, std_v = 0.25, min_std_p = 0, max_std_p = .5, interval, cellsize, ppoints, cpoints){
  
  
  p_vals <- seq(from = min_std_p, to = max_std_p, by = interval)
  list <- as.list(p_vals)
  
  #preallocate n x ___ matrix
  gridmatrix <- matrix(nrow = n, ncol = length(list))
  countymatrix <- matrix(nrow = n, ncol = length(list))
  propertymatrix <- matrix(nrow = n, ncol = length(list))
  pixelmatrix <- matrix(nrow = n, ncol = length(list))
  
  clustergrid <- matrix(nrow = n, ncol = length(list))
  clusterprop <- matrix(nrow = n, ncol = length(list))
  clustercounty <- matrix(nrow = n, ncol = length(list))
  
  covergrid <- matrix(nrow = n, ncol = length(list))
  coverprop <- matrix(nrow = n, ncol = length(list))
  covercounty <- matrix(nrow = n, ncol = length(list))
  
  
  for(k in list){
    tic("loop")
    std_p <- k
    
    b0 = qnorm(base_0, mean = 0, sd = std_avp)
    b1 = qnorm(base_1, mean = 0, sd = std_avp) - b0
    b2 = qnorm(trend + base_0, mean = 0, sd = std_avp) - b0
    b3 = qnorm( pnorm(b0+b1+b2, mean = 0, sd = std_avp) + ATT , mean = 0, sd = std_avp) - (b0 + b1 + b2)
    
    DID_estimand <- (pnorm(b0+b1+b2+b3, 0, (std_a^2+std_v^2+std_p^2)^.5) - pnorm(b0+b1, 0, (std_a^2+std_v^2+std_p^2)^.5)
                     - (pnorm(b0+b2, 0, (std_a^2+std_v^2+std_p^2)^.5) - pnorm(b0, 0, (std_a^2+std_v^2+std_p^2)^.5)) )
    
    countyscape = county_scapegen(nobs, cellsize, ppoints, cpoints)
    pixloc_df = countyscape$pixloc_df
    
    pixloc <- pixloc_df[order(pixloc_df$pixels),]
    
    
    
    for(i in 1:n){
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
        mutate(defor = ifelse(indic > 0, 1, y)) %>%
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
      
      DID4 <- lm(y_it ~  post*treat, 
                  data   = panels
      )
      
      place <- which(list == k)
      
      #calculating bias from each aggregation method
      gridmatrix[i,place] <- DID1$coefficients - DID_estimand
      countymatrix[i,place] <- DID3$coefficients - DID_estimand
      propertymatrix[i,place] <- DID2$coefficients - DID_estimand
      pixelmatrix[i,place] <- DID4$coefficients[4] - DID_estimand
      
      cluster_se1    <- sqrt(diag(vcovHC(DID1, type = "HC0", cluster = "group")))
      cluster_se2    <- sqrt(diag(vcovHC(DID2, type = "HC0", cluster = "group")))
      cluster_se3    <- sqrt(diag(vcovHC(DID3, type = "HC0", cluster = "group")))
      clustergrid[i,place] <- between(DID_estimand, DID1$coefficients - 1.96 * cluster_se1, DID1$coefficients + 1.96 * cluster_se1)*1
      clusterprop[i,place] <- between(DID_estimand, DID2$coefficients - 1.96 * cluster_se2, DID2$coefficients + 1.96 * cluster_se2)*1
      clustercounty[i,place] <- between(DID_estimand, DID3$coefficients - 1.96 * cluster_se3, DID3$coefficients + 1.96 * cluster_se3)*1
      
      se1 <- sqrt(DID1$vcov)
      se2 <- sqrt(DID2$vcov)
      se3 <- sqrt(DID3$vcov)
      covergrid[i,place] <- between(DID_estimand, DID1$coefficients - 1.96 * se1, DID1$coefficients + 1.96 * se1)*1
      coverprop[i,place] <- between(DID_estimand, DID2$coefficients - 1.96 * se2, DID2$coefficients + 1.96 * se2)*1
      covercounty[i,place] <- between(DID_estimand, DID3$coefficients - 1.96 * se3, DID3$coefficients + 1.96 * se3)*1
      
    }
    print(k)
    toc()
  }
  
  gridbias_df <- as.data.frame(cbind(colMeans(gridmatrix), p_vals)) 
  countybias_df <- as.data.frame(cbind(colMeans(countymatrix), p_vals)) 
  propbias_df <- as.data.frame(cbind(colMeans(propertymatrix), p_vals)) 
  pixbias_df <- as.data.frame(cbind(colMeans(pixelmatrix), p_vals)) 
  
  gridcluster_df <- as.data.frame(cbind(colMeans(clustergrid), p_vals)) 
  countycluster_df <- as.data.frame(cbind(colMeans(clustercounty), p_vals)) 
  propcluster_df <- as.data.frame(cbind(colMeans(clusterprop), p_vals)) 
  
  covergrid_df <- as.data.frame(cbind(colMeans(covergrid), p_vals)) 
  covercounty_df <- as.data.frame(cbind(colMeans(covercounty), p_vals)) 
  coverprop_df <- as.data.frame(cbind(colMeans(coverprop), p_vals)) 
  
  plot <- ggplot() + 
    geom_line(data = gridbias_df, aes(x = p_vals, y = V1, color="aggregated to grid"), size =1) +
    geom_line(data = countybias_df, aes(x = p_vals, y = V1, color="aggregated to county"), size =1) +
    geom_line(data = propbias_df, aes(x = p_vals, y = V1, color="aggregated to property"), size =1) +
    geom_line(data = pixbias_df, aes(x = p_vals, y = V1, color="pixel"), size =1) +
    scale_color_manual(values = c(
      'aggregated to grid' = 'green',
      'aggregated to county' = 'red',
      'aggregated to property' = 'blue',
      'pixel' = 'yellow'))+
    labs(x = "property level std. error", y = "Bias", caption = paste("Bias as a function of property perturbations depends on aggregation method")) + 
    geom_hline(yintercept = 0, linetype = "dashed")
  
  plot2 <- ggplot() + 
    geom_line(data = gridcluster_df, aes(x = p_vals, y = V1, color = "grid clustered"), size =1, linetype = "dashed") +
    geom_line(data = countycluster_df, aes(x = p_vals, y = V1, color = "county clustered"), size =1, linetype = "dashed") +
    geom_line(data = propcluster_df, aes(x = p_vals, y = V1, color = "property clustered"), size =1, linetype = "dashed") +
    geom_line(data = covergrid_df, aes(x = p_vals, y = V1, color = "grid classical"), size =1) +
    geom_line(data = covercounty_df, aes(x = p_vals, y = V1, color = "county classical"), size =1) +
    geom_line(data = coverprop_df, aes(x = p_vals, y = V1, color = "property classical"), size =1) +
    scale_color_manual(values = c(
      'grid clustered' = 'green',
      'county clustered' = 'red',
      'property clustered' = 'blue',
      'grid classical' = 'darkgreen',
      'county classical' = 'darkred',
      'property classical' = 'darkblue'))+
    labs(x = "property level std. error", y = "coverage probability", caption = paste("Bias as a function of property perturbations depends on aggregation method")) + 
    geom_hline(yintercept = 0.95, linetype = "dashed")
  
  
  outputs = list("plot" = plot, "plot2" = plot2)
  
  return(outputs)
  
}



